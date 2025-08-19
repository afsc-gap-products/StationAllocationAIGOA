#' Simulate a simple random design Gulf of Alaska survey
#'
#' @author Zack Oyafuso \email{zack.oyafuso@@noaa.gov}
#' @description Internal function used in simulate_goa_strs()
#'
#' @param d_gc matrix with dimensions n_cells x n_spp where n_cells is the
#'             total number of cells and n_spp is the total number of species
#' @param sample_size integer, total number of random samples
#' @returns dataframe of mean, sd
#' @export
#'
sample_srs = function(d_gc = NULL,
                      sample_size=550){

  ss1 <- sample(1:nrow(x = d_gc), size = sample_size)
  species_names <- colnames(x = d_gc)

  #Calculate the mean of the random sample
  srs_mean <- apply(X = d_gc[ss1, ],
                    MARGIN = 2,
                    FUN = mean)
  srs_var <-apply(X = d_gc[ss1, ],
                  MARGIN = 2,
                  FUN = var)

  #Output: dataframe mean, sd, and CV of random sample for each species
  return(data.frame(species_names, srs_mean, srs_var))
}

#' Simulate a stratified random design Gulf of Alaska survey
#'
#' @author Zack Oyafuso \email{zack.oyafuso@@noaa.gov}
#'
#' @param d_gct array with dimensions n_cells x n_spp, n_years where n_cells
#'             is the total number of cells and n_spp is the total number of
#'             species, and n_years is the total number of survey years
#' @param grid_df data.frame with n_cells records, contains at least a column
#'                called "STRATUM" that has the stratum information
#' @param stn_allocation data.frame with n_strata records, where n_strata
#'                       is the total number of survey strata, contains
#'                       columns "stratum" and "n" for the station allocation.
#' @returns named list of stratum_info, stratified random means and cvs by
#'          species and year.
#' @export
#'
simulate_goa_strs <- function(d_gct = NULL,
                              grid_df = NULL,
                              stn_allocation = goa_stn_all) {

  n_strata <- nrow(x = stn_allocation)
  n_cells <- nrow(x = grid_df)
  n_species <- dim(x = d_gct)[2]
  n_years <- dim(x = d_gct)[3]
  species_names <- dimnames(x = d_gct)[[2]]

  ## Create a 20 x 6 matrix (WIS, 59x6 for TSRS) with 20 strata.
  ## The 6 cols: NNh,NN,nh,mean,sd,cv PvS 3/25/24
  stratum_info <-
    data.frame(stratum = names(x = table(grid_df$STRATUM)),
               NNh = as.numeric(table(grid_df$STRATUM)),
               WWh = as.numeric(table(grid_df$STRATUM) / n_cells)) |>
    merge(y = stn_allocation,
          by = "stratum")

  strata_names <- stratum_info$stratum

  ## finite-pop'n correction factor
  stratum_info$fin_corr <- with(stratum_info, (NNh - n) / NNh)

  mean_stats <- var_stats <-
    array(dim = c(n_strata, n_species, n_years),
          dimnames = list(strata_names, species_names, NULL))

  for (iyear in 1:n_years) {
    for (ii in strata_names) { #Loop over strata -- start

      stratum_idx <- grid_df$id[grid_df$STRATUM == ii]
      # plot(Lat ~ Lon, data = grid_df[grid_df$id %in% stratum_idx, ],
      # main = paste("Stratum", stratum_info$stratum[ii]),
      # pch = 16, cex = 0.5)

      #For each stratum ii, randomly sample stations based on the allocation
      #and populate "strata.stats" with the SRS mean, sd, CV
      stratum_stats <-
        sample_srs(d_gc = d_gct[stratum_idx,
                                ,
                                iyear],
                   sample_size = stratum_info$n[stratum_info$stratum == ii])
      mean_stats[ii, , iyear] <- stratum_stats$srs_mean
      var_stats[ii, , iyear] <- stratum_stats$srs_var
    }
  }

  #Calculate STRS mean, which is a weighted mean across all strata
  mean_all <- apply(X = mean_stats,
                    MARGIN = 2:3,
                    FUN = function(x) sum(x * stratum_info$WWh))

  #Calculate STRS standard error of the mean with finite-pop'n correction factor
  se_all <- apply(X = var_stats,
                  MARGIN = 2:3,
                  FUN = function(x)
                    sqrt(sum(stratum_info$WWh^2*stratum_info$fin_corr*(x)/stratum_info$n))
  )

  ## Calculate STRS CV
  cv_all <- round(x = se_all / mean_all, digits = 4)

  return(list(stratum_info = stratum_info,
              mean_all = mean_all,
              cv_all = cv_all))
}
