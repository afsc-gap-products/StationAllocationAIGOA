#' Allocate stations across strata
#'
#' @description
#' Allocates stations across strata based on: stratum size, historical standard
#' deviation of cpue, ex vessel price, and abundance of a preselected subset of
#' species (reference a tech memo for how exactly Neyman allocations are handled?).
#'
#' Created 12/18/13
#'
#' @author Ned Laman \email{ned.laman@@noaa.gov}
#'
#' @param channel open connection object. Created from gapindex::get_connected()
#' @param species_weightings dataframe of species codes and weightings
#' @param stratum_stats dataframe of stratum means and sds for each species
#' @param n_total integer. Number of total stations
#' @param min_n_stratum integer. Minimum number of stations in a stratum
#'
#' @return A named vector of allocated stations across strata
#'
#' @export
#'

allocate.effort <- function(channel = NULL,
                            species_weightings = NULL,
                            stratum_stats = NULL,
                            n_total = 400,
                            min_n_stratum = 2){

  ## Connect to Oracle if not already
  if (is.null(x = channel)) channel <- gapindex::get_connected()
  obs_years <- sort(x = unique(x = stratum_stats$YEAR))
  obs_species <- sort(x = unique(x = stratum_stats$SPECIES_CODE))

  ## Query stratum information from Oracle
  strata <- RODBC::sqlQuery(channel = channel,
                            query = "SELECT AREA_ID AS STRATUM,
                                   AREA_KM2 AS AREA
                                   FROM GAP_PRODUCTS.AREA
                                   WHERE SURVEY_DEFINITION_ID = 52
                                   AND DESIGN_YEAR = 1980
                                   AND AREA_TYPE = 'STRATUM'",
                            believeNRows = FALSE)

  ## Attach the stratum area data from GAP_PRODUCTS.AREA to the stratum_stats
  ## df using "STRATUM" as the key.
  stratum_stats <- merge(x = stratum_stats,
                         y = strata,
                         by = "STRATUM",
                         all = TRUE)

  ## Loop through each year and species and calculate the neyman allocation
  ## for that combination, append to n_opt.
  n_opt <- data.frame()

  for (iyear in obs_years) { ## Loop over years -- start
    for (ispp in obs_species) { ## Loop over species -- start

      ## subset stratum_stats to only iyear and ispp
      subdf = subset(x = stratum_stats,
                     subset = YEAR == iyear & SPECIES_CODE == ispp)

      ## For some species like rougheye, dusky, etc., their time series starts
      ## after 1991, so they won't have biomass values in the earlier years.
      if (sum(subdf$CPUE_KGKM2_MEAN) == 0) next

      ## For strata with only one station, there is no sd calculation,
      ## so we impute the variance with the mean sd across strata
      subdf$CPUE_KGKM2_SD[is.na(x = subdf$CPUE_KGKM2_SD)] <-
        mean(subdf$CPUE_KGKM2_SD, na.rm = T)

      ## Calculate the Neyman allocation with the caveat that each stratum
      ## will have at least the min_n_stratum in each stratum
      n_by_stratum <- n_total * (subdf$CPUE_KGKM2_SD * subdf$AREA) /
        max(1, sum(subdf$CPUE_KGKM2_SD * subdf$AREA))
      n_by_stratum <- sapply(X = n_by_stratum,
                             FUN = function(x)
                               ifelse(test = x >= min_n_stratum,
                                      yes = x,
                                      no = min_n_stratum))
      n_by_stratum <- round(x = n_by_stratum, digits = 0)
      temp_total <- sum(n_by_stratum)

      ## Usually when the minimum stratum effort constraint occurs, this bumps
      ## the total sample size to above the requested total sample size
      ## (e.g., print(temp_total)). At least for now, there is not an elegant
      ## way to include the minimum stratum effort as a formal constraint in
      ## the Neyman allocation for the AI BTS survey station allocation.
      ##
      ## This next iteration step gradually reduces the total sample size in the
      ## Neyman allocation calculation until the 'effective' sample size is
      ## equal to the requested sample size.
      initial_n <- n_total
      iter = 1

      while((temp_total != n_total) & (iter != 1000) ) {
        initial_n <- ifelse(temp_total < n_total,
                            initial_n + 1,
                            initial_n - 1)

        n_by_stratum <- initial_n * (subdf$CPUE_KGKM2_SD * subdf$AREA) /
          max(1, sum(subdf$CPUE_KGKM2_SD * subdf$AREA))
        n_by_stratum <- sapply(X = n_by_stratum,
                               FUN = function(x)
                                 ifelse(test = x >= min_n_stratum,
                                        yes = x,
                                        no = min_n_stratum))
        n_by_stratum <- round(x = n_by_stratum, digits = 0)

        temp_total <- sum(n_by_stratum)
        iter <- iter + 1
      }

      subdf$N_OPT <- n_by_stratum

      n_opt <- rbind(n_opt, subdf)
    } ## Loop over species -- end
  } ## Loop over years -- end

  ## Calculate each species' mean stratum effort allocation across years
  mean_allocation_spp <-
    do.call(what = rbind,
            args = lapply(X = split(x = n_opt,
                                    f = list(n_opt$STRATUM,
                                             n_opt$SPECIES_CODE)),
                          FUN = function(subdf) {
                            subdf$MEAN_N_OPT <- mean(subdf$N_OPT)
                            return(unique(subdf[, c("STRATUM",
                                                    "SPECIES_CODE",
                                                    "MEAN_N_OPT")]))
                          }))

  ## Merge the species weighting df to mean_allocation_spp using SPECIES_CODE
  ## as the key
  mean_allocation_spp <- merge(x = mean_allocation_spp,
                               y = species_weightings,
                               by = "SPECIES_CODE")

  ## For each stratum the effort allocated will be a weighted average across
  ## species, weighted by the mean economic value of the species.
  weighted_allocation <-
    do.call(what = c,
            args = lapply(X = split(x = mean_allocation_spp,
                                    f = mean_allocation_spp$STRATUM),
                          FUN = function(subdf) {
                            weighted.mean(x = subdf$MEAN_N_OPT,
                                          w = subdf$econ_value)
                          }))

  ## Round and return return
  return(round(x = weighted_allocation, digits = 0))
}
