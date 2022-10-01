#' Gulf of Alaska Bottom Trawl Survey: Allocate stations across strata
#'
#' @description
#' To be filled in
#'
#' @author Zack Oyafuso \email{zack.oyafuso@@noaa.gov}
#'
#' @param n integer.Total number of stations
#' @param species character vector. See details for full list of species
#' @param min_n_per_stratum integer. Minimum number of stations to assign a stratum
#' @param max_iter integer. Maximum number of times for bethel algorithm to run to converge to `n`.
#'
#' @return A named vector of stations across strata
#'

goa_allocate_stations <-
  function(n = 550,
           min_n_per_stratum = 4,
           species = c("arrowtooth flounder", ## Atherestes stomias
                       "Pacific cod", ## Gadus macrocephalus
                       "walleye pollock", ## Gadus chalcogrammus
                       "rex sole", ## Glyptocephalus zachirus
                       "flathead sole", ## Hippoglossoides elassodon
                       "Pacific halibut", ## Hippoglossus stenolepis
                       "southern rock sole", ## Lepidopsetta bilineata
                       "northern rock sole", ## Lepidopsetta polyxystra
                       "Pacific ocean perch", ## Sebastes alutus
                       "silvergray rockfish", ## Sebastes brevispinis
                       "northern rockfish", ## Sebastes polyspinis
                       "dusky rockfish", ## Sebastes variabilis
                       "BS and RE rockfishes", ## Sebastes aleutianus and S. melanostictus
                       "Dover sole", ## Microstomus pacificus
                       "shortspine thornyhead" ## Sebastolobus alascanus)
           ),
           max_iter = 5000)
  {

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##   Load data
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    data(frame)
    data(grid_goa_sp)

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##   Check that species list matches current species list
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(!all(species %in% attributes(frame)$sp))
      stop(paste("Argument `species` contains names that are not currently",
                 "in the list of included species. See",
                 "?AIGOASurveyPlanning::goa_allocate_stations for full list"))

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##   Constants
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    spp_idx <- match(species, attributes(frame)$sp )
    ns_opt <- length(spp_idx)
    n_cells <- nrow(frame)
    n_years <- unique(frame$WEIGHT)[1]

    strata_names <- levels(factor(grid_goa_sp@data$STRATUM))

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##   Subset frame dataset to the specified species set
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    frame <- frame[, c('domainvalue', "id", "WEIGHT",
                       paste0("Y", spp_idx), paste0("Y", spp_idx, "_SQ_SUM") )]
    names(frame)[-1:-3] <- c(paste0("Y", 1:ns_opt),
                             paste0("Y", 1:ns_opt, "_SQ_SUM"))

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##   CVs under SRS: upper CV bounds for MS allocation
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## Initiate CVs to be those calculated under SRS
    ## buildStrataDF calculates the stratum means and variances, X1 = 1
    ##     means to calculate those statics on the whole domain
    frame$X1 <- 1
    srs_stats <- SamplingStrata::buildStrataDF(dataset = frame)
    srs_n <- n
    srs_var <- as.matrix(srs_stats[, paste0("S", 1:ns_opt)])^2

    srs_var <- sweep(x = srs_var,
                     MARGIN = 1,
                     STATS = (1 - srs_n / n_cells) / srs_n,
                     FUN = "*")

    srs_cv <- sqrt(srs_var) / srs_stats[, paste0("M", 1:ns_opt)]

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##   Single species CV: lower CV bounds for MS allocation
    ##   Set X1 to strata
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    frame$X1 <- as.integer(factor(grid_goa_sp@data$STRATUM))
    strs_stats <- SamplingStrata::buildStrataDF(dataset = frame)
    strs_stats <- strs_stats[order(as.integer(strs_stats$STRATO)), ]
    strs_stats$N <- strs_stats$N / n_years

    ss_cv <- data.frame()
    ss_allocations <- matrix(nrow = nrow(strs_stats), ncol = ns_opt,
                             dimnames = list(strata_names, species))

    message("")

    for (ispp in 1:ns_opt) {
      error_df <- data.frame("DOM" = "DOM1",
                             "CV1" = as.numeric(srs_cv[ispp]),
                             "domainvalue"  = 1)
      temp_stratif <- cbind(DOM1 = 1,
                            CENS = 0,
                            strs_stats[, c(names(strs_stats)[1:2],
                                           paste0(c("M", "S"), ispp)) ])
      names(temp_stratif)[-c(1:4)] <- paste0(c("M", "S"), 1)

      temp_bethel <- SamplingStrata::bethel(
        errors = error_df,
        stratif = temp_stratif,
        realAllocation = T,
        printa = T,
        minnumstrat = 4)
      temp_n <- as.integer(sum(temp_bethel))
      temp_cv <- as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "])

      iter = 1
      while (temp_n != n & iter != max_iter){
        over_under <- temp_n > n
        CV_adj <- ifelse(over_under == TRUE,
                         yes = 1.01,
                         no = 0.99)

        error_df$CV1 <- temp_cv* CV_adj
        temp_bethel <- SamplingStrata::bethel(stratif = temp_stratif,
                                              errors = error_df,
                                              printa = TRUE,
                                              minnumstrat = min_n_per_stratum)

        temp_n <- sum(as.numeric(temp_bethel))
        temp_cv <- as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "])

        iter = 1 + iter
      }

      ss_cv <- rbind(ss_cv, data.frame(species = species[ispp],
                                       ss_cv = temp_cv))
      ss_allocations[, ispp] <- temp_bethel
      message(paste0("Incorporating ", species[ispp]) )
    }

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##   multi species CV:
    ##   Start at SRS CV --> bethel algorithm --> optimal sampling size (n_srs)
    ##   If n_srs < n stations, reduce CV constraints across species by 0.1%
    ##   using the CVs optimized for each species (ss_cv) as a lower bound.
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    message("Now running multispecies allocation" )

    error_df <- cbind("DOM" = "DOM1",
                      srs_cv,
                      "domainvalue"  = 1)
    names(error_df)[2:(1 + ns_opt)] <- paste0("CV", 1:ns_opt)

    temp_stratif <- cbind(DOM1 = 1,
                          CENS = 0,
                          strs_stats[, c(names(strs_stats)[1:2],
                                         paste0("M", 1:ns_opt),
                                         paste0("S", 1:ns_opt)) ])

    temp_bethel <- SamplingStrata::bethel(
      errors = error_df,
      stratif = temp_stratif,
      realAllocation = T,
      printa = T,
      minnumstrat = min_n_per_stratum)
    temp_n <- sum(ceiling(temp_bethel))

    iter = 1
    while (temp_n != n & iter != max_iter){
      over_under <- temp_n > n
      CV_adj <- ifelse(over_under == TRUE,
                       yes = 1.001,
                       no = 0.999)

      updated_cv_constraint <-
        as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "]) * (CV_adj) +
        ss_cv$ss_cv * (1  - CV_adj)

      error_df[, paste0("CV", 1:ns_opt)] <- as.numeric(updated_cv_constraint)

      temp_bethel <- SamplingStrata::bethel(stratif = temp_stratif,
                                            errors = error_df,
                                            printa = TRUE,
                                            minnumstrat = min_n_per_stratum)

      temp_n <- sum(as.numeric(temp_bethel))
      iter <- iter + 1
    }

    ms_cv <- as.numeric(updated_cv_constraint)
    ms_allocation <- as.integer(ceiling(temp_bethel))
    names(ms_allocation) <- levels(factor(grid_goa_sp@data$STRATUM))

    expected_cvs <- cbind( t(srs_cv), ss_cv$ss_cv, ms_cv )
    dimnames(expected_cvs) <- list(species, c("srs_cv", "ss_cv", "ms_cv"))

    return(list(ss_allocations = ss_allocations,
                ms_allocation = ms_allocation,
                expected_cvs = as.data.frame(expected_cvs)))
  }
