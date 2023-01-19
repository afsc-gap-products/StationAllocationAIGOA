#' Gulf of Alaska Bottom Trawl Survey: Allocate stations across strata
#'
#' @description
#' Allocates stations under the Gulf of Alaska stratified random design, last
#' restratified in 2023. Species information is incorporated from the survey
#' years between 1996-2021. Bethel algorithm is used to conduct the multispecies
#' station allocation.
#'
#' @author Zack Oyafuso \email{zack.oyafuso@@noaa.gov}
#'
#' @param n integer.Total number of stations
#' @param species character vector. See details for full list of species
#' @param min_n_per_stratum integer. Minimum number of stations to assign a
#'                          stratum
#' @param max_iter integer. Maximum number of times for bethel algorithm to run
#'                 to converge to `n`.
#' @param year integer. The year as an integer is used as a seed for the random
#'             number generator when drawing random stations.
#'
#' @return A named list with elements:
#' 1) ms_allocation: vector, allocation of `n` stations across strata using
#'                   information from the `species` vector.
#' 2) drawn_stations: dataframe with `n` row and 15 fields. Notable fields are:
#' \describe{
#' \item{STATIONID}{Station ID of the grid the station resides. This format of the station id is XX-YY where XX is the grid location on the W-E axis and YY is the grid location on the N-S axis.}
#' \item{STRATUM}{Stratum ID. The 2025 restratification stratum id is formatted as the year the stratification was created (i.e., 2025), followed by "_" and then a three digit code in the format ABC similar to the format used in the historical GOA survey design. The second digit (B) refers to the managment area (1: Shumagin, 2: Chirikof, 3: Kodiak, 4: Yakutat, and 5: Southeastern). The first digit (A) denotes the relative depth of the stratum, 0 being the shallowest stratum and 5 being the deepest stratum. Because each management area has different depth strata, this digit is not comparable among management areas, just within a management area, e.g., 2023_010 is shallower than 2023_110. The third digit C is a replicate counter in case there are multiple strata within a management area x depth bin combination (there are no replicates in the 2023 GOA stratification, so this digit is always zero). This last digit is consistent with the historical GOA strata ids. }
#' \item{AREA_KM2}{area of station in km^2}
#' \item{TRAWLABLE_AREA_KM2}{trawlable area within station in km^2}
#' \item{PERIMETER_KM}{perimeter of station in km}
#' \item{CENTER_LAT, CENTER_LONG}{Latitude and longitude of station centroid}
#' }
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
           max_iter = 5000,
           year = 2025){

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##   Check that species list matches current species list
    ##   Check that year is an integer
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (!all(species %in% dimnames(StationAllocationAIGOA::D_gct)[[2]]))
      stop(paste("message from StationAllocationAIGOA::goa_allocate_stations:",
                 "Argument `species` contains names that are not currently",
                 "in the list of included species. See",
                 "?StationAllocationAIGOA::goa_allocate_stations",
                 "for full species list"))

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##   Constants
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    dens <- StationAllocationAIGOA::D_gct
    spp_idx <- match(species, dimnames(dens)[[2]] )
    ns_opt <- length(spp_idx)
    n_cells <- dim(dens)[1]
    n_years <- dim(dens)[3]

    strata_names <- with(StationAllocationAIGOA::depth_mods_2023,
                         stratum[used])
    NMFS_area <- with(StationAllocationAIGOA::depth_mods_2023,
                      manage_area[used])
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##   Single species CV: lower CV bounds for MS allocation
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    srs_var <- apply(X = dens[, spp_idx, ],
                     MARGIN = 2,
                     FUN = function(x) var(as.vector(x)))
    srs_var <- sweep(x = as.matrix(srs_var),
                     MARGIN = 1,
                     STATS = (1 - n / n_cells) / n,
                     FUN = "*")
    srs_mean <- apply(X = dens[, spp_idx, ],
                      MARGIN = 2,
                      FUN = function(x) mean(as.vector(x)))
    srs_cv <- sqrt(srs_var) / srs_mean

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##   Single species CV: lower CV bounds for MS allocation
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    strs_sd <- apply(X = dens[, spp_idx, ],
                     MARGIN = 2,
                     FUN = function(x)
                       sapply(X = split(x = x,
                                        f = optim_df$STRATUM),
                              FUN = function(xx) sd(as.vector(xx))))
    strs_mean <- apply(X = dens[, spp_idx, ],
                       MARGIN = 2,
                       FUN = function(x)
                         sapply(X = split(x = x,
                                          f = optim_df$STRATUM),
                                FUN = function(xx) mean(as.vector(xx))))

    strs_stats <- cbind(data.frame(
      STRATO = 1:length(strata_names),
      N = as.numeric(table(optim_df$STRATUM)),
      COST = 1, CENS = 0, DOM1 = 1,
      X1 = 1:length(strata_names)#,
      # row.names = strata_names
    ),
    matrix(data = as.vector(strs_mean), ncol = ns_opt,
           dimnames = list(NULL, paste0("M", 1:ns_opt))),
    matrix(data = as.vector(strs_sd), ncol = ns_opt,
           dimnames = list(NULL, paste0("S", 1:ns_opt))) )

    ## Data objects for ss allocations
    ss_cv <- data.frame()
    ss_allocations <- matrix(nrow = length(strata_names),
                             ncol = ns_opt,
                             dimnames = list(strata_names, species))

    for (ispp in 1:ns_opt) { ## Loop over species -- start

      ## Set up bethel algorithm @ cv under srs
      error_df <- data.frame("DOM" = "DOM1",
                             "CV1" = as.numeric(srs_cv[ispp]),
                             "domainvalue"  = 1)
      temp_stratif <- cbind(DOM1 = 1,
                            CENS = 0,
                            strs_stats[, c(names(strs_stats)[1:2],
                                           paste0(c("M", "S"), ispp)) ])
      names(temp_stratif)[-c(1:4)] <- paste0(c("M", "S"), 1)

      # Run initial allocation and record the total effort and cv
      temp_bethel <- SamplingStrata::bethel(
        errors = error_df,
        stratif = temp_stratif,
        realAllocation = T,
        printa = T,
        minnumstrat = 4)
      temp_n <- as.integer(sum(temp_bethel))
      temp_cv <- as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "])

      ## iteratively adjust temp_cv and rerun bethel until the total effort = n
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

      ## If max_iter it hit and temp_n != n, shoot out a warning
      if (temp_n != n & iter != max_iter)
        warning(paste0("Maximum iteration has been reached for single-species ",
                       "allocation of ", species[ispp], ". Total stations ",
                       "allocated may not be equal to ", n,
                       ". Consider increasing the `max_iter` ",
                       "argument in goa_allocate_stations()."))

      ## record cv and station allocation
      ss_cv <- rbind(ss_cv, data.frame(species = species[ispp],
                                       ss_cv = temp_cv))
      ss_allocations[, ispp] <- temp_bethel
      message(paste0("Incorporating ", species[ispp]) )
    } ## Loop over species -- end

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##   multi species CV:
    ##   Start at SRS CV --> bethel algorithm --> optimal sampling size (n_srs)
    ##   If n_srs < n stations, reduce CV constraints across species by 0.1%
    ##   using the CVs optimized for each species (ss_cv) as a lower bound.
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    message("Now running multispecies allocation" )

    ## set up bethel algorithm @ cv under srs
    error_df <- cbind(data.frame("DOM" = "DOM1"),
                      matrix(srs_cv, nrow = 1,
                             dimnames = list(NULL, paste0("CV", 1:ns_opt))),
                      data.frame("domainvalue"  = 1))

    temp_stratif <- cbind(DOM1 = 1,
                          CENS = 0,
                          STRATO = strs_stats$STRATO,
                          N = strs_stats$N,
                          strs_stats[, c(paste0("M", 1:ns_opt),
                                         paste0("S", 1:ns_opt)) ])

    ## Run initial allocation and record the total effort and cv
    temp_bethel <- SamplingStrata::bethel(
      errors = error_df,
      stratif = temp_stratif,
      realAllocation = T,
      printa = T,
      minnumstrat = min_n_per_stratum)
    temp_n <- sum(ceiling(temp_bethel))

    ## iteratively adjust temp_cv and rerun bethel until the total effort = n
    ## Use the single-species cvs as a lower bound when adjusting temp_cv
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

    ## If max_iter it hit and temp_n != n, shoot out a warning
    if (temp_n != n & iter != max_iter)
      warning(paste0("Maximum iteration has been reached for ",
                     "MS allocation. Total stations allocated may not be ",
                     "equal to ", n, ". Consider increasing the `max_iter` ",
                     "argument in goa_allocate_stations()."))

    ## Record Multispcies allocation and cvs
    ms_cv <- as.numeric(updated_cv_constraint)
    ms_allocation <- as.integer(ceiling(temp_bethel))
    names(ms_allocation) <- rownames(strs_sd)
    ms_allocation <- ms_allocation[strata_names]

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## Randomly drawn stations
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    drawn_stations <- c()

    for (i in 1:length(ms_allocation)) { ## Loop over strata -- start

      #Set seed
      set.seed(year)
      istratum <- names(ms_allocation)[i]

      ## available stations are those that are trawlable and > 5 km^2
      available_stations <- with(stations_2023,
                                 which(STRATUM == istratum &
                                         TRAWLABLE != "N" &
                                         AREA_KM2 >= 5.00 |
                                         TRAWLABLE == "N" &
                                         TRAWLABLE_AREA_KM2 >= 5.00))
      temp_samples <- sample(x = available_stations,
                             size = ms_allocation[i],
                             prob = stations_2023$AREA_KM2[available_stations],
                             replace = FALSE)
      drawn_stations <- c(drawn_stations, temp_samples)
    } ## Loop over strata -- end

    drawn_stations <- stations_2023[drawn_stations, ]
    drawn_stations$vessel <- vessel_names

    return(list(ms_allocation = data.frame(stratum = strata_names,
                                           nmfs_area = NMFS_area,
                                           ms_allocation = ms_allocation),
                drawn_stations = drawn_stations))
  }
