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
#' @param vessel_names character vector.
#'
#' @return A named list with elements:
#' 1) ms_allocation: vector, allocation of `n` stations across strata using
#'                   information from the `species` vector.
#' 2) drawn_stations: dataframe with `n` row and 15 fields. Notable fields are:
#' \describe{
#' \item{STATIONID}{Station ID of the grid the station resides. This format of the station id is XX-YY where XX is the grid location on the W-E axis and YY is the grid location on the N-S axis.}
#' \item{STRATUM}{Stratum ID. The 2023 restratification stratum id is formatted as the year the stratification was created (i.e., 2023), followed by "_" and then a three digit code in the format ABC similar to the format used in the historical GOA survey design. The second digit (B) refers to the managment area (1: Shumagin, 2: Chirikof, 3: Kodiak, 4: Yakutat, and 5: Southeastern). The first digit (A) denotes the relative depth of the stratum, 0 being the shallowest stratum and 5 being the deepest stratum. Because each management area has different depth strata, this digit is not comparable among management areas, just within a management area, e.g., 2023_010 is shallower than 2023_110. The third digit C is a replicate counter in case there are multiple strata within a management area x depth bin combination (there are no replicates in the 2023 GOA stratification, so this digit is always zero). This last digit is consistent with the historical GOA strata ids. }
#' \item{AREA_KM2}{area of station in km^2}
#' \item{TRAWLABLE_AREA_KM2}{trawlable area within station in km^2}
#' \item{PERIMETER_KM}{perimeter of station in km}
#' \item{CENTER_LAT, CENTER_LONG}{Latitude and longitude of station centroid}
#' \item{vessel}{vessel name assigned to station}
#' }
#'
#'@examples
#'
#' spp_list <- c("walleye pollock", "Pacific cod", "arrowtooth flounder",
#' "flathead sole", "rex sole", "northern rock sole",
#' "southern rock sole", "Dover sole", "Pacific halibut",
#' "Pacific ocean perch", "BS and RE rockfishes",
#' "silvergray rockfish", "dusky rockfish", "northern rockfish",
#' "shortspine thornyhead")
#'
#' GOA_allocation_2023 <-
#'   StationAllocationAIGOA::goa_allocate_stations(
#'     species = spp_list,
#'     n = 550,
#'     min_n_per_stratum = 4,
#'     year = 2023,
#'     vessel_names = c("vessel_1",  "vessel_2"))
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
           year = 2023,
           vessel_names = c("vessel_1", "vessel_2"))
  {

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##   Check that species list matches current species list
    ##   Check that year is an integer
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (!all(species %in% attributes(StationAllocationAIGOA::frame)$sp))
      stop(paste("message from StationAllocationAIGOA::goa_allocate_stations:",
                 "Argument `species` contains names that are not currently",
                 "in the list of included species. See",
                 "?StationAllocationAIGOA::goa_allocate_stations",
                 "for full species list"))

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##   Constants
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    spp_idx <- match(species,
                     attributes(StationAllocationAIGOA::frame)$sp )
    ns_opt <- length(spp_idx)
    n_cells <- nrow(StationAllocationAIGOA::frame)
    n_years <- unique(StationAllocationAIGOA::frame$WEIGHT)[1]

    strata_names <- levels(factor(StationAllocationAIGOA::frame$stratum))

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##   Subset frame dataset to the specified species set
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    frame_subdf <-
      StationAllocationAIGOA::frame[, c('domainvalue', "id", "WEIGHT",
                                        "stratum",
                                        paste0("Y", spp_idx),
                                        paste0("Y", spp_idx, "_SQ_SUM") )]
    names(frame_subdf)[-1:-4] <- c(paste0("Y", 1:ns_opt),
                                   paste0("Y", 1:ns_opt, "_SQ_SUM"))

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##   CVs under SRS: upper CV bounds for MS allocation
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## Initiate CVs to be those calculated under SRS
    ## buildStrataDF calculates the stratum means and variances, X1 = 1
    ##     means to calculate those statics on the whole domain
    frame_subdf$X1 <- 1
    srs_stats <- SamplingStrata::buildStrataDF(dataset = frame_subdf)
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
    frame_subdf$X1 <- as.integer(factor(frame_subdf$stratum))
    strs_stats <- SamplingStrata::buildStrataDF(dataset = frame_subdf)
    strs_stats <- strs_stats[order(as.integer(strs_stats$STRATO)), ]
    strs_stats$N <- strs_stats$N / n_years

    ss_cv <- data.frame()
    ss_allocations <- matrix(nrow = nrow(strs_stats), ncol = ns_opt,
                             dimnames = list(strata_names, species))

    message("")

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

      ## record cv and station allocation
      ss_cv <- rbind(ss_cv, data.frame(species = species[ispp],
                                       ss_cv = temp_cv))
      ss_allocations[, ispp] <- temp_bethel
      message(paste0("Incorporating ", species[ispp]) )
    } ## Loop over species -- start

    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ##   multi species CV:
    ##   Start at SRS CV --> bethel algorithm --> optimal sampling size (n_srs)
    ##   If n_srs < n stations, reduce CV constraints across species by 0.1%
    ##   using the CVs optimized for each species (ss_cv) as a lower bound.
    ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    message("Now running multispecies allocation" )

    ## set up bethel algorithm @ cv under srs
    error_df <- cbind("DOM" = "DOM1",
                      srs_cv,
                      "domainvalue"  = 1)
    names(error_df)[2:(1 + ns_opt)] <- paste0("CV", 1:ns_opt)

    temp_stratif <- cbind(DOM1 = 1,
                          CENS = 0,
                          strs_stats[, c(names(strs_stats)[1:2],
                                         paste0("M", 1:ns_opt),
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

    ## Record Multispcies allocation and cvs
    ms_cv <- as.numeric(updated_cv_constraint)
    ms_allocation <- as.integer(ceiling(temp_bethel))
    names(ms_allocation) <-
      levels(factor(StationAllocationAIGOA::frame$stratum))

    ## Randomly drawn stations
    drawn_stations <- c()

    for (i in 1:length(ms_allocation)) { ## Loop over strata -- start

      #Set seed
      set.seed(year)
      istratum <- names(ms_allocation)[i]

      ## available stations are those that are trawlable and > 5 km^2
      available_stations <- with(stations_2023,
                                 which(STRATUM == istratum &
                                         TRAWLABLE == T &
                                         AREA_KM2 >= 5.00))
      temp_samples <- sample(x = available_stations,
                             size = ms_allocation[i],
                             prob = stations_2023$AREA_KM2[available_stations],
                             replace = FALSE)
      drawn_stations <- c(drawn_stations, temp_samples)
    } ## Loop over strata -- end

    drawn_stations <- stations_2023[drawn_stations, ]
    drawn_stations$vessel <- vessel_names

    ## Expected CVs
    # expected_cvs <- cbind( t(srs_cv), ss_cv$ss_cv, ms_cv )
    # dimnames(expected_cvs) <- list(species, c("srs_cv", "ss_cv", "ms_cv"))

    return(list(#ss_allocations = ss_allocations,
                ms_allocation = ms_allocation,
                #expected_cvs = as.data.frame(expected_cvs),
                drawn_stations = drawn_stations))
  }
