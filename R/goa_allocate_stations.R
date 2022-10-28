#' Gulf of Alaska Bottom Trawl Survey: Allocate stations across strata
#'
#' @description
#' To be filled in
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
#' @param vesssel_names character vector.
#' @param output_dir character string. Path for outputs.
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
           max_iter = 5000,
           year = 2023,
           vessel_names = c("vessel_1", "vessel_2"),
           output_dir = NULL)
  {
# return(StationAllocationAIGOA::goa_grid_2023)
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

    if (!is.null(output_dir))
      if (!dir.exists(output_dir))
        stop(paste("message from StationAllocationAIGOA::goa_allocate_stations:",
                   "directory path `output_dir` does not exist, please provide",
                   "a path for output products"))

    # ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ##   Constants
    # ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # spp_idx <- match(species, attributes(frame)$sp )
    # ns_opt <- length(spp_idx)
    # n_cells <- nrow(frame)
    # n_years <- unique(frame$WEIGHT)[1]
    #
    # strata_names <- levels(factor(grid_goa_sp@data$STRATUM))
    #
    #
    #
    # ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ##   Subset frame dataset to the specified species set
    # ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # frame <- frame[, c('domainvalue', "id", "WEIGHT",
    #                    paste0("Y", spp_idx), paste0("Y", spp_idx, "_SQ_SUM") )]
    # names(frame)[-1:-3] <- c(paste0("Y", 1:ns_opt),
    #                          paste0("Y", 1:ns_opt, "_SQ_SUM"))
    #
    # ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ##   CVs under SRS: upper CV bounds for MS allocation
    # ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ## Initiate CVs to be those calculated under SRS
    # ## buildStrataDF calculates the stratum means and variances, X1 = 1
    # ##     means to calculate those statics on the whole domain
    # frame$X1 <- 1
    # srs_stats <- SamplingStrata::buildStrataDF(dataset = frame)
    # srs_n <- n
    # srs_var <- as.matrix(srs_stats[, paste0("S", 1:ns_opt)])^2
    #
    # srs_var <- sweep(x = srs_var,
    #                  MARGIN = 1,
    #                  STATS = (1 - srs_n / n_cells) / srs_n,
    #                  FUN = "*")
    #
    # srs_cv <- sqrt(srs_var) / srs_stats[, paste0("M", 1:ns_opt)]
    #
    # ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ##   Single species CV: lower CV bounds for MS allocation
    # ##   Set X1 to strata
    # ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # frame$X1 <- as.integer(factor(grid_goa_sp@data$STRATUM))
    # strs_stats <- SamplingStrata::buildStrataDF(dataset = frame)
    # strs_stats <- strs_stats[order(as.integer(strs_stats$STRATO)), ]
    # strs_stats$N <- strs_stats$N / n_years
    #
    # ss_cv <- data.frame()
    # ss_allocations <- matrix(nrow = nrow(strs_stats), ncol = ns_opt,
    #                          dimnames = list(strata_names, species))
    #
    # message("")
    #
    # for (ispp in 1:ns_opt) {
    #   error_df <- data.frame("DOM" = "DOM1",
    #                          "CV1" = as.numeric(srs_cv[ispp]),
    #                          "domainvalue"  = 1)
    #   temp_stratif <- cbind(DOM1 = 1,
    #                         CENS = 0,
    #                         strs_stats[, c(names(strs_stats)[1:2],
    #                                        paste0(c("M", "S"), ispp)) ])
    #   names(temp_stratif)[-c(1:4)] <- paste0(c("M", "S"), 1)
    #
    #   temp_bethel <- SamplingStrata::bethel(
    #     errors = error_df,
    #     stratif = temp_stratif,
    #     realAllocation = T,
    #     printa = T,
    #     minnumstrat = 4)
    #   temp_n <- as.integer(sum(temp_bethel))
    #   temp_cv <- as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "])
    #
    #   iter = 1
    #   while (temp_n != n & iter != max_iter){
    #     over_under <- temp_n > n
    #     CV_adj <- ifelse(over_under == TRUE,
    #                      yes = 1.01,
    #                      no = 0.99)
    #
    #     error_df$CV1 <- temp_cv* CV_adj
    #     temp_bethel <- SamplingStrata::bethel(stratif = temp_stratif,
    #                                           errors = error_df,
    #                                           printa = TRUE,
    #                                           minnumstrat = min_n_per_stratum)
    #
    #     temp_n <- sum(as.numeric(temp_bethel))
    #     temp_cv <- as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "])
    #
    #     iter = 1 + iter
    #   }
    #
    #   ss_cv <- rbind(ss_cv, data.frame(species = species[ispp],
    #                                    ss_cv = temp_cv))
    #   ss_allocations[, ispp] <- temp_bethel
    #   message(paste0("Incorporating ", species[ispp]) )
    # }
    #
    # ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ##   multi species CV:
    # ##   Start at SRS CV --> bethel algorithm --> optimal sampling size (n_srs)
    # ##   If n_srs < n stations, reduce CV constraints across species by 0.1%
    # ##   using the CVs optimized for each species (ss_cv) as a lower bound.
    # ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # message("Now running multispecies allocation" )
    #
    # error_df <- cbind("DOM" = "DOM1",
    #                   srs_cv,
    #                   "domainvalue"  = 1)
    # names(error_df)[2:(1 + ns_opt)] <- paste0("CV", 1:ns_opt)
    #
    # temp_stratif <- cbind(DOM1 = 1,
    #                       CENS = 0,
    #                       strs_stats[, c(names(strs_stats)[1:2],
    #                                      paste0("M", 1:ns_opt),
    #                                      paste0("S", 1:ns_opt)) ])
    #
    # temp_bethel <- SamplingStrata::bethel(
    #   errors = error_df,
    #   stratif = temp_stratif,
    #   realAllocation = T,
    #   printa = T,
    #   minnumstrat = min_n_per_stratum)
    # temp_n <- sum(ceiling(temp_bethel))
    #
    # iter = 1
    # while (temp_n != n & iter != max_iter){
    #   over_under <- temp_n > n
    #   CV_adj <- ifelse(over_under == TRUE,
    #                    yes = 1.001,
    #                    no = 0.999)
    #
    #   updated_cv_constraint <-
    #     as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "]) * (CV_adj) +
    #     ss_cv$ss_cv * (1  - CV_adj)
    #
    #   error_df[, paste0("CV", 1:ns_opt)] <- as.numeric(updated_cv_constraint)
    #
    #   temp_bethel <- SamplingStrata::bethel(stratif = temp_stratif,
    #                                         errors = error_df,
    #                                         printa = TRUE,
    #                                         minnumstrat = min_n_per_stratum)
    #
    #   temp_n <- sum(as.numeric(temp_bethel))
    #   iter <- iter + 1
    # }
    #
    # ## Multispcies allocation and cvs
    # ms_cv <- as.numeric(updated_cv_constraint)
    # ms_allocation <- as.integer(ceiling(temp_bethel))
    # names(ms_allocation) <- levels(factor(grid_goa_sp@data$STRATUM))
    #
    # # ## Randomly drawn stations
    # drawn_stations <- c()
    #
    # for (i in 1:length(ms_allocation)) {
    #   set.seed(year)
    #   istratum <- names(ms_allocation)[i]
    #   available_stations <- with(stations@data,
    #                              which(STRATUM == istratum & TRAWL == T))
    #   temp_samples <- sample(x = available_stations,
    #                          size = ms_allocation[i],
    #                          prob = stations$AREA_KM2[available_stations],
    #                          replace = FALSE)
    #   drawn_stations <- c(drawn_stations, temp_samples)
    # }
    #
    # drawn_stations <- stations[drawn_stations, ]
    # drawn_stations$vessel <- vessel_names
    #
    # ## Expected CVs
    # expected_cvs <- cbind( t(srs_cv), ss_cv$ss_cv, ms_cv )
    # dimnames(expected_cvs) <- list(species, c("srs_cv", "ss_cv", "ms_cv"))
    #
    #
    # if (!is.null(output_dir)) {
    #   pdf(file = paste0(output_dir, "/goa_station_allocation_", year, ".pdf"),
    #       width = 10, height = 7, onefile = TRUE)
    #
    #   for (iarea in c("Shumagin", "Chirikof", "Kodiak", "Yakutat", "Southeast")) {
    #
    #     temp_strata <-
    #       sort(unique(grid_goa_sp@data$STRATUM[grid_goa_sp$MGT_AREA == iarea]))
    #
    #     par(mfrow = switch(paste(length(temp_strata)),
    #                        "5" = c(2, 3),
    #                        "4" =  c(2, 2)),
    #         mar = c(0,0,0,0))
    #
    #     for (istratum in temp_strata) {
    #       plot(subset(strata_list, MGT_AREA == iarea), border = F)
    #       plot(ak_land, col = "tan", add = T, border = "tan")
    #       plot(subset(x = strata_list, subset = STRATUM == istratum), add = TRUE,
    #            col = "red", border = F)
    #       plot(subset(x = stations, subset = STRATUM == istratum & TRAWL == F),
    #            add = TRUE, col = "grey", border = "grey")
    #
    #       for (ivessel in 1:length(vessel_names)) {
    #         plot(subset(x = drawn_stations,
    #                     subset = STRATUM == istratum &
    #                       vessel == vessel_names[ivessel]),
    #              col = c("black", "blue")[ivessel],
    #              border = c("black", "blue")[ivessel],
    #              add = TRUE)
    #       }
    #
    #       vessel_n <-
    #         table(drawn_stations$vessel[drawn_stations$STRATUM == istratum])
    #
    #       box()
    #       stratum_description <-
    #         with(subset(depth_mods, stratum == istratum),
    #              paste0(lower_depth_m, " - ", upper_depth_m, " m"))
    #
    #       legend(switch(iarea,
    #                     "Southeast" = "bottomleft",
    #                     "Yakutat" = "bottom",
    #                     "Kodiak" = "right",
    #                     "Chirikof" = "bottomright",
    #                     "Shumagin"= "bottom"),
    #              legend = paste0(vessel_names, ": ", vessel_n, " stations"),
    #              fill = c("black", "blue"),
    #              title = paste0(istratum, ": ", stratum_description))
    #     }
    #
    #   }
    #   dev.off()
    # }
    #
    # return(list(ss_allocations = ss_allocations,
    #             ms_allocation = ms_allocation,
    #             expected_cvs = as.data.frame(expected_cvs),
    #             drawn_stations = drawn_stations))
  }
