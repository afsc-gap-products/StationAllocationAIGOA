#' Randomly draw stations under a stratified random design
#'
#' @param station_allocation dataframe
#' @param survey_year integer. The year as an integer is used as a seed for the random
#'             number generator when drawing random stations.
#' @param trawl character vector. Options are "Y" for trawlable, "N" for
#'              untrawlable, or "UNK" for unknown trawlablity.
#' @return A dataframe with `n` row and 15 fields. Notable fields are:
#' \describe{
#' \item{STATIONID}{Station ID of the grid the station resides. This format of the station id is XX-YY where XX is the grid location on the W-E axis and YY is the grid location on the N-S axis.}
#' \item{STRATUM}{Stratum ID. The 2025 restratification stratum id is formatted as the year the stratification was created (i.e., 2025), followed by "_" and then a three digit code in the format ABC similar to the format used in the historical GOA survey design. The second digit (B) refers to the managment area (1: Shumagin, 2: Chirikof, 3: Kodiak, 4: Yakutat, and 5: Southeastern). The first digit (A) denotes the relative depth of the stratum, 0 being the shallowest stratum and 5 being the deepest stratum. Because each management area has different depth strata, this digit is not comparable among management areas, just within a management area, e.g., 2023_010 is shallower than 2023_110. The third digit C is a replicate counter in case there are multiple strata within a management area x depth bin combination (there are no replicates in the 2023 GOA stratification, so this digit is always zero). This last digit is consistent with the historical GOA strata ids. }
#' @export
#'

draw_strs_stations <- function(stn_allocation,
                               survey_year,
                               stn_df,
                               trawl = c("Y", "N", "UNK")[c(1, 3)],
                               min_area_km2 = 5) {
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Randomly drawn stations
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  drawn_stations <- c()
  set.seed(survey_year)

  for (i in 1:nrow(x = stn_allocation)) { ## Loop over strata -- start

    #Set seed

    istratum <- stn_allocation$stratum[i]

    ## available stations are those that are trawlable and > min_area_km2
    available_stations <- with(stn_df,
                               which(STRATUM == istratum &
                                       TRAWLABLE %in% trawl &
                                       AREA_KM2 >= min_area_km2))
    temp_samples <- sample(x = available_stations,
                           size = stn_allocation$ms_allocation[i],
                           prob = stn_df$AREA_KM2[available_stations],
                           replace = FALSE)
    drawn_stations <- c(drawn_stations, temp_samples)
  } ## Loop over strata -- end

  return(stn_df[drawn_stations, ])

}
