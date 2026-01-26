#' Query available Aleutian Islands stations
#'
#' @description
#' This script is a translation of the PL/SQL procedure stored on the AI schema
#' in Oracle. It utilizes the functionality of creating/dropping tables in
#' Oracle using sqlSave in place of the Global Temporary table strategy
#' utilized in the PL/SQL procedure. This function also updates the selection
#' of AI candidate stations to incorporate designations of trawlability stored
#' in the AI.AIGRID_GIS table and originating with the electronic station logs
#' kept aboard the vessel by the FPC.
#'
#' Date last modified: 06/19/2019
#'
#' @author Ned Laman \email{ned.laman@@noaa.gov}
#'
#' @param channel open connection object. Created from gapindex::get_connected()
#'
#' @import \pkg{RODBC}
#'
#' @return A dataframe of candidate stations (field stationid) for Aleutian
#'  Islands bottom trawl surveys which includes the station centroid (fields
#'  longitude & latitude in decimal degrees) and stratum id.
#'

get.ai.stations <- function(channel = NULL) {

  ## Connect to Oracle if not already
  if(is.na(x = channel)) channel <- gapindex::get_connected()

  dat <- RODBC::sqlQuery(
    channel = channel,
    query = read_sql(filepath =
                       system.file("SQL/get_ai_stations.sql",
                                   package = "StationAllocationAIGOA")))

  return(dat)

}
