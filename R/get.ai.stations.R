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

  TEMPSTNS <- RODBC::sqlQuery(
    channel = channel,
    query = read_sql(filepath =
                       system.file("SQL/get_ai_stations.sql",
                                   package = "StationAllocationAIGOA")))

  cat("\nPopulating candidate station list...\n")
  ## Create a new TEMPSTNS table with sqlSave
  RODBC::sqlSave(channel = channel,
                 dat = TEMPSTNS,
                 rownames = F,
                 fast = F,
                 nastring = NULL)

  # Query to get area and position of all station segments for all
  # candidate stations
  GRDPTS <- RODBC::sqlQuery(channel = channel,
                            query = "SELECT
                                     A.STATIONID,
                                     A.STRATUM,
                                     TO_NUMBER(B.CENTER_LONG) LONGITUDE,
	                                   TO_NUMBER(B.CENTER_LAT) LATITUDE,
	                                   TO_NUMBER(B.AREA_KM2) AREA
	                                   FROM
	                                   TEMPSTNS A,
	                                   AI.AIGRID_GIS B
	                                   WHERE A.STATIONID = B.STATIONID
	                                   AND A.STRATUM = TO_NUMBER(B.STRATUM)")

  # Create a new GRDPTS table in Oracle with sqlSave
  RODBC::sqlSave(channel = channel,
                 dat = GRDPTS,
                 rownames = F,
                 fast = F,
                 nastring = NULL)

  # Query to get largest continuous area within a station for assigning the
  # station centroid
  MAXAREA <- RODBC::sqlQuery(channel = channel,
                             query = "SELECT
                                      STATIONID,
                                      STRATUM,
                                      MAX(AREA) MAX_AREA
                                      FROM GRDPTS
                                      GROUP BY STATIONID, STRATUM")

  # Create a new MAXAREA table in Oracle with sqlSave
  RODBC::sqlSave(channel = channel,
                 dat = MAXAREA,
                 rownames = F,
                 fast = F,
                 nastring = NULL)

  cat("\nAssigning coordinates for station centroids...\n")
  # Assign station centroids to stations identified by stationid & stratum
  dat <- RODBC::sqlQuery(channel = channel,
                         query =  "SELECT
                                   A.STATIONID,
                                   A.STRATUM,
                                   LONGITUDE,
                                   LATITUDE
                                   FROM GRDPTS A,
                                   MAXAREA B
                                   WHERE A.STATIONID = B.STATIONID
                                   AND A.STRATUM = B.STRATUM
                                   AND NVL(AREA,-1) = NVL(MAX_AREA, -1)")
  names(dat) <- casefold(names(x = dat))

  # clean up by dropping temporary tables
  RODBC::sqlDrop(channel, "TEMPSTNS")
  RODBC::sqlDrop(channel, "GRDPTS")
  RODBC::sqlDrop(channel, "MAXAREA")

  return(dat)

}
