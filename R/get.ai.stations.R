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
#' @param channel open connection object. Created from RODBC::odbcConnect()
#' @param schema character string. SQL schema name. Defaults to AIGOA_WORK_DATA
#' @param pwrd character string. SQL password
#' @param survey character string. Defaults to "AI" for Aleutian Islands
#'
#' @import \pkg{RODBC}
#'
#' @return A dataframe of candidate stations (field stationid) for Aleutian
#'  Islands bottom trawl surveys which includes the station centroid (fields
#'  longitude & latitude in decimal degrees) and stratum id.
#'
#'


get.ai.stations <- function(channel = NA,
                            schema = "AIGOA_WORK_DATA",
                            pwrd = NA,
                            survey = 'AI'){

  if(is.na(channel)){
    require(RODBC)
    channel <- RODBC::odbcConnect(dsn = "AFSC", uid = schema,
                           pwd = pwrd, believeNRows = FALSE)
    close.channel = TRUE
  }else{
    close.channel <- FALSE
  }

  # Query to identify candidate AI stations
  stn.qry <- paste0("(select distinct stationid, stratum from racebase.haul
		where region = '", survey, "' and cruise >= 199401 and performance >= 0
		and stationid is not null and stratum is not null and (stationid, stratum)
		in (select distinct stationid,to_number(stratum) stratum
		from ", survey, ".", survey, "_gridpoints) group by stationid, stratum
		union select stationid, stratum from ", survey, ".", survey, "grid_gis
		where (stratum, stationid) in (select distinct  STRATUM, STATIONID
		from RACEBASE.HAUL where REGION = '", survey, "' and PERFORMANCE >= 0 and CRUISE = 199101
		and HAUL_TYPE = 3 and STATIONID is not null and STRATUM is not null
		and VESSEL = 85 minus select distinct STRATUM, STATIONID from RACEBASE.HAUL
		where PERFORMANCE >= 0 and HAUL_TYPE = 3 and STRATUM is not null and
		STATIONID is not null and CRUISE > 199400 and REGION = '", survey, "')
		and trawlable = 'Y') minus select distinct stationid, stratum
		from ", survey, ".", survey, "grid_gis where trawlable = 'N'")
  cat("\nPopulating candidate station list...\n")
  # run query
  TEMPSTNS <- sqlQuery(channel = channel, query = stn.qry)
  ## Create a new TEMPSTNS table with sqlSave
  sqlSave(channel = channel, dat = TEMPSTNS, rownames = F, fast = F, nastring = NULL)

  # Query to get area and position of all station segments for all candidate stations
  grid.qry <- paste0("select a.stationid, a.stratum, to_number(b.longitude) longitude,
		to_number(b.latitude) latitude, to_number(b.area) area from tempstns a,
		", survey, ".", survey, "_gridpoints b where a.stationid = b.stationid and a.stratum =
		to_number(b.stratum)")
  # run query
  GRDPTS <- sqlQuery(channel = channel, query = grid.qry)
  # Create a new GRDPTS table in Oracle with sqlSave
  sqlSave(channel = channel, dat = GRDPTS, rownames = F, fast = F, nastring = NULL)

  # Query to get largest continuous area within a station for assigning the station centroid
  maxarea.qry <- "select stationid, stratum, max(area) max_area from grdpts group by stationid, stratum"
  # run query
  MAXAREA <- sqlQuery(channel = channel, query = maxarea.qry)
  # Create a new MAXAREA table in Oracle with sqlSave
  sqlSave(channel = channel, dat = MAXAREA, rownames = F, fast = F, nastring = NULL)

  # assign station centroids to stations identified by stationid & stratum
  centroid.qry <- "select a.stationid, a.stratum, longitude, latitude
		from grdpts a, maxarea b where a.stationid = b.stationid and
		a.stratum = b.stratum and nvl(area,-1) = nvl(max_area, -1)"

  cat("\nAssigning coordinates for station centroids...\n")
  # run query
  dat <- sqlQuery(channel = channel, query = centroid.qry)
  names(dat) <- casefold(names(dat))

  # clean up by dropping temporary tables
  sqlDrop(channel, "TEMPSTNS")
  sqlDrop(channel, "GRDPTS")
  sqlDrop(channel, "MAXAREA")

  if(close.channel)close(channel)

  return(dat)

}
