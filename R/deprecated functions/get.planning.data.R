#' Stratum-level catch stats
#'
#' @description
#' Planning data were generated from a schema-dependent stored procedure
#' (STATION_PLANNING) in Oracle. Average and standard deviation of CPUE hauls
#' for a subset of preselected species were calculated for each year and
#' stratum id.
#' Created 06/19/2019
#'
#' @param channel open connection object. Created from RODBC::odbcConnect()
#' @param schema character string. SQL schema name. Defaults to AIGOA_WORK_DATA
#' @param pwrd character string. SQL password
#' @param survey character string. Defaults to "AI" for Aleutian Islands
#'
#' @author Ned Laman \email{ned.laman@@noaa.gov}
#'
#' @return A dataframe containing:\tabular{ll}{
#'    \code{year} \tab integer. Year of survey. \cr \tab \cr
#'    \code{stratum} \tab integer. Stratum id. \cr \tab \cr
#'    \code{species_code} \tab integer. Species code. \cr \tab \cr
#'    \code{haul_count} \tab integer. Total number of hauls. \cr \tab \cr
#'    \code{cpue} \tab numeric. Average CPUE. \cr \tab \cr
#'    \code{cpue_std} \tab numeric. Standard deviation of CPUE. \cr \tab \cr
#' }
#'
#' These values inform the allocation of stations (effort) to survey strata.
#'

get.planning.data <- function(channel = NA, survey = "AI",
                              schema = "AIGOA_WORK_DATA", pwrd = pwrd) {

  # test channel to Oracle and open if necessary
  if(is.na(channel)){
    require(RODBC)
    channel <- odbcConnect(dsn = "AFSC", uid = schema, pwd = pwrd, believeNRows = FALSE)
    close.channel = TRUE
  }else{
    close.channel <- FALSE
  }

  # SQLPLUS query translated from STATION_PLANNING procedure stored on AI schema in Oracle
  plandata.qry <- paste0("select Y.YEAR,Y.STRATUM,Y.SPECIES_CODE,nvl(HAUL_COUNT,0) HAUL_COUNT,
		nvl(CATCH_COUNT,0) CATCH_COUNT,MEAN_WGT_CPUE CPUE, sqrt(VAR_WGT_CPUE) CPUE_STD
		from ", survey, ".BIOMASS_STRATUM B, (select  A.YEAR, B.STRATUM, C.SPECIES_CODE from
		(select distinct YEAR from ", survey, ".BIOMASS_STRATUM where SURVEY = '", survey, "') A,
		(select STRATUM from GOA.GOA_STRATA where SURVEY = '", survey, "') B, (select SPECIES_CODE
		from GOA.ANALYSIS_SPECIES where BIOMASS_FLAG in ('", survey, "','BOTH')) C) Y where
		B.YEAR(+)=Y.YEAR and B.STRATUM(+)=Y.STRATUM and B.SPECIES_CODE(+)=Y.SPECIES_CODE")
  cat("\nGetting stratum CPUE and variance for analysis species...\n")
  # run query
  dat <- sqlQuery(channel = channel, query = plandata.qry)
  names(dat) <- casefold(names(dat))

  # convert NA CPUEs and SD of CPUEs to 0s
  dat$cpue[is.na(dat$cpue)] = 0
  dat$cpue_std[is.na(dat$cpue_std)] = 0

  if(close.channel)close(channel)

  return(dat)

}
