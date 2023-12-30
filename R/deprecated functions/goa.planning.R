#' Main wrapper function for planning GAP STRS survey in the AI-GOA
#'
#' @description
#' Modifications as of 6/19/2019:
#' Encapsulated creation of the AI_STATIONS list of candidate stations and
#' creating the PLANNING_DATA table in R functions called from goa.planning()
#' and allocate.effort, respectively. This entailed intake of uid, pwd
#' parameters and passing those along with with the survey parameter to those
#' functions. Also changed data constraints for AI stations to include good
#' Ocean Hope stations from 1991 that have yet to be replicated since and
#' conditioned the general rule of using "previously successful tows" with
#' consideration of trawlability field populated in AI.AIGRID_GIS from
#' electronic station logs maintained on the vessel. The routine as written in
#' this project is functional only for the Aleutian Islands. Gulf of Alaska
#' station allocation originates from old SPLUS code curated by Paul von Szalay.
#'
#' Modifications as of 1/20/2016
#' Added a csv output in addition to the xlsx workbook and identified it with
#' number.of.tows
#'
#' Modification as of 12/19/2013
#' Added readline and choose.dir functions to accept parameter inputs. Changed
#' write.table for reporting allocation table to writing an XLSX file
#'
#' Modifications as of 12/18/2013
#' Initially just changing indentation and line breaks to get preferred code
#' form and added ", believeNRows = FALSE" to sqlQuery's to deal with version
#' 3.0.1 of R.
#'
#' @author Ned Laman \email{ned.laman@@noaa.gov}
#'
#'

goa.planning <- function(selection.type = "random", min.size = 5){

  # set options for session
  options(digits = 10)
  options(stringsAsFactors = F)

  if(survey == "AI") {
    ## As of June2019 candidate stations are chosen from previously
    ## successfully trawled stations in RACEBASE.HAUL (including Ocean Hope
    ## stations from 1991) combined with station trawlability recorded
    ## in AI.AIGRID_GIS
    points <- get.ai.stations(channel, schema, pwrd, survey)
    names(points) <- casefold(names(points))
    stations.available <- tapply(points$stationid, points$stratum, length)
  }

  allocations <- allocate.effort(survey,
                                 channel = channel,
                                 points = points,
                                 stations.available = stations.available,
                                 number.of.tows = number.of.tows,
                                 pwrd = pwrd, schema = schema,
                                 output.file = output.file)
  strata <- allocations[[2]]
  allocations <- allocations[[1]]

  stations <- pick.gridpoints(allocations, channel, survey, points)

  survey.plan <- vessel.allocation(stations, strata, survey)
  names(survey.plan) <- toupper(names(survey.plan))
  varTypes <- c("number","varchar2(10)","number",rep("float",2))
  names(varTypes) <- names(survey.plan)

  csv.path <- paste(output.file, "\\", survey, "allocation",
                    number.of.tows, ".csv", sep = "")
  ## cannot get rid of column names when writing a csv unless row.names = T
  write.csv(survey.plan, csv.path, row.names = F)

  ## Save station allocation as an Excel workbook (.xlsx)
  xl.wb.path <- paste(output.file, "\\", survey, "allocation.xlsx", sep = "")
  xl.wb <- loadWorkbook(xl.wb.path, create = TRUE)
  createSheet(xl.wb, name = "StationAllocation")
  writeWorksheet(xl.wb, survey.plan, sheet = "StationAllocation")
  saveWorkbook(xl.wb)

  # assign("default.output.file", output.file, env = .GlobalEnv)

  ## Write copy of survey.plan to an Oracle table for Region & Year Station
  ## Allocation if table does not exist already, a vector of error language
  ## is returned from Oracle in the first if statement if the table does exist
  ## you are prompted to overwrite or not terminus probably needs work
  O.table <- paste0(survey, "_", year, "_STATION_ALLOCATION")
  D <- sqlQuery(channel, paste0("select count(*) n_records from ",
                                survey, "_", year, "_STATION_ALLOCATION"))

  if(is.vector(D)){
    sqlSave(channel, dat = survey.plan, tablename = O.table,
            varTypes = varTypes, rownames = FALSE)
  }else if(is.data.frame(D)){
    qna <- readline(paste0("Do you want to overwrite the", year,
                           " station allocation (enter y or n)? "))
    qna <- (toupper(qna))
    if(qna == "Y"){
      sqlDrop(channel, O.table)
      sqlSave(channel, dat = survey.plan, tablename = O.table,
              varTypes = varTypes, rownames = FALSE)
    }else{
      print("You either opted not to overwrite or entered something other than Y or N above.  Try again!")
      stop()
    }
  }

  odbcClose(channel)

  cat(paste0("Done!  \nResults are found at ", xl.wb.path,
             "\nand in AIGOA_WORK_DATA.", O.table, " in Oracle.\n"))

  invisible()

}
