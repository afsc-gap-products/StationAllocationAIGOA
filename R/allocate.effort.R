#' Allocate stations across strata
#'
#' @description
#' Allocates stations across strata based on: stratum size, historical standard
#' devation of cpue, ex vessel price, and abundance of a preselected subset of
#' species, and the time it takes to complete tows. (Reference a tech memo for
#' how exactly Neyman Allocations are handled?).
#'
#' Updates of code from Ned:
#' Initial changes are to reformat code into my preferred configuration
#' Added ", believeNRows = FALSE" to odbc calls for version 3.0.1 of R
#'
#' Added package XLConnect to provide functionality to to read Excel
#' workbooks under version 3.0.1 of R in 64-bit environment
#'
#' Created 12/18/13
#'
#' @author Ned Laman \email{ned.laman@@noaa.gov}
#' @import tcltk XLConnect
#'
#' @param survey character string. Defaults to "AI" for Aleutian Islands
#' @param channel open connection object. Created from RODBC::odbcConnect()
#' @param schema character string. SQL schema name. Defaults to AIGOA_WORK_DATA
#' @param pwrd character string. SQL password
#' @param points dataframe of potential station locations queried from SQL,
#' called in AIGOASurveyPlanning::goa.planning().
#' @param stations.available integer. Vector of available stations in each
#' stratum.
#' @param number.of.tows integer. Number of total stations.
#' @param output.file character string. Output directory.
#'
#' @return A list with a vector of allocated stations across strata
#' and a dataframe with relevant fields:\tabular{ll}{
#'    \code{survey} \tab character. Location of survey. \cr \tab \cr
#'    \code{stratum} \tab integer. Stratum id. \cr \tab \cr
#'    \code{area} \tab numeric. Total area (km2?). \cr \tab \cr
#'    \code{inpfc_area} \tab character string. Management area as defined by the
#'     International North Pacific Fisheries Commission. \cr \tab \cr
#'     \code{inpfc_area} \tab character string. Additional information about
#'      the stratum. \cr \tab \cr
#'    \code{min_depth} \tab integer. minimium depth of the stratum. \cr \tab \cr
#'    \code{max_depth} \tab integer. maximium depth of the stratum. \cr \tab \cr
#'    }
#' }

allocate.effort <- function(survey,
                            channel = NA,
                            schema,
                            pwrd,

                            points,
                            stations.available,
                            number.of.tows = 741,

                            output.file){

  if(is.na(channel)){
    require(RODBC)
    channel <- odbcConnect(dsn = "AFSC", uid = schema, pwd = pwrd, believeNRows = FALSE)
    close.channel = TRUE
  }else{
    close.channel <- FALSE
  }

  # require(tcltk)
  # require(XLConnect)

  ## directing traffic to subtask level Survey Planning folder
  if(survey == "GOA")
    survey.name <- survey
  if(survey == "AI")
    survey.name <- "ALEUTIAN"

  cat("\nGetting planning data from Oracle...\n")
  planning.data <- get.planning.data(channel = channel,
                                     survey = survey,
                                     pwrd = pwrd,
                                     schema = schema)

  ## Import the vessel strata assignment data from the Excel spreadsheet
  ## Error here in 2013 allocation where xls worksheet had an extra column
  ## that was blank in the last position yielding a vector of NAs that blew
  ## stuff up downstream
  vessel.strata.assignments.filename <-
    loadWorkbook(paste("g:\\", survey.name, "\\survey planning\\",
                       survey, ".vessel.strata.assignments.xlsx",
                       sep = ""))
  vessel.strata.assignments <-
    readWorksheet(vessel.strata.assignments.filename,
                  sheet="Sheet1")

  planning.strata <- vessel.strata.assignments$stratum[apply(vessel.strata.assignments[4:dim(vessel.strata.assignments)[2]], 1, sum) > 0]
  planning.data <- planning.data[!is.na(match(planning.data$stratum, planning.strata)),  ]
  #
  # Figure out what years are available in the Oracle planning data and present
  # the user with the choice of using these years and accept the input.
  #
  potential.years <- rev(sort(unique(planning.data$year)))
  years <- potential.years[potential.years > 1989]
  years.vector <- rev(sort(years))
  #
  # Invoke user input for the year weighting factor.
  #
  year.weights <- rep(1, length(years))
  names(year.weights) <- as.character(years.vector)
  #
  # Number of survey days is now implicit in number of stations which is based on average tows per day
  #

  cat(paste("\nThe estimated number of tows is", number.of.tows, "\n\n"))
  #
  # Get strata information from goa.goa_strata.
  #
  cat("\nGetting strata information from Oracle...\n")

  if(survey == "AI" | survey == "GOA") {
    strata <- sqlQuery(channel = channel,
                       query = paste("select * from goa.goa_strata where survey ='",
                                     survey, "'", sep = ""),
                       believeNRows = FALSE)
  }

  if(survey == "EBSSlope") {
    strata <- sqlQuery(channel = channel,
                       query = "select * from ebsslope.slope_strata",
                       believeNRows = FALSE)
  }

  names(strata) <- casefold(names(strata))
  strata <- strata[!is.na(match(strata$stratum, planning.strata)),  ]

  # for GOA there are two tow.costs columns but for AI there is presently (12/2013) only one
  if(survey == "GOA") {
    tow.costs <- time.costs(survey = survey, channel = channel)
    tow.costs <- tow.costs[!is.na(match(tow.costs$stratum, planning.strata)),  ]
    strata <- data.frame(strata[order(strata$stratum),  ], tow.costs$cost[order(tow.costs$stratum)])
  }
  # assigns equal weighting = 1 to each stratum
  if(survey == "AI" | survey == "EBSSlope") {
    strata <- data.frame(strata[order(strata$stratum),  ],
                         rep(1, length(strata$stratum)),
                         rep(1, length(strata$stratum)))
  }

  names(strata)[length(names(strata))] <- "tow.cost"

  #
  # Import excel spreadsheet controlling what species to include in the analysis
  #
  cat("\nRetrieving species to use in analysis from excel spreadsheet...\n")

  # planning.species.filename <- paste("g:\\", survey.name, "\\survey planning\\",
  #                                    survey, ".planning.species.xlsx", sep = "")
  planning.species.filename <- paste0("data/", survey, ".planning.species.xlsx")

  planning.species.channel <- loadWorkbook(planning.species.filename)
  planning.species <- readWorksheet(planning.species.channel, sheet = "Sheet1")

  names(planning.species) <- gsub("#", ".", names(planning.species))
  planning.data <- planning.data[!is.na(match(planning.data$species_code, planning.species$species.code[planning.species$include == "Y"])),]
  planning.data$cpue[is.na(planning.data$cpue)] <- 0
  planning.data$cpue_std[is.na(planning.data$cpue_std)] <- 0
  #
  # Create answer matrices and arrays
  #
  answer.matrix <- matrix(nrow = length(strata$stratum), ncol = length(sort(unique(planning.data$species_code))),
                          dimnames = list(strata$stratum, sort(unique(planning.data$species_code))))
  answer.years <- array(dim = c(length(strata$stratum), length(sort(unique(planning.data$species_code))), length(years.vector)),
                        dimnames = list(strata$stratum, sort(unique(planning.data$species_code)), years.vector))
  answer <- matrix(nrow = length(strata$stratum), ncol = length(years.vector), dimnames = list(strata$stratum, years.vector))

  cat("\nAllocating effort between strata...\n")
  #
  # Loop through each year
  #
  for(year in 1.:length(years.vector)) {

    year <- years.vector[year]
    cat(paste("Now processing year ", year, "...\n", sep = ""))
    year.data <- planning.data[planning.data$year == year,  ]
    #
    # Loop through each species within each year
    #
    weight <- 1.
    weight.factor <- vector()
    #
    # For each species within each year, calculate an optimal survey plan
    #
    for(spp in (sort(unique(planning.data$species_code)))) {

      species.data <- planning.data[planning.data$species_code == spp,  ]
      year.species.data <- year.data[year.data$species_code == spp,  ]
      year.species.data <- year.species.data[order(year.species.data$stratum),  ]
      #
      # Figure out which strata have no tows and assign mean cpue and std to them
      #
      # as of 12/24/13 the min.depths definition was incorrect (strata$min.depth instead of strata$min_depth) and returned a null
      # suspect that min.depths is not used again but maybe it shows up somewhere else.
      min.depths <- strata$min_depth[match(year.species.data$stratum, strata$stratum)]
      no.cpue.strata <- year.species.data$stratum[is.na(year.species.data$cpue)]

      if(length(no.cpue.strata) > 0.) {
        # stratum average of CPUE standard deviation
        year.species.data$cpue_std[is.na(year.species.data$cpue)] <- tapply(species.data$cpue_std, species.data$stratum, mean,
                                                                            na.rm = T)[as.character(no.cpue.strata)]
        # stratum average of CPUE
        year.species.data$cpue[is.na(year.species.data$cpue)] <- tapply(species.data$cpue, species.data$stratum, mean,
                                                                        na.rm = T)[as.character(no.cpue.strata)]
      }

      no.cpue.std.strata <- year.species.data$stratum[is.na(year.species.data$cpue_std)]

      if(length(no.cpue.std.strata) > 0.) {
        # stratum average of CPUE standard deviation
        year.species.data$cpue_std[is.na(year.species.data$cpue_std)] <- tapply(species.data$cpue_std, species.data$stratum, mean,
                                                                                na.rm = T)[as.character(no.cpue.std.strata)]
        # if stratum CPUE is 0 but stratum CPUE STD is > 0, then CPUE = mean CPUE STD for that stratum
        if(length(year.species.data$cpue[year.species.data$cpue == 0. & year.species.data$cpue_std > 0.]) > 0.) {
          no.cpue.std.strata <- year.species.data$stratum[year.species.data$cpue == 0. & year.species.data$cpue_std > 0.]
          year.species.data$cpue[year.species.data$cpue == 0. & year.species.data$cpue_std > 0.] <-
            tapply(species.data$cpue[!is.na(species.data$cpue_std)], species.data$stratum[!is.na(species.data$cpue_std)], mean,
                   na.rm = T)[as.character(no.cpue.std.strata)]
          year.species.data$cpue.std[is.na(year.species.data$cpue_std)] <- 0.
        }
      }

      year.species.areas <- strata$area[match(year.species.data$stratum, strata$stratum)]
      #
      # Create a cost vector, given cost per stratum
      # for AI the stratum cost is 1.0 as of 12/24/13
      costs <- sqrt(strata$tow.cost[match(year.species.data$stratum, strata$stratum)])
      #
      # Calculate the optimal allocation as a fraction of the total survey effort
      # for this species, given this year's data.
      #

      # (stratum area (km2) * stratum mean cpue standard deviation) divided by the stratum costs
      var.numerator <- ((year.species.areas * year.species.data$cpue_std)/costs)
      # total standard deviation by area and stratum over costs
      var.denominator <- sum(var.numerator)

      if(var.denominator == 0)
        var.denominator <- 1
      #
      # Express this as a number of tows and store answer in answer.matrix.
      #
      # essentially applying tow allocation on the basis of proportional contribution by area to overall variance of the data
      answer.matrix[as.character(year.species.data$stratum), as.character(spp)] <- (var.numerator/var.denominator) * number.of.tows
      #
      # Calculate the species weighting factor (ex vessel price * abundance is used here, but can be anything)
      #

      ex.vessel.price <- planning.species$ex.vessel.price[planning.species$species.code == spp]
      year.species.data$cpue[is.na(year.species.data$cpue)] <- 0
      # weighting factor is essentially the sum of the mean biomass cpue scaled up to area and weighted by the ex.vessel.price/1000
      weight.factor[weight] <- sum(year.species.areas * year.species.data$cpue) * (ex.vessel.price/1000)
      # code above throws warnings because ex.vessel.price is a vector of length n that typically has one value and a bunch of NAs
      # code below eliminates the NAs and warnings, but at this point I do not know the cascading implications of making this change
      # weight.factor[weight] <- sum(year.species.areas * year.species.data$cpue) * (ex.vessel.price[!is.na(ex.vessel.price)]/1000)
      weight <- weight + 1.
    }

    #
    # Weight the results and calculate the weighted optimal allocation for that year
    #
    answer.matrix <- t(weight.factor * t(answer.matrix))
    answer.matrix <- (answer.matrix/sum(answer.matrix)) * number.of.tows
    #
    # Write the answer to the arrays.
    #
    answer[as.character(year.species.data$stratum), as.character(year)] <- apply(answer.matrix, 1., sum)
    answer.years[,  , as.character(year)] <- answer.matrix
  }

  answer.years <- (answer.years/sum(answer.years)) * number.of.tows
  number.tows <- apply(answer, 1, weighted.mean, year.weights)
  # previous to 12/24/13 sample density returned a null because strata$area was spelled strata$areas
  sample.density <- number.tows/(strata$area/1000.)
  #
  # Present the optimized allocation results to the user.
  #
  cat("\nHere are the optimized sampling densities (tows/1000km2) by stratum:\n")
  print(round(number.tows/(year.species.areas/1000.), 1.))
  cat("\nHere are the optimized number of tows by stratum:\n")
  print(round(number.tows, 0))

  stations.unavailable <- number.tows - stations.available
  #
  # Make sure the number of allocated tows doesn't exceed the number of available tows by more than 2.
  # If it does, set the allocated number to the available number + 2 and redistribute effort to other strata
  #
  if(any(stations.unavailable > 0.5))
    cat(paste("Unadjusted allocation would require", round(stations.unavailable[stations.unavailable > 0.5], 0),
              "additional tows from stratum", names(stations.unavailable[stations.unavailable > 0.5]), "\n"))

  while(any(stations.unavailable > 2.5)) {
    cat("\nAdjusting allocation so number of stations allocated does not exceed the number of available stations + 2\n")
    number.tows[(number.tows - stations.available) > 2.5] <- stations.available[(number.tows - stations.available) > 2.5] + 2
    number.tows <- number.tows * (number.of.tows/sum(number.tows))
    stations.unavailable <- number.tows - stations.available
  }
  #
  # Make sure no stratum will get less than 2 tows
  #
  if(any(number.tows < 2)) {
    cat("\nAdjusting allocation so no strata receives less than 2 tows.\n")
    number.tows[number.tows < 1.5] <- 2.
  }
  #
  # Make sure the total number of tows equals what is requested
  #
  number.tows <- number.tows * (number.of.tows/sum(number.tows))
  most.deserving <- number.tows - round(number.tows, 0)

  while(sum(round(number.tows, 0)) < number.of.tows) {
    cat("\nAdding station in stratum", names(number.tows[most.deserving == max(most.deserving)]), "to account for rounding.\n")
    number.tows[most.deserving == max(most.deserving)] <- number.tows[most.deserving == max(most.deserving)] + 1
    most.deserving[most.deserving == max(most.deserving)] <- most.deserving[most.deserving == max(most.deserving)] - 1
  }

  most.deserving[number.tows < 2.5] <- 0

  while(sum(round(number.tows, 0)) > number.of.tows) {
    cat("\nRemoving station in stratum", names(number.tows[most.deserving == min(most.deserving)]), "to account for rounding.\n")
    number.tows[most.deserving == min(most.deserving)] <- number.tows[most.deserving == min(most.deserving)] - 1
    most.deserving[most.deserving == min(most.deserving)] <- most.deserving[most.deserving == min(most.deserving)] + 1
  }
  #
  # Present the optimized, adjusted allocation results to the user.
  #
  cat("\nHere are the optimized and adjusted sampling densities (tows/1000km2) by stratum:\n")
  print(round(number.tows/(year.species.areas/1000.), 1.))
  cat("\nHere are the optimized and adjusted number of tows per stratum:\n")
  print(round(number.tows, 0))
  #
  # Plot the data, if desired, and output the final answer.
  #
  #
  # cat("\nWould you like to plot the data? ")
  # plot.answer <- readline()
  # if(plot.answer == "Y" || plot.answer == "y")
  plot.allocations.by.stratum(answer.years, number.tows, strata, channel, survey, output.file)
  list(round(number.tows, 0), strata)
}
