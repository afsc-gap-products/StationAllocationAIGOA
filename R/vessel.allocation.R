#' Allocate stations across vessels
#'
#' @author Ned Laman \email{ned.laman@@noaa.gov}
#'
#' @param stations dataframe of sampled stations created interanlly in
#' AIGOASurveyPlanning::goa.planning()
#' @param strata dataframe of strata characteristics created by
#' AIGOASurveyPlanning::allocate.effort()
#' @param survey character string. Defaults to "AI" for Aleutian Islands

vessel.allocation <- function(stations, strata, survey){
  ## 12-18-13 N Laman
  ## altered odbcExcel calls to run from XLConnect package

  ## directing traffic to subtask level Survey Planning folder
  if(survey == "GOA")
    survey.name <- survey
  if(survey == "AI")
    survey.name <- "ALEUTIAN"

  cat("Allocating selected gridcells between vessels...\n")

  # Import the vessel data from the Excel spreadsheet
  # vessel.names.filename <- loadWorkbook(paste("g:\\", survey.name, "\\survey planning\\", survey, ".vessel.names.xlsx", sep = ""))
  vessel.names.filename <-
    loadWorkbook(paste0("data/", survey, ".vessel.names.xlsx"))

  vessels <- readWorksheet(vessel.names.filename, sheet = "Sheet1")
  names(vessels) <- gsub("#", ".", names(vessels))

  # Import the vessel strata assignment data from the Excel spreadsheet
  # vessel.strata.assignments.filename <-
  #   loadWorkbook(paste("g:\\", survey.name, "\\survey planning\\", survey,
  #                      ".vessel.strata.assignments.xlsx", sep = ""))
  vessel.strata.assignments.filename <-
    loadWorkbook(paste0("data/", survey, ".vessel.strata.assignments.xlsx"))

  vessel.strata.assignments <- readWorksheet(vessel.strata.assignments.filename, sheet = "Sheet1")
  names(vessel.strata.assignments) <- gsub("#", ".", names(vessel.strata.assignments))

  if(all(vessel.strata.assignments$allocation.type == "random")) {
    allocated.stations.final <- data.frame(rep(vessels$vessel_code, length.out = length(stations$stationid)), stations)
    names(allocated.stations.final) <- c("vessel", names(stations))

  }else{
    num.vessels <- length(vessels$vessel_code)
    # Add the tow cost data to the chosen station data
    stations.names <- names(stations)
    tow.costs <- strata$tow.cost[match(stations$stratum, strata$stratum,  )]
    stations <- data.frame(stations, tow.costs)
    names(stations) <- c(stations.names, "costs")
    # First assign the deepest tows (> 700 m) to vessel 3 and split the 500 - 700 m
    # tows between vessels 2 and 3 and calculate the total costs (in terms of time
    # to each vessel that result from this allocation.
    # Allocate by INPFC area in ascending area for best results.
    # For all the other areas, split the tows equally between the vessels assigned
    # that stratum in "vessels.strata.assignments and assign the tows
    # from north to south within each stratum always starting with the
    # "deeper" boat" and keeping track of costs.
    # For Shumagin and Southeast areas, do a random allocation of stations for
    # all depths where available wire is not an issue, taking into account the costs
    # already incurred by the deep and middle boats above.  Vessel 3 should always
    # be the "deep" boat.
    assignments <- vessel.strata.assignments[match(unique(stations$stratum), vessel.strata.assignments$stratum, nomatch = 0),]
    first.allocation <- data.frame(matrix(nrow = sum(assignments[4:dim(vessel.strata.assignments)[2]]), ncol = 5))
    names(first.allocation) <- c("vessel", "stratum", "stations", "tow.cost", "stratum.cost")
    rownum <- 1

    for(st in unique(stations$stratum)){
      st.stations <- stations[stations$stratum == st,  ]
      num.tows <- length(st.stations[, 1.])
      num.ves <- 1:num.vessels
      vessels.in.stratum <- num.ves[unlist(vessel.strata.assignments[vessel.strata.assignments$stratum == st,
                                                                     4:dim(vessel.strata.assignments)[2]]) == 1]
      tows.needed <- rep(round(num.tows/length(vessels.in.stratum), 0), length(vessels.in.stratum))

      while(sum(tows.needed) != num.tows){

        if(sum(tows.needed) > num.tows)
          tows.needed[length(tows.needed)] <- tows.needed[1] - 1

        if(sum(tows.needed) < num.tows)
          tows.needed[1] <- tows.needed[1] + 1

      }

      first.allocation[rownum:(rownum + length(tows.needed) - 1),  ] <-
        cbind(vessels.in.stratum, st, tows.needed,
              unique(st.stations$cost), tows.needed * unique(st.stations$cost))
      rownum <- rownum + length(tows.needed)

    }

    first.allocation$tow.cost <- as.numeric(first.allocation$tow.cost)
    first.allocation$stations <- as.numeric(first.allocation$stations)
    # This section adjusts the vessel allocations by stratum in
    # an attempt to keep the total "cost" equal among vessels.

    vessels.stratum.costs <- tapply(first.allocation$tow.cost, list(first.allocation$vessel, first.allocation$stratum),sum)
    vessels.costs <- tapply((first.allocation$tow.cost * first.allocation$stations), first.allocation$vessel, sum)
    mean.cost <- mean(vessels.costs)
    cost.differential <- vessels.costs - mean.cost
    recipients.done <- vector()

    while(max(abs(cost.differential)) > min(vessels.stratum.costs, na.rm = T)){
      least.tows <- ceiling(tapply(first.allocation$stations, first.allocation$stratum, sum) * 0.2)
      # tows.available <- t(t(tapply(first.allocation$stations, list(first.allocation$vessel, first.allocation$stratum),
      #	sum)) - least.tows)
      ## tows.available generates non conformable array error
      df <- as.data.frame(t(tapply(first.allocation$stations, list(first.allocation$vessel, first.allocation$stratum),
                                   default = 0, sum)))
      df <- df-least.tows
      tows.available <- t(df)
      total.tows <- tapply(first.allocation$station, list(first.allocation$vessel, first.allocation$stratum),
                           default = 0, sum)

      tows.available[total.tows <= 2] <- 0
      recipient <- as.numeric(names(cost.differential[cost.differential == min(cost.differential)]))

      if(length(recipients.done > 0))
        recipient <- as.numeric(names(cost.differential[cost.differential == min(cost.differential[match(names(cost.differential),
                                                                                                         recipients.done, nomatch = 0) == 0])]))

      potential.donors <- vessels$vessel.number[vessels$vessel.number != recipient]
      bad.donors <- vector()
      next.recipient <- 0
      donor <- as.numeric(names(cost.differential[cost.differential == max(cost.differential)]))
      x <- tows.available[donor,  ][!is.na(tows.available[recipient,  ]) & !is.na(tows.available[donor,  ])]

      while(sum(x) == 0){
        print("bad donor added")
        bad.donors <- c(bad.donors, donor)
        good.donors <- potential.donors[match(potential.donors, bad.donors, nomatch = 0) == 0]

        if(length(good.donors) == 0) {
          next.recipient <- 1
          recipients.done <- c(recipients.done, recipient)
          print(recipients.done)
        }

        break

        donor <- as.numeric(names(cost.differential[good.donors][cost.differential[good.donors] ==
                                                                   max(cost.differential[good.donors])]))
        x <- tows.available[donor,  ][!is.na(tows.available[recipient,  ]) & !is.na(tows.available[donor,  ])]
      }

      if(next.recipient)
        next
      transfer.stratum <- sample(as.numeric(rep(names(x), x)), 1)

      if(first.allocation$tow.cost[first.allocation$vessel == donor & first.allocation$stratum == transfer.stratum] >
         max(cost.differential))
        break

      first.allocation$stations[first.allocation$vessel == donor & first.allocation$stratum == transfer.stratum] <-
        first.allocation$stations[first.allocation$vessel == donor & first.allocation$stratum ==
                                    transfer.stratum] - 1

      first.allocation$stations[first.allocation$vessel == recipient & first.allocation$stratum == transfer.stratum] <-
        first.allocation$stations[first.allocation$vessel == recipient & first.allocation$stratum ==
                                    transfer.stratum] + 1

      vessels.costs <- tapply(X = (first.allocation$tow.cost * first.allocation$stations),
                              INDEX = first.allocation$vessel,
                              FUN = sum)

      cost.differential <- vessels.costs - mean.cost

      if(sum(tows.available[donor,  ], na.rm = T) == 0){
        break
      }
    }

    # Now take the adjusted vessel allocations and assign them to vessels
    # north to south with the shallowest boat getting the northernmost
    # tows.
    for(st in sort(unique(first.allocation$stratum))){
      allocation.type <- vessel.strata.assignments$allocation.type[vessel.strata.assignments$stratum == st]
      st.allocations <- first.allocation[first.allocation$stratum == st,  ]
      num.tows.v <- tapply(st.allocations$stations, st.allocations$vessel, sum)
      st.stations <- stations[stations$stratum == st,  ]

      if(allocation.type == "random"){
        st.stations <- data.frame(rep(st.allocations$vessel_code,
                                      st.allocations$stations),
                                  st.stations)
      }

      if(allocation.type == "east-west"){
        survey.supsmu <- supsmu(as.numeric(stations$center.long),
                                as.numeric(stations$center.lat))
        st.stations <- st.stations[rev(order(as.numeric(st.stations$center.lat) - approx(survey.supsmu$x,
                                                                                         survey.supsmu$y, st.stations$center.long)$y)),  ]
        assigned.vessels <- rep(st.allocations$vessel_code, st.allocations$stations)
        st.stations <- cbind(assigned.vessels, st.stations)
      }

      if(allocation.type == "north-south"){
        st.stations <- north.to.south(st.stations, num.tows.v, names(allocated.stations.final))
      }

      names(st.stations) <- c("vessel", "id", "stratum", "latitude",
                              "longitude", "area", "cost")

      if(exists("allocated.stations.final")){
        allocated.stations.final <- rbind(allocated.stations.final, st.stations)
      }else{
        allocated.stations.final <- st.stations
        names(allocated.stations.final) <- names(st.stations)
      }
    }
  }

  allocated.stations.final

}
