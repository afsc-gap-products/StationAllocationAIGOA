#' Randomly sample stations
#'
#' @description
#' Based on a particular allocation of stations across strata, randomly sample
#' station locations.
#'
#' @author Ned Laman \email{ned.laman@@noaa.gov}
#'
#' @param allocations integer. Named vector of station allocations across strata..
#' @param survey character string. Defaults to "AI" for Aleutian Islands
#' @param points dataframe of potential station locations queried from SQL,
#' called in AIGOASurveyPlanning::goa.planning().
#'

pick.gridpoints <- function(allocations,
                            user.name.password,
                            survey,
                            points,
                            min.size = 5)
{
  points <- points[!is.na(match(points$stratum, names(allocations))),  ]
  row.names(points) <- 1:length(points$stationid)
  cat("Randomly picking gridcells for each stratum...\n")
  all.samples <- vector()

  # Use gridcells that are larger than min.size (default 5 km2) and have not
  # been identified as untrawlable or in the wrong stratum.
  if(survey == "GOA") points$area_km2 <- as.numeric(points$area_km2)

  # For each stratum, randomly select the allocated number of stations from
  # the available pool with the probability of selection proportional to the
  # gridcell area.

  additional.rows <- vector()
  add.row <- 1
  for(i in as.numeric(names(allocations))) {

    stratum.data <- points[points$stratum == i,  ]
    point.rows <- as.numeric(row.names(stratum.data))
    sample.size <- allocations[as.character(i)]

    if(survey == "GOA")
      stratum.sample <- sample(x = point.rows,
                               size = sample.size,
                               prob = stratum.data$area_km2)

    if(survey == "AI" | survey == "EBSSlope") {
      if(length(point.rows) < sample.size) {
        tows.missing <- sample.size - length(point.rows)
        additional.rows[add.row:(add.row + tows.missing - 1)] <-
          rep(i, tows.missing)
        add.row <- add.row + tows.missing
        cat(paste("You will need to locate", tows.missing,
                  "extra tow(s) in stratum", i, "\n"))
        sample.size <- length(point.rows)
      }
      if(sample.size == 0)
        next
      stratum.sample <- sample(point.rows, sample.size)
    }
    all.samples <- c(all.samples, stratum.sample)
  }

  all.samples <- points[all.samples,  ]
  additional.nas <- rep(NA, length(additional.rows))
  additional.rows <- data.frame(additional.nas, additional.rows, additional.nas, additional.nas, additional.nas)
  names(additional.rows) <- names(all.samples)
  all.samples <- rbind(all.samples, additional.rows)
  all.samples[order(all.samples$stratum),  ]
}
