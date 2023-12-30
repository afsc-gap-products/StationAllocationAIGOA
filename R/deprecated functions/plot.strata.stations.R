#' Plot effort across stratum
#'
#' @description
#' Creates three barplots:
#'
#' 1) Number of stations allocated to each stratum
#'
#' 2) Density of station allocations across strata
#'
#' 3) Density of station allocations across depth strata
#'
#' @author Ned Laman \email{ned.laman@@noaa.gov}
#'
#' @param number.tows numeric. Vector output from
#' AIGOASurveyPlanning::allocate.effort
#' @param strata Dataframe. Output from AIGOASurveyPlanning::allocate.effort
#' @param survey character string. Defaults to "AI" for Aleutian Islands
#'

plot.strata.stations <- function (number.tows, strata, survey) {

  barplot(height = number.tows,
          names.arg = names(number.tows),
          ylab = "Number of Stations",
          main = "Stations allocated per Stratum",
          las = 2, cex.names = 0.7)
  box()
  # locator(n = 1)
  barplot(height = number.tows/(strata$area[match(names(number.tows),
                                                  strata$stratum)]/1000),
          names.arg = names(number.tows),
          las = 2, main = "Stations allocated/1000 km sq", cex.names = 0.7,
          xlab = "Stratum", ylab = "Stations allocated/1000 km sq")
  box()

  # locator(n = 1)
  allocation.density.by.area <-
    tapply(X = number.tows,
           INDEX = strata$inpfc_area[match(names(number.tows),
                                           strata$stratum)],
           FUN = sum)/tapply(X = strata$area/1000,
                             INDEX = strata$inpfc_area[match(names(number.tows),
                                                             strata$stratum)],
                             FUN = sum)
  if (survey == "GOA")
    plot.order <- c(3, 1, 2, 5, 4)
  if (survey == "AI")
    plot.order <- c(3, 2, 1, 4)
  if (survey == "EBSSlope")
    plot.order <- c(1, 2, 3, 4, 5, 6)
  barplot(height = allocation.density.by.area[plot.order],
          names.arg = names(allocation.density.by.area[plot.order]),
          main = "Stations allocated/1000 km sq by INPFC area",
          xlab = "INPFC Area", ylab = "Stations allocated/1000 km sq")
  box()
  # locator(n = 1)
  allocation.density.by.depth <-
    tapply(X = number.tows,
           INDEX = strata$summary_depth[match(names(number.tows),
                                              strata$stratum)],
           FUN = sum)/tapply(X = strata$area/1000,
                             INDEX = strata$summary_depth[match(names(number.tows),                             strata$stratum)],
                             FUN = sum)
  barplot(height = allocation.density.by.depth,
          names.arg = names(allocation.density.by.depth),
          las = 2, main = "Stations allocated/1000 km sq by depth",
          xlab = "Depth", ylab = "Stations allocated/1000 km sq")
  box()
  # locator(n = 1)
}
