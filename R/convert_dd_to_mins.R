#' Converts Latitude or Longitude from Decimal Degrees to Minutes?
#'
#' @description Station lat/lon are converted from decimal degrees to
#' minutes for the skipper paper station logs
#'
#' @author Zack Oyafuso \email{zack.oyafuso@@noaa.gov}
#'
#' @param x vector of latitudes or longitudes to convert.
#' @param digits integer indicating the number of decimal places to be used.
#'
#' @export
#'

convert_dd_to_mins <- function(x, digits = 2) {
  round(x = (x - floor(x) * 60 / 100) + floor(x) * 100,
        digits = digits)
}
