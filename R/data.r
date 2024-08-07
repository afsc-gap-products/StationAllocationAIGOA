
#' Predicted densities: Gulf of Alaska Groundfish, 1996-present
#'
#' Predicted densities (kg/km2) from VAST for 15 species/complexes from
#' 1996-present (excluding 2001).
#'
#' @format A 3-dimensional array of predicted densities, units kg/km2
#' \describe{
#'   \item{grid cell}{VAST interpolation grid cell, denoted by subscript g}
#'   \item{taxon}{species/complex, denoted by subscript c}
#'   \item{time}{year, denoted by subscript t}
#'   }
#'
"D_gct"

#' Gulf of Alaska (GOA) Bottom Trawl Survey Stations
#'
#' Includes trawlability information from 2023
#'
#' @format A data frame with 10000 rows and 2 variables:
#' \describe{
#'   \item{GRIDID}{Position of the 5-km grid cell }
#'   \item{NMFS_AREA}{NMFS statistical area name.}
#'   \item{REP_AREA}{Numeric code for NMFS statistical areas.}
#'   \item{STRATUM}{Numeric code for GOA stratum.}
#'   \item{TRAWLABLE}{Character string denoting trawlability status. Valid
#'   values: 'Y' for trawlable, 'N' for untrawlable, and 'UNK' for unknown.}
#'   \item{AREA_KM2}{Numeric, station area in square kilometers}
#' }
#'
"goa_stations"

#' VAST interpolation grid
#'
#' data frame of VAST interpolation grid cells used in the optimization of
#' station effort across Gulf of Alaska bottom trawl survey strata (2025-on).
#'
#' @format A data frame with fields:
#' \describe{
#'   \item{ID}{Integer. Unique cell id.}
#'   \item{Area_km2}{Numeric. total area of the cell in square kilometers.}
#'   \item{Depth_m}{Numeric. Depth of the cell in meters}
#'   \item{Eastings}{Numeric. Eastings in meters using NAD83 / Alaska Albers (EPSG:3338)}
#'   \item{Northings}{Numeric. Northings in meters using NAD83 / Alaska Albers (EPSG:3338)}
#'   \item{NMFS_AREA}{NMFS statistical area name.}
#'   \item{STRATUM}{Numeric code for GOA stratum.}
#' }
"optim_df"

#' 2025 Gulf of Alaska Stratum Boundaries
#'
#' data frame loop up table of depth stratum boundaries by NMFS management area used for the Gulf of Alaska bottom trawl survey design (2025-on).
#'
#' @format A data frame with fields:
#' \describe{
#'   \item{NMFS_AREA}{NMFS statistical area name.}
#'   \item{REP_AREA}{Numeric code for NMFS statistical areas.}
#'   \item{STRATUM}{Numeric code for GOA stratum.}
#'   \item{DEPTH_MIN_M}{Numeric. Minimum boundary of the stratum in meters.}
#'   \item{DEPTH_MAX_M}{Numeric. Maximum boundary of the stratum in meters.}
#'   \item{AREA_KM2}{Numeric. Total area of the stratum in square kilometers.}
#'   \item{PERIM_KM}{Numeric. Total perimeter of the stratum in kilometers.}
#'   \item{USED}{Boolean. Is the stratum currently used in the GOA survey footprint? This field is relevant when considering the GOA survey footprint (up to 700 m or 1000 m?).}
#' }
"goa_stratum_boundaries"
