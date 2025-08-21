
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
"pred_grid"

