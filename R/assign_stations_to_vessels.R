#' Assign stations to vessels
#'
#' @description Assign stations to vessels
#' @param drawn_stns A dataframe containing picked stations
#' @param order_by The column string to order drawn_stns before assignment
#' @param vessel_ids vector of vessel codes to assign
#'
#' @export
#'

assign_stations_to_vessels <- function(
    drawn_stns,
    order_by = "STRATUM",
    vessel_ids = c(148, 176)
) {

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##  Assign stations to vessels w/o any laning first
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  drawn_stns <- drawn_stns[order(drawn_stns[, order_by]), ]
  drawn_stns$VESSEL <- vessel_ids

  return(drawn_stns)
}
