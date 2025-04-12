#' Create a function to insert a header row at a set interval
#'
#' @description Inserts running breaks for header and title rows during
#' the production of the paper FPC and Skipper Station Log Sheets
#'
#' @author Zack Oyafuso \email{zack.oyafuso@@noaa.gov}
#'
#' @param df dataframe of the station log
#' @param n interger. At which inverval should a page break be made? Default
#' is 13, means that title and header headers are inserted every 13 records.
#' @param title_column_idx integer. Column index where the title should be placed
#' @param title string. Running title of the log, e.g.,
#' "Alaska Provider FPC Station Log 202501"
#'
#' @export
#'

insert_running_headers <-
  function(df,
           n = 13,
           title_column_idx = 4,
           title = "Alaska Provider FPC Station Log 202501") {

    out <- data.frame()
    total <- nrow(x = df)

    ## Put the title sort of in the middle of the sheet
    title_row <- paper_station_log[0, ]
    title_row[1, title_column_idx] = title
    header_row <- as.list(x = names(x = df))

    for (i in seq(from = 1, to = total, by = n)) {
      chunk <- df[i:min(i + n - 1, total), ]
      out <- rbind(out, title_row, header_row, chunk)
    }
    names(x = out) <- names(x = df)
    return(out[-1:-2, ])
  }
