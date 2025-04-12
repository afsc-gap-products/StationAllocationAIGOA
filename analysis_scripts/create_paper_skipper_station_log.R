##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Create Paper Station Logs for Skipper
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

library(openxlsx)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import a given allocation and order by longitude.
##   GOA: ascending ordering by longitude means West -> East
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
year <- 2025
survey <- "GOA"
total_n <- 400

goa_allocated_stations <- openxlsx::read.xlsx(
  xlsxFile = paste0("G:/", survey, "/",
                    survey, " ", year, "/Station Allocation/",
                    tolower(x = survey), "_", year,
                    "_station_allocation_", total_n, ".xlsx"),
  sheet = "Station Allocation"
)

goa_allocated_stations <-
  goa_allocated_stations[order(goa_allocated_stations$LONGITUDE), ]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create a function to insert a header row at a set interval
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
insert_running_headers <-
  function(df,
           n = 13,
           title = "Alaska Provider Station Log 202501") {

    out <- data.frame()
    total <- nrow(x = df)

    ## Put the title sort of in the middle of the sheet
    title_row <- paper_station_log[0, ]
    title_row[1, "Latitude"] = title
    header_row <- as.list(x = names(x = df))

    for (i in seq(from = 1, to = total, by = n)) {
      chunk <- df[i:min(i + n - 1, total), ]
      out <- rbind(out, title_row, header_row, chunk)
    }
    names(x = out) <- names(x = df)
    return(out[-1:-2, ])
  }

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create a new workbook
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wb <- openxlsx::createWorkbook()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Loop through vessels and create separate logs
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vessel_names <- c("148" = "OEX", "176" = "AKP")

for (ipage in 1:length(x = vessel_names)) { ## Loop over vessels -- start

  ivessel <- vessel_names[ipage]

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   Take the allocation table and format into the form of the station log
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  paper_station_log <-
    data.frame(subset(x = goa_allocated_stations,
                      select = c("VESSEL", "STATION", "STRATUM",
                                 "LATITUDE", "LONGITUDE")),
               "Comments" = "",
               check.names = F )
  paper_station_log$LATITUDE <-
    StationAllocationAIGOA::convert_dd_to_mins(x = paper_station_log$LATITUDE)
  paper_station_log$LONGITUDE <-
    StationAllocationAIGOA::convert_dd_to_mins(x = paper_station_log$LONGITUDE)

  paper_station_log$VESSEL <- vessel_names[paste(paper_station_log$VESSEL)]

  names(x = paper_station_log) <- c("Vessel", "ID", "Stratum",
                                    "Latitude", "Longitude", "Comments")

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   Insert a header row every at a set interval
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  page_break_interval <- 27
  paper_station_log <-
    insert_running_headers(df = paper_station_log,
                           n = page_break_interval,
                           title = paste0("Skipper's Station Log ", year, "-01"))

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   Add Sheet
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sheetname <- paste(ivessel, year)
  openxlsx::addWorksheet(wb = wb,sheetName = sheetname)

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   Add data
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # Write the station log first without headers
  openxlsx::writeData(wb = wb,
                      sheet = sheetname,
                      x = paper_station_log,
                      startRow = 3,
                      colNames = FALSE)

  # Add subtitle
  openxlsx::writeData(wb = wb,
                      sheet = sheetname,
                      x = paper_station_log[page_break_interval + 1,],
                      startRow = 1,
                      colNames = FALSE)

  # Add header row with rotated headers
  openxlsx::writeData(wb = wb,
                      sheet = sheetname,
                      x = as.list(x = names(x = paper_station_log)),
                      startRow = 2,
                      colNames = FALSE)

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   Format cells
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## Set font and text size of the station records, add thin borders
  openxlsx::addStyle(
    wb = wb,
    sheet = sheetname,
    style = createStyle(fontName = "Arial", fontSize = 10,
                        halign = "center", valign = "center",
                        border = "TopBottomLeftRight", borderStyle = "thin"),
    rows = 1:(nrow(x = paper_station_log) + 30),
    cols = 1:ncol(x = paper_station_log),
    gridExpand = TRUE
  )

  ## Rotate headers, set font and text side, bolden text
  openxlsx::addStyle(
    wb = wb,
    sheet = sheetname,
    style = createStyle(textDecoration = "bold",
                        fontName = "Arial", fontSize = 10,
                        halign = "center", valign = "center", wrapText = TRUE,
                        border = "TopBottomLeftRight", borderStyle = "medium"),
    rows = seq(from = 2,
               to = nrow(x = paper_station_log) + 15,
               by = page_break_interval + 2),
    cols = 1:ncol(x = paper_station_log),
    gridExpand = TRUE
  )

  ## Highlight station specific to the vessel in yellow
  openxlsx::addStyle(
    wb = wb,
    sheet = sheetname,
    style = createStyle(
      fontName = "Arial", fontSize = 10,
      halign = "center", valign = "center",
      border = "TopBottomLeftRight", borderStyle = "thin",
      fgFill = "#FFFF00"),
    rows = which(x = paper_station_log$Vessel == ivessel) + 2,
    cols = 1:ncol(x = paper_station_log),
    gridExpand = TRUE)

  ## Highlight bonus stations in blue
  openxlsx::addStyle(
    wb = wb,
    sheet = sheetname,
    style = createStyle(
      fontName = "Arial", fontSize = 10,
      halign = "center", valign = "center",
      border = "TopBottomLeftRight", borderStyle = "thin",
      fgFill = "#90E0EF"),
    rows = which(x = is.na(x = paper_station_log$Latitude)) + 2,
    cols = 1:ncol(x = paper_station_log),
    gridExpand = TRUE)

  ## Remove borders for the title headers
  openxlsx::addStyle(
    wb = wb,
    sheet = sheetname,
    style = createStyle(fontName = "Arial", fontSize = 10,
                        textDecoration = "bold", wrapText = FALSE,
                        halign = "center", valign = "bottom",
                        borderStyle = "medium", border = "TopBottom"),
    rows = seq(from = 1,
               to = nrow(x = paper_station_log) + page_break_interval,
               by = page_break_interval + 2),
    cols = 1:ncol(x = paper_station_log),
    gridExpand = TRUE
  )

  ## Remove the top border from the title header
  openxlsx::addStyle(
    wb = wb,
    sheet = sheetname,
    style = createStyle(fontName = "Arial", fontSize = 10,
                        halign = "center",  valign = "bottom",
                        borderStyle = "medium", border = "Bottom",
                        wrapText = FALSE, textDecoration = "bold"),
    rows = seq(from = 1,
               to = nrow(x = paper_station_log) + page_break_interval,
               by = page_break_interval + 2),
    cols = 1:ncol(x = paper_station_log),
    gridExpand = TRUE
  )

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   Set column widths and row heights. Adjust as needed.
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## Header columns
  openxlsx::setColWidths(wb = wb,
                         sheet = sheetname,
                         cols = 1:ncol(x = paper_station_log),
                         widths = c(10, 10, 10, 10, 10, 30))

  ## Station records
  openxlsx::setRowHeights(wb = wb,
                          sheet = sheetname,
                          rows = 3:(nrow(x = paper_station_log) +
                                      page_break_interval),
                          heights = 25)

  ## Title headers
  openxlsx::setRowHeights(wb = wb,
                          sheet = sheetname,
                          rows = seq(from = 1,
                                     to = nrow(x = paper_station_log) +
                                       page_break_interval,
                                     by = page_break_interval + 2),
                          heights = 25)

  ## Column field headers
  openxlsx::setRowHeights(wb = wb,
                          sheet = sheetname,
                          rows = seq(from = 2,
                                     to = nrow(x = paper_station_log) +
                                       page_break_interval,
                                     by = page_break_interval + 2),
                          heights = 25)

} ## Loop over vessels -- end

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save workbook
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
openxlsx::saveWorkbook(
  wb = wb,
  file = paste0("analysis_scripts/",
                year, " ", survey, " Skipper Station Logs.xlsx"),
  overwrite = TRUE
)
