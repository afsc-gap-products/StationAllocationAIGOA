##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Create Paper Station Logs
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

library(openxlsx)
library(StationAllocationAIGOA)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import a given allocation and order by longitude.
##   GOA: ascending ordering by longitude means West -> East
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
year <- 2025
survey <- "GOA"
total_n <- 400

station_allocation_path <- paste0("G:/", survey, "/",
                                  survey, " ", year, "/Station Allocation/",
                                  tolower(x = survey), "_", year,
                                  "_station_allocation_", total_n, ".xlsx")
output_path <- paste0("G:/RACE_Survey_App/files/Station info/AI_GOA/",
                      "Station logs/Paper logs/",
                      year, " ", survey, " FPC Station Logs.xlsx")

goa_allocated_stations <- openxlsx::read.xlsx(
  xlsxFile = station_allocation_path,
  sheet = "Station Allocation"
)

goa_allocated_stations <-
  goa_allocated_stations[order(goa_allocated_stations$LONGITUDE), ]

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

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   Take the allocation table and format into the form of the station log
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  paper_station_log <-
    data.frame(
      rbind(
        subset(x = goa_allocated_stations,
               subset = VESSEL == names(x = ivessel)
               & STATION_TYPE == "prescribed",
               select = c("STATION", "STRATUM")),
        subset(x = goa_allocated_stations,
               subset = STATION_TYPE == "bonus_stn",
               select = c("STATION", "STRATUM"))
      ),
      "Haul" = "",
      "Trawl with survey gear?" = "",
      "Trawl with tire gear?" = "",
      "New station?" = "",
      "Too steep?" = "",
      "Hard/Rocky" = "",
      "Rolly" = "",
      "Pinnacles" = "",
      "Unnavigable" = "",
      "Snags" = "",
      "Ledges" = "",
      "Cable" = "",
      "Sand waves" = "",
      "Fixed gear" = "",
      "Comments" = "",
      check.names = F )

  names(x = paper_station_log)[names(x = paper_station_log) %in%
                                 c("STATION", "STRATUM")] <-
    c("Station ID", "Stratum")

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   Insert a header row every at a set interval
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  page_break_interval <- 13
  paper_station_log <-
    insert_running_headers(df = paper_station_log,
                           n = page_break_interval,
                           title_column_idx = 9,
                           title = paste(ivessel, year, "FPC Station Log"))

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   Add Sheet and make landscape
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  sheetname <- paste(ivessel, year)
  openxlsx::addWorksheet(wb, sheetname)
  openxlsx::pageSetup(wb, sheet = sheetname, orientation = "landscape")

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   Add data
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   Format cells
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## Set font and text size of the station records, add thin borders
  openxlsx::addStyle(
    wb = wb,
    sheet = sheetname,
    style = createStyle(fontName = "Arial", fontSize = 10,
                        halign = "center", valign = "center",
                        border = "TopBottomLeftRight", borderStyle = "thin"),
    rows = 1:(nrow(x = paper_station_log) + 15),
    cols = 1:ncol(x = paper_station_log),
    gridExpand = TRUE
  )

  ## Rotate headers, set font and text side, bolden text
  openxlsx::addStyle(
    wb = wb,
    sheet = sheetname,
    style = createStyle(textRotation = 90, textDecoration = "bold",
                        fontName = "Arial", fontSize = 10,
                        halign = "center", valign = "center", wrapText = TRUE,
                        border = "TopBottomLeftRight", borderStyle = "medium"),
    rows = seq(from = 2,
               to = nrow(x = paper_station_log) + 15,
               by = page_break_interval + 2),
    cols = 1:(ncol(x = paper_station_log) - 1),
    gridExpand = TRUE
  )

  ## Re-rotate "Comments" header
  openxlsx::addStyle(
    wb = wb,
    sheet = sheetname,
    style = createStyle(fontName = "Arial", fontSize = 10,
                        textDecoration = "bold",
                        halign = "center", valign = "bottom", wrapText = TRUE,
                        border = "TopBottomLeftRight", borderStyle = "medium"),
    rows = seq(from = 2,
               to = nrow(x = paper_station_log) + 15,
               by = page_break_interval + 2),
    cols = ncol(x = paper_station_log),
    gridExpand = TRUE
  )

  ## Fill some of the cell header colors gray
  openxlsx::addStyle(
    wb = wb,
    sheet = sheetname,
    style = createStyle(
      textRotation = 90, fontName = "Arial", fontSize = 10,
      halign = "center", valign = "center", wrapText = TRUE,
      border = "TopBottomLeftRight", borderStyle = "medium",
      textDecoration = "bold", fgFill = "#999DA0"),
    rows = seq(from = 2,
               to = nrow(paper_station_log) + 15,
               by = page_break_interval + 2),
    cols = 4:16,
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
               to = nrow(x = paper_station_log) + 15,
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
               to = nrow(x = paper_station_log) + 15,
               by = page_break_interval + 2),
    cols = 1:ncol(x = paper_station_log),
    gridExpand = TRUE
  )

  ## Highlight bonus stations in blue
  openxlsx::addStyle(
    wb = wb,
    sheet = sheetname,
    style = createStyle(
      fontName = "Arial", fontSize = 10,
      halign = "center", valign = "center",
      border = "TopBottomLeftRight", borderStyle = "thin",
      fgFill = "#90E0EF"),
    rows = which(x = is.na(x = paper_station_log$`Station ID`) &
                   !is.na(x = paper_station_log$Stratum)) + 2,
    cols = 1:ncol(x = paper_station_log),
    gridExpand = TRUE)

  ##~~~~~~~~~~~~~~~~~~~~~~
  ##   Set column widths and row heights. Adjust as needed.
  ##~~~~~~~~~~~~~~~~~~~~~~

  ## Header columns
  openxlsx::setColWidths(wb = wb,
                         sheet = sheetname,
                         cols = 1:ncol(x = paper_station_log),
                         widths = c(12, 5, 7, 7, 7,
                                    rep(3.25,11), 40))

  ## Station records
  openxlsx::setRowHeights(wb = wb,
                          sheet = sheetname,
                          rows = 3:(nrow(x = paper_station_log) + 15),
                          heights = 30)

  ## Title headers
  openxlsx::setRowHeights(wb = wb,
                          sheet = sheetname,
                          rows = seq(from = 1,
                                     to = nrow(x = paper_station_log) + 15,
                                     by = page_break_interval + 2),
                          heights = 15)

  ## Column field headers
  openxlsx::setRowHeights(wb = wb,
                          sheet = sheetname,
                          rows = seq(from = 2,
                                     to = nrow(x = paper_station_log) + 15,
                                     by = page_break_interval + 2),
                          heights = 70)

} ## Loop over vessels -- end

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save workbook
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
openxlsx::saveWorkbook(wb = wb, file = output_path, overwrite = TRUE)
