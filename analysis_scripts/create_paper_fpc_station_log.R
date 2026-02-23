##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Create Paper Station Logs
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
##                Print single-sided
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

library(openxlsx)
library(StationAllocationAIGOA)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import a given allocation and order by longitude.
##   GOA: ascending ordering by longitude means West -> East
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
year <- 2026
survey <- "ALEUTIAN"
survey_short <- "AI"
total_n <- 400

station_allocation_path <- paste0("Y:/RACE_GF/", survey, "/",
                                  survey_short, " ", year, "/Station Allocation/",
                                  tolower(x = survey_short), "_", year,
                                  "_station_allocation_", total_n, "stn.xlsx")
output_path <- paste0("Y:/RACE_GF/RACE_Survey_App/files/Station info/AI_GOA/",
                      "Station logs/Paper logs/", survey_short, "/",
                      year, " ", survey_short, " FPC Station Logs.xlsx")

goa_allocated_stations <- openxlsx::read.xlsx(
  xlsxFile = station_allocation_path,
  sheet = "Station Allocation"
)


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

  subset_allocation <-
    goa_allocated_stations[
      ## order by longitude, convert to degrees W if LONGITUDE > 0
      order(ifelse(test = goa_allocated_stations$LONGITUDE > 0,
                   yes = goa_allocated_stations$LONGITUDE - 360,
                   no = goa_allocated_stations$LONGITUDE),
            decreasing = ifelse(test = survey_short == "AI",
                                yes = FALSE,
                                no = FALSE)),
    ] |> subset(VESSEL == names(x = ivessel) | STATION_TYPE == "bonus")

  ## Reorder stations so that the bonus and new stations are interspersed
  # Use ave() to get the row count within each STRATUM
  row_idx <- ave(seq_len(nrow(subset_allocation)),
                 subset_allocation$STRATUM,
                 FUN = seq_along)

  # Define the priority logic
  priority <- ifelse(row_idx == 1,
                     yes = 1,
                     no = ifelse(
                       test = subset_allocation$STATION_TYPE %in%
                         c("bonus", "new"),
                       yes = 2,
                       no = 3))

  #Reorder the dataframe: sort by STRATUM first, then by our custom priority
  subset_allocation <- subset_allocation[
    order(subset_allocation$STRATUM, priority),
  ]

  #Clear row names
  rownames(x = subset_allocation) <- NULL

  # Add empty rows for alternate stations
  alt_rows <- subset_allocation
  alt_rows[] <- NA

  subset_allocation <- rbind(subset_allocation,
                             alt_rows)
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   Take the allocation table and format into the form of the station log
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  paper_station_log <-
    data.frame(
      subset_allocation |> subset(select = c(STATION, STRATUM)),
      "Haul" = "",
      # "Trawl with survey gear?" = "",
      # "Trawl with tire gear?" = "",
      # "New station?" = "",
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
      "Comments" = ifelse(
        test = is.na(x = subset_allocation$STATION_TYPE),
        yes = "Alt. station for",
        no = ifelse(test = subset_allocation$STATION_TYPE == "assigned",
                    yes = "",
                    no = paste0(subset_allocation$STATION_TYPE, " station:")
        )
      ),
      check.names = F )

  names(x = paper_station_log)[names(x = paper_station_log) %in%
                                 c("STATION", "STRATUM")] <-
    c("Station ID", "Stratum")

  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##   Insert a header row every at a set interval
  ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  page_break_interval <- 13
  paper_station_log <-
    StationAllocationAIGOA::insert_running_headers(
      df = paper_station_log,
      n = page_break_interval,
      title_column_idx = 13,
      title = paste(ivessel, year, "FPC Station Log")
    )

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
    rows = 1:(nrow(x = paper_station_log)+ 2),
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
               to = nrow(x = paper_station_log)+ 2,
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
               to = nrow(x = paper_station_log)+ 2,
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
               to = nrow(paper_station_log)+ 2,
               by = page_break_interval + 2),
    cols = 4:13,
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
               to = nrow(x = paper_station_log)+ 2,
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
               to = nrow(x = paper_station_log)+ 2,
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
    rows = which(x = paper_station_log$Comments == "bonus station:") + 2,
    cols = 1:ncol(x = paper_station_log),
    gridExpand = TRUE)

  openxlsx::addStyle(
    wb = wb,
    sheet = sheetname,
    style = createStyle(
      fontName = "Arial", fontSize = 10,
      halign = "left", valign = "center",
      border = "TopBottomLeftRight", borderStyle = "thin",
      fgFill = "#90E0EF"),
    rows = which(x = paper_station_log$Comments == "bonus station:") + 2,
    cols = ncol(x = paper_station_log),
    gridExpand = TRUE)

  ## For AI surveys, highlight new station in green
  if (survey_short == "AI")
    openxlsx::addStyle(
      wb = wb,
      sheet = sheetname,
      style = createStyle(
        fontName = "Arial", fontSize = 10,
        halign = "center", valign = "center",
        border = "TopBottomLeftRight", borderStyle = "thin",
        fgFill = "#84EAB3"),
      rows = which(x = paper_station_log$Comments == "new station:") + 2,
      cols = 1:ncol(x = paper_station_log),
      gridExpand = TRUE)

  if (survey_short == "AI")
    openxlsx::addStyle(
      wb = wb,
      sheet = sheetname,
      style = createStyle(
        fontName = "Arial", fontSize = 10,
        halign = "left", valign = "center",
        border = "TopBottomLeftRight", borderStyle = "thin",
        fgFill = "#84EAB3"),
      rows = which(x = paper_station_log$Comments == "new station:") + 2,
      cols = ncol(x = paper_station_log),
      gridExpand = TRUE)

  openxlsx::addStyle(
    wb = wb,
    sheet = sheetname,
    style = createStyle(
      fontName = "Arial", fontSize = 10,
      halign = "left", valign = "center",
      border = "TopBottomLeftRight", borderStyle = "thin"),
    rows = which(x = paper_station_log$Comments == "Alt. station for") + 2,
    cols = ncol(x = paper_station_log),
    gridExpand = TRUE)
  ##~~~~~~~~~~~~~~~~~~~~~~
  ##   Set column widths and row heights. Adjust as needed.
  ##~~~~~~~~~~~~~~~~~~~~~~

  ## Header columns
  openxlsx::setColWidths(wb = wb,
                         sheet = sheetname,
                         cols = 1:ncol(x = paper_station_log),
                         widths = c(12, 5, 7,
                                    rep(3.25,10), 63))

  ## Station records
  openxlsx::setRowHeights(wb = wb,
                          sheet = sheetname,
                          rows = 3:(nrow(x = paper_station_log)+ 2),
                          heights = 30)

  ## Title headers
  openxlsx::setRowHeights(wb = wb,
                          sheet = sheetname,
                          rows = seq(from = 1,
                                     to = nrow(x = paper_station_log)+ 2,
                                     by = page_break_interval + 2),
                          heights = 15)

  ## Column field headers
  openxlsx::setRowHeights(wb = wb,
                          sheet = sheetname,
                          rows = seq(from = 2,
                                     to = nrow(x = paper_station_log)+ 2,
                                     by = page_break_interval + 2),
                          heights = 70)

  # Page X of Y in center footer
  setHeaderFooter(wb = wb,
                  sheet = sheetname,
                  footer = c(NA, "Page &[Page] of &[Pages]", NA))
} ## Loop over vessels -- end

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save workbook
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
openxlsx::saveWorkbook(wb = wb, file = output_path, overwrite = TRUE)
