##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Inventory stations by trawlability status
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import libraries and connect to Oracle
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Import Packages
library(terra)
library(gapindex)

sql_channel <- gapindex::get_connected()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Data
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Import station trawlability info. This shapefile contains both the
## historical and 2025 stations
station_trawl_info <-
  terra::vect(x = "data/GOA/processed_shapefiles/station_trawl_info.shp")

## Import historical tow path lines
towpaths <- terra::vect(x = paste0(dirname(path = getwd()),
                                   "/Globe/outputs/towpaths/towpaths.shp"))
towpaths <- terra::project(x = towpaths, terra::crs(x = station_trawl_info))

## Import haul data from 1990. This is the start of the time series so
## we will only use towpaths from then.
goa_hauls_from_1990 <-
  RODBC::sqlQuery(channel = sql_channel,
                  query = "SELECT FLOOR(CRUISE / 100) AS YEAR,
                           HAULJOIN, STATIONID,
                           CASE
                            WHEN PERFORMANCE >= 0 THEN 'TRUE'
                            WHEN PERFORMANCE < 0 THEN 'FALSE'
                           END AS PERFORMANCE
                           FROM RACEBASE.HAUL
                           WHERE REGION = 'GOA'
                           --AND CRUISE >= 199000
                           ORDER BY YEAR")
names(towpaths) <- "HAULJOIN"
towpaths <- merge(x = towpaths,
                  y = goa_hauls_from_1990,
                  by = "HAULJOIN")

## How many unique stations are in the historical vs 2025 stations
tapply(X = station_trawl_info$GOAGRID_ID,
       INDEX = station_trawl_info$DESIGN_YEA,
       FUN = function(x) length(unique(x)))
## There are 21587 historical stations and 24333 new stations

## For historical stations that consist of the entire grid cell, we can
## assume that the new stations that are contained within those grid cells will
## inherit the trawlability information of those cells.

historical_single_stn_cells <-
  with(subset(x = as.data.frame(station_trawl_info),
              subset = DESIGN_YEA == 1984),
       names(x = which(x = table(STATIONID) == 1)))
length(x = historical_single_stn_cells)
#There 9141 grid cells that only contain 1 historical station

## The remaining grid cells were historically characterized by multiple stations
historical_multi_stn_cells <-
  with(subset(x = as.data.frame(station_trawl_info),
              subset = DESIGN_YEA == 1984),
       names(x = which(x = table(STATIONID) != 1)))
length(x = historical_multi_stn_cells)
#There are 5181 grid cells that were historically characterized by mulitple
#stations

## Of these 5181 grid cells with multiple historical stations, there will be
## instances where all of the stations within those cells have the same
## trawlability information.
multi_stn_cells_same_trawl_info <-
  with(subset(x = as.data.frame(station_trawl_info),
              subset = DESIGN_YEA == 1984 &
                STATIONID %in% historical_multi_stn_cells),
       names(x = which(tapply(X = TRAWLABLE,
                              INDEX = STATIONID,
                              FUN = function(x)
                                length(x = unique(x = x))) == 1) ) )

multi_stn_cells_mixed_trawl_info <-
  with(subset(x = as.data.frame(station_trawl_info),
              subset = DESIGN_YEA == 1984 &
                STATIONID %in% historical_multi_stn_cells),
       names(x = which(tapply(X = TRAWLABLE,
                              INDEX = STATIONID,
                              FUN = function(x)
                                length(x = unique(x = x))) != 1) ) )

length(x = multi_stn_cells_same_trawl_info)
length(x = multi_stn_cells_mixed_trawl_info)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Speck dissoluion: for grid cells that are > 5 km2, dissolve the
##   trawlability information of the portions of stations that are < 1 km2 into
##   the larger portion of the station with the different trawlability info.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
new_station_trawl_info <-
  station_trawl_info[station_trawl_info$DESIGN_YEA == 2025 &
                       station_trawl_info$AREA_KM2 > 0]
## The default scenario for a new station is no change, Scenario 1
new_station_trawl_info$FLAG <- 1

for (icell in multi_stn_cells_mixed_trawl_info) { ## loop over cells -- start

  ## Subset stations within icell
  temp_cell <- subset(x = new_station_trawl_info,
                      subset = new_station_trawl_info$STATIONID == icell)

  ## By tabulating stratum records within x, we're essentially
  ## querying stations with multiple records, each record holds
  ## a different trawlability datum.
  mixed_stns <-
    as.numeric(x = names(x = which(x = table(temp_cell$GOAGRID_ID) > 1)))

  for (istation in mixed_stns) { ## Loop over mixed stations -- start
    temp_stn <- subset(x = temp_cell,
                       subset = temp_cell$GOAGRID_ID == istation)
    ## Scenario 2: If there is a speck...
    if (any(temp_stn$AREA_KM2 < 1)) {
      ## Otherwise, subset the station from x
      major_stn <- temp_stn[which.max(x = temp_stn$AREA_KM2)]
      ## and subset the speck station
      speck_stn <- temp_stn[temp_stn$AREA_KM2 < 1]
      if (nrow(x = speck_stn) > 1)
        speck_stn <- terra::combineGeoms(x = speck_stn[1], y = speck_stn[2])
      ## then absorb the speck station into the bigger station
      temp_combined_geo <- terra::combineGeoms(x = major_stn, y = speck_stn)
      temp_combined_geo$FLAG <- 2

      ## and then replace the merged station in new_station_trawl_info
      new_station_trawl_info <-
        new_station_trawl_info[new_station_trawl_info$GOAGRID_ID != istation]
      new_station_trawl_info <- rbind(new_station_trawl_info,
                                      temp_combined_geo)
      print(paste("Station", istation, "in grid cell", icell, "had a speck"))


    } else({
      ## Scenario 3: station is a mixture of T area (with good tows paths)
      ## and either UKN or UT area. Since there is a good tow in the station,
      ## the whole station is turned to T.
      good_tows <- subset(x = towpaths,
                          subset = towpaths$PERFORMANCE == T &
                            towpaths$STATIONID %in%
                            unique(temp_stn$STATIONID))

      ## Query whether there are any good tows in the mixed station
      good_tow_in_station <-
        any(terra::relate(x = good_tows,
                          y = temp_stn,
                          relation = "intersects"))

      ## If so, convert the non-T area in the station as T
      if (good_tow_in_station & any(temp_stn$TRAWLABLE == "Y")) {
        trawl_area <- subset(x = temp_stn,
                             subset = temp_stn$TRAWLABLE == "Y")
        non_trawl_area <- subset(x = temp_stn,
                                 subset = temp_stn$TRAWLABLE != "Y")
        ## if the non-trawlable area consists of UKN and UT areas,
        ## then merge and dissolve into one geometry
        if (nrow(x = non_trawl_area) > 1)
          non_trawl_area <-
          terra::aggregate(x = non_trawl_area)
        temp_combined_geo <-
          terra::combineGeoms(x = trawl_area,
                              y =  non_trawl_area)
        temp_combined_geo$FLAG <- 3
        ## and then replace the merged station in new_station_trawl_info
        new_station_trawl_info <-
          new_station_trawl_info[new_station_trawl_info$GOAGRID_ID != istation]
        new_station_trawl_info <- rbind(new_station_trawl_info,
                                        temp_combined_geo)
        print(paste("Station", istation, "in grid cell", icell,
                    "converted to TRAWLABLE"))
      }
      ## Scenario 4-10: if there are no tows that
      if ((!good_tow_in_station) |
          (good_tow_in_station & !any(temp_stn$TRAWLABLE %in% "Y")) ) {
        larger_area <- temp_stn[which.max(x = temp_stn$AREA_KM2)]
        other_area <- temp_stn[-which.max(x = temp_stn$AREA_KM2)]
        if (nrow(other_area) > 1) {
          other_area$TRAWLABLE <-
            other_area$TRAWLABLE[which.max(x = other_area$AREA_KM2)]
          other_area <- terra::combineGeoms(x = other_area[1],
                                            y =  other_area[2])
          other_area$AREA_KM2 <- terra::expanse(other_area) / 1000 / 1000
        }

        temp_combined_geo <- terra::combineGeoms(x = larger_area,
                                                 y =  other_area)

        ## Scenario 4: if the lesser area has an area > 5 km2, then station is
        ## turned to unknown
        if (other_area$AREA_KM2 > 5) {
          temp_combined_geo$TRAWLABLE <- "UNK"
          temp_combined_geo$FLAG <- 4
        }

        ## Scenario 5:
        if (larger_area$TRAWLABLE =="N" & other_area$AREA_KM2 < 5) {
          temp_combined_geo$TRAWLABLE <- "N"
          temp_combined_geo$FLAG <- 5
        }

        ## Scenario 6:
        if (larger_area$TRAWLABLE == "UNK" & other_area$TRAWLABLE == "N") {
          temp_combined_geo$TRAWLABLE <- "N"
          temp_combined_geo$FLAG <- 6
        }

        ## Scenario 7
        if (larger_area$TRAWLABLE == "UNK" & other_area$TRAWLABLE == "Y") {
          temp_combined_geo$TRAWLABLE <- "UNK"
          temp_combined_geo$FLAG <- 7
        }

        ## Scenario 8:
        if (larger_area$TRAWLABLE == "Y" & other_area$TRAWLABLE == "UNK") {
          temp_combined_geo$TRAWLABLE <- "UNK"
          temp_combined_geo$FLAG <- 8
        }

        ## Scenarios 9 and 10:
        if (larger_area$TRAWLABLE == "Y" & other_area$TRAWLABLE == "N") {
          temp_combined_geo$TRAWLABLE <-
            ifelse(test = larger_area$AREA_KM2 < 5,
                   yes = "N",
                   no = "UNK")
          temp_combined_geo$FLAG <-
            ifelse(test = larger_area$AREA_KM2 < 5,
                   yes = 9,
                   no = 10)
        }

        ## and then replace the merged station in new_station_trawl_info
        new_station_trawl_info <-
          new_station_trawl_info[new_station_trawl_info$GOAGRID_ID != istation]
        new_station_trawl_info <- rbind(new_station_trawl_info,
                                        temp_combined_geo)
        print(paste("Station", istation, "in grid cell", icell,
                    "converted to", temp_combined_geo$TRAWLABLE))
      }
    })
  } ## Loop over mixed stations -- end

} ## Loop over cells -- end

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Recalculate area and perimeter of each new station
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
new_station_trawl_info$AREA_KM2 <-
  terra::expanse(x = new_station_trawl_info) / 1000 / 1000
new_station_trawl_info$PERIMETER_ <-
  terra::perim(x = new_station_trawl_info) / 1000

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
terra::writeVector(x = new_station_trawl_info,
                   filename = paste0("data/GOA/processed_shapefiles/",
                                     "goa_stations_2025_modified.shp"),
                   filetype = "ESRI Shapefile",
                   overwrite = T)
terra::writeVector(x = towpaths,
                   filename = paste0("data/GOA/processed_shapefiles/",
                                     "towpaths.shp"),
                   overwrite = T)
