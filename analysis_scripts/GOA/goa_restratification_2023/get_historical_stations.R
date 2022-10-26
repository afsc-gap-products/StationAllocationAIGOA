##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Packages
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(terra)
library(RODBC)
library(getPass)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Connect to VPN!
##   Connect to oracle
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
assign(x = "channel",
       value = RODBC::odbcConnect("AFSC",
                                  uid = getPass::getPass("uid"),
                                  pwd = getPass::getPass("pwd")),
       envir = .GlobalEnv)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Query Station data from Oracle
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grid_q <- "select * from GOA.GOAGRID_GIS"
goa_grid_2021 <-  RODBC::sqlQuery(channel = channel, query = grid_q)
attributes(goa_grid_2021)$date.accessed <- Sys.Date()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Query Strata data from Oracle
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
strata_q <- "select * from GOA.GOA_STRATA where SURVEY = 'GOA'"
goa_strata_2021 <-  RODBC::sqlQuery(channel = channel, query = strata_q)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import shapefile of untrawlable areas as of 2019, Land, and sandman reef
##   that potentially coincide with the goa grid
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
untrawl_2019_shp <- terra::vect(x = paste0("data/GOA/shapefiles_from_GDrive/",
                                           "goagrid2019_landuntrawlsndmn.shp"))

stratum_0_IDs <- sort(subset(x = as.data.frame(untrawl_2019_shp),
                             subset = STRATUM == 0)$GOAGRID_ID)

untrawlable_stations_2019 <- sort(subset(x = as.data.frame(untrawl_2019_shp),
                                         subset = STRATUM != 0)$GOAGRID_ID)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Identify untrawlable stations as of 2021
##   Create a shapefile of untrawlable stations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
untrawlable_stations_2021 <-
  sort(subset(x = goa_grid_2021, subset = TRAWLABLE == "N")$GOAGRID_ID)

goa_grid_2021_shp <- terra::vect("data/GOA/shapefiles_from_GDrive/goagrid.shp")

untrawl_2021_shp <-
  goa_grid_2021_shp[goa_grid_2021_shp$GOAGRID_ID %in% untrawlable_stations_2021]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Merge 2019 shapefile of untrawlable areas and
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
untrawl_2021_shp <- terra::union(untrawl_2019_shp[untrawl_2019_shp$STRATUM == 0],
                                 untrawl_2021_shp)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Plot
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_grid <- terra::vect(x = "data/GOA/shapefiles_from_GDrive/goagrid.shp")
goa_domain <- terra::aggregate(x = goa_grid)

pdf(file = "data/GOA/processed_shapefiles/UT_areas_2021.pdf",
    width = 10, height = 8)
plot(untrawl_2019_shp, col = "cyan", border = F, axes = F)
plot(untrawl_2021_shp[!untrawl_2021_shp$GOAGRID_ID %in%
                        c(stratum_0_IDs, untrawlable_stations_2019)],
     col = "red", add = TRUE, border = F)
plot(goa_domain, add = TRUE, lwd = 0.5)
legend("topleft", legend = c("Untrawlable as of 2019,\nLand, and Sandman Reef",
                             "Added in 2021"),
       fill = c("cyan", "red"))
dev.off()


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Upload dfs to repo
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(goa_grid_2021, file = "data/GOA/grid_goa_2021.rds")
saveRDS(goa_strata_2021, file = "data/GOA/grid_strata_2021.rds")

if (!dir.exists("data/GOA/processed_shapefiles/"))
  dir.create("data/GOA/processed_shapefiles/")
terra::writeVector(x = untrawl_2021_shp,
                   filename = paste0("data/GOA/processed_shapefiles/",
                                     "GOA_untrawl_2021.shp"),
                   overwrite=TRUE)
