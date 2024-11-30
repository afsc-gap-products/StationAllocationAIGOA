##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       GOA Station Allocation
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   For year 2025
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Packages, laning areas around Kodiak
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# devtools::install_github(repo = "afsc-gap-products/StationAllocationAIGOA")
library(StationAllocationAIGOA)
library(terra)
library(akgfmaps)
library(RColorBrewer)

laning_area <-
  terra::vect(x = "data/GOA/shapefiles_akgfmaps/goa_laning_area.shp")

goa_stations <-
  terra::vect(x = "data/GOA/shapefiles_akgfmaps/goa_stations_2025.gpkg")
goa_strata <-
  terra::vect(x = "data/GOA/shapefiles_akgfmaps/goa_strata_2025.gpkg")
goa_stations[, c("x", "y")] <-
  terra::crds(x = terra::centroids(x = goa_stations,
                                   inside = TRUE))

## `goa_base` are basic shape layers from the akgfmaps package
goa_base <- akgfmaps::get_base_layers(select.region = "goa",
                                      set.crs = "EPSG:3338")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Calculate a 550 station allocation
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
shallow_boat <- 176 # Alaska Provider
deep_boat <- 148    # Ocean Explorer
current_year <- 2025

goa_stn_allocation <- StationAllocationAIGOA::goa_allocate_stations(
  n = 550,
  min_n_per_stratum = 4,
  species = c(
    "arrowtooth flounder", ## Atherestes stomias
    "Pacific cod", ## Gadus macrocephalus
    "walleye pollock", ## Gadus chalcogrammus
    "rex sole", ## Glyptocephalus zachirus
    "flathead sole", ## Hippoglossoides elassodon
    "Pacific halibut", ## Hippoglossus stenolepis
    "southern rock sole", ## Lepidopsetta bilineata
    "northern rock sole", ## Lepidopsetta polyxystra
    "Pacific ocean perch", ## Sebastes alutus
    "silvergray rockfish", ## Sebastes brevispinis
    "northern rockfish", ## Sebastes polyspinis
    "dusky rockfish", ## Sebastes variabilis
    "REBS rockfish", ## Sebastes aleutianus and S. melanostictus
    "Dover sole", ## Microstomus pacificus
    "shortspine thornyhead" ## Sebastolobus alascanus
  ),
  max_iter = 5000,
  trawl = c("Y", "N", "UNK")[c(1, 3)],
  survey_year = current_year,
  planning_years = c(1996, 1999, seq(from = 2003, to = 2023, by = 2))
)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stn_allocation <- goa_stn_allocation$drawn_stations
stn_allocation$ID <- 1:nrow(x = stn_allocation)

## Attach easting and northings to the allocated stations
stn_allocation <- merge(x = stn_allocation,
                        y = goa_stations[, c("STATION", "x", "y")],
                        by = c("STATION"))

stns_in_lane <- terra::relate(x = terra::vect(x = stn_allocation,
                                              geom = c("x", "y")),
                              y = laning_area,
                              relation = "intersects",
                              pairs = T)
strata_in_lane <- unique(stn_allocation[stns_in_lane[, "id.x"],]$STRATUM)
strata_out_lane <- unique(stn_allocation[-stns_in_lane[, "id.x"],]$STRATUM)

stn_allocation[stns_in_lane[, "id.x"], "VESSEL"] <-
  c(deep_boat, shallow_boat)[stns_in_lane[, "id.y"]]

for (istratum in goa_strata$STRATUM) {
  n_in_stratum <- sum(stn_allocation$STRATUM == istratum &
                        is.na(x = stn_allocation$VESSEL))

  vessel_assignments <- vector(length = n_in_stratum)
  vessel_assignments[1:n_in_stratum] <- sample(x = c(deep_boat, shallow_boat))

  stn_allocation[stn_allocation$STRATUM == istratum &
                   is.na(x = stn_allocation$VESSEL), "VESSEL"] <-
    vessel_assignments

}

diff_in_stns <- abs(x = diff(x = table(stn_allocation$VESSEL)))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Rebalance stations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
larger_boat <- names(x = table(stn_allocation$VESSEL))[
  which.max(table(stn_allocation$VESSEL))
]
smaller_boat <- names(x = table(stn_allocation$VESSEL))[
  which.min(table(stn_allocation$VESSEL))
]

set.seed(1 + current_year)
stns_switch <-
  table(sample(x = strata_out_lane,
               prob = goa_stn_allocation$ms_allocation$ms_allocation /
                 sum(goa_stn_allocation$ms_allocation$ms_allocation),
               replace = TRUE, size = diff_in_stns / 2))

for (istratum in as.numeric(names(x = stns_switch)) ) {
  stn_allocation$VESSEL[
    sample(x = with(stn_allocation,
                    which(x = VESSEL == larger_boat
                          & STRATUM == istratum)),
           size = stns_switch[paste(istratum)])
  ] <- as.numeric(smaller_boat)

}
table(stn_allocation$VESSEL)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Plot ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pdf(file = "analysis_scripts/GOA/goa_restratification_2025/goa_stations_allocation.pdf",
    width = 8, height = 6, onefile = TRUE)
for (iarea in c(610, 620, 630, 640, 650)) { ## Loop over area -- start

  ## temporary objects
  temp_strata <- goa_strata[goa_strata$REP_AREA == iarea]
  n_strata <- nrow(x = temp_strata)
  stratum_cols <- c(RColorBrewer::brewer.pal(name = "Set1",
                                             n = n_strata - 1), "cyan")

  ## Plot blank
  plot(temp_strata, border = F, axes = F, col = "white", lwd = 0.1)

  ## add land
  plot(goa_base$akland, add = TRUE, col = "tan", border = T, lwd = 0.1)

  ## add strata
  plot(temp_strata, border = TRUE, axes = F, col = stratum_cols,
       lwd = 0.1, add = TRUE)

  ## Add outlines of stations
  # plot(goa_stations[goa_stations$REP_AREA == iarea],
  #      border = T, lwd = 0.05, add = TRUE )

  ## Add chosen stations
  with(stn_allocation[stn_allocation$REP_AREA == iarea, ],
       points(x, y, pch = c("176" = 1, "148" = 2)[paste(VESSEL)],
              cex = 0.75, lwd = 0.5)
  )

  ## Legend
  legend_labels <- with(as.data.frame(temp_strata),
                        paste0("Stratum ", STRATUM, ": ",
                               DEPTH_MIN_M , " - ", DEPTH_MAX_M, " m"))
  # legend_labels <- c(legend_labels, "Untrawlable")

  legend(c("610" = "topleft", "620" = "topleft", "630" = "bottomright",
           "640" = "bottom", "650" = "bottomleft")[paste(iarea)],
         legend = legend_labels,
         title = "Stratum Legend",
         fill = stratum_cols,
         xpd = NA)

  legend(c("610" = "bottom", "620" = "bottom", "630" = "bottom",
           "640" = "bottomleft", "650" = "bottom")[paste(iarea)],
         legend = c(176, 148), pch = 1:2, ncol = 2, xpd = NA)

  mtext(side = 3, line = -2, text = iarea, font = 2, cex = 1.5)
}  ## Loop over area -- end
dev.off()

## Open pdf file
if (file.exists("analysis_scripts/GOA/goa_restratification_2025/goa_stations_allocation.pdf")) {
  system('open  analysis_scripts/GOA/goa_restratification_2025/goa_stations_allocation.pdf')
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Calculate priority strata for bonus stations and cutting stations
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contingency_table <- goa_stn_allocation$ms_allocation
names(x = contingency_table)[names(x = contingency_table) == "ms_allocation"] <-
  "ms_allocation_550"

for (ieffort in seq(from = 540, to = 450, by = -10)) {
  temp_allocation <- StationAllocationAIGOA::goa_allocate_stations(
    n = ieffort,
    min_n_per_stratum = 4,
    survey_year = current_year,
    planning_years = c(1996, 1999, seq(from = 2003, to = 2023, by = 2))
  )
  contingency_table[, paste0("ms_allocation_", ieffort)] <-
    temp_allocation$ms_allocation$ms_allocation
}

apply(X = contingency_table[, paste0("ms_allocation_",
                                     seq(from = 550, to = 450, by = -10))],
      MARGIN = 2,
      FUN = function(x) tapply(X = x,
                               INDEX = contingency_table$nmfs_area,
                               FUN = sum))
