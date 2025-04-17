##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create Stratum and Station Map for 2025 GOA Survey
##   Zack Oyafuso (zack.oyafuso@noaa.gov)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import libraries and connect to Oracle
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(akgfmaps)
library(gapindex)
library(RColorBrewer)
library(ggplot2)

channel <- gapindex::get_connected(check_access = FALSE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import StratumxNMFS combinations and station allocation
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nmfs <- RODBC::sqlQuery(
  channel = channel,
  query = "SELECT STRATUM, AREA_ID AS NMFS FROM GAP_PRODUCTS.STRATUM_GROUPS
           WHERE SURVEY_DEFINITION_ID = 47
           AND AREA_ID IN (610, 620, 630, 640, 650)")

goa_2025_allocation <-
  sf::st_read(dsn = "G:/GOA/GOA 2025/Station Allocation/goa_2025_station_allocation_450_aea.gpkg")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import GOA base layers from akgfmaps package
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_base_layers <- akgfmaps::get_base_layers(select.region = "goa",
                                             set.crs = "auto",
                                             design.year = 2025)
## Attach NMFS area to the stratum slot
goa_base_layers$survey.strata <- merge(x = goa_base_layers$survey.strata,
                                       y = nmfs,
                                       by = "STRATUM")

nmfs_colors <- RColorBrewer::brewer.pal(n = 5, name = "Set1")
names(x = nmfs_colors) <- c(610, 620, 630, 640, 650)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create Stratum Map
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ggplot() +
  ## Plot land
  geom_sf(data = goa_base_layers$akland, fill = "tan", color = NA) +
  ## Plot strata color-filled by NMFS area
  geom_sf(data = goa_base_layers$survey.strata,
          fill = nmfs_colors[paste(goa_base_layers$survey.strata$NMFS)],
          lwd = 0.1, color = "black", alpha = 0.5  ) +
  ## Plot graticules
  geom_sf(data = goa_base_layers$graticule,
          color = "black", lwd = 0.05, cex = 1.5) +
  ## Constrain plot to the GOA survey area boundaries
  coord_sf(xlim = goa_base_layers$plot.boundary$x,
           ylim = goa_base_layers$plot.boundary$y) +
  ## Plot NMFS area labels
  annotate("text",
           x = c(-653732, -103732, 270004, 600004, 980004),
           y = c(400000, 478000, 748000, 970000, 770000),
           cex = 2.5,
           label = c("Shumagin (610)",
                     "Chirikof\n(620)",
                     "Kodiak\n(630)",
                     "West Yakutat\n(640)",
                     "Southeast\nOutside\n(650)") ) +
  ## Add x and y axes labels (empty)
  scale_x_continuous(name = "",
                     breaks = goa_base_layers$lon.breaks) +
  scale_y_continuous(name = "",
                     breaks = goa_base_layers$lat.breaks) +
  ggspatial::annotation_scale(location = 'tl', line_width = 0.5) +
  theme_bw()
## Save
ggsave(path = "analysis_scripts/GOA/goa_restratification_2025/",
       filename = "stratum_map.png",
       width = 6, height = 3, units = "in", dpi = 750 )

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create Station Map
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ggplot() +
  ## Plot land
  geom_sf(data = goa_base_layers$akland, fill = "tan", color = NA) +
  ## Plot strata color-filled by NMFS area
  geom_sf(data = goa_base_layers$survey.strata,
          fill = nmfs_colors[paste(goa_base_layers$survey.strata$NMFS)],
          lwd = 0.1, color = "grey50", alpha = 0.5 ) +
  geom_sf(data = goa_2025_allocation,
          aes(shape = factor(VESSEL)),  ## Map shape to vessel type
          cex = 0.75) +
  ## Define manual shape scale for legend
  scale_shape_manual(name = "Charter Vessel",
                     values = c("148" = 4, "176" = 16),
                     labels = c("Ocean Explorer", "Alaska Provider")) +
  ## Plot graticules
  geom_sf(data = goa_base_layers$graticule,
          color = "black", lwd = 0.05, cex = 1.5) +
  ## Constrain plot to the GOA survey area boundaries
  coord_sf(xlim = goa_base_layers$plot.boundary$x,
           ylim = goa_base_layers$plot.boundary$y) +
  ## Plot NMFS area labels
  annotate("text",
           x = c(-653732, -103732, 270000, 600000, 980000),
           y = c(400000, 478000, 748000, 970000, 770000),
           cex = 2.5,
           label = c("Shumagin (610)",
                     "Chirikof\n(620)",
                     "Kodiak\n(630)",
                     "West Yakutat\n(640)",
                     "Southeast\nOutside\n(650)") ) +
  ## Add x and y axes labels (empty)
  scale_x_continuous(name = "",
                     breaks = goa_base_layers$lon.breaks) +
  scale_y_continuous(name = "",
                     breaks = goa_base_layers$lat.breaks) +
  ggspatial::annotation_scale(location = 'tl', line_width = 0.5) +
  theme_bw() +
  theme(
    legend.position = c(0.675, 0.25), ## Move legend to bottom-right
    legend.text = element_text(size = 5), ## Reduce text size
    legend.title = element_text(size = 5) ## Reduce title size
  )
ggsave(filename = paste0("analysis_scripts/GOA/goa_restratification_2025/",
                         "goa_2025_station_map.png"),
       width = 6, height = 3, units = "in", dpi = 500)
