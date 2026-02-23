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
# library(RColorBrewer)
library(ggplot2)

channel <- gapindex::get_connected(check_access = FALSE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Stratum x inpfc combinations and station allocation
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ai_allocation <- ## change file name with new survey year
  sf::st_read(dsn = paste0("Y:/RACE_GF/ALEUTIAN/AI 2026/Station Allocation/",
                           "ai_2026_station_allocation_400stn.gpkg"))

inpfc <- RODBC::sqlQuery(
  channel = channel,
  query = "SELECT STRATUM, AREA_ID AS inpfc FROM GAP_PRODUCTS.STRATUM_GROUPS
           WHERE SURVEY_DEFINITION_ID = 52
           -- INPFC Areas for WAI, SBS, CAI, and EAI
           AND AREA_ID IN (299, 799, 3499, 5699)")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import GOA base layers from akgfmaps package (DESIGN_YEAR 1991)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ai_base_layers <- akgfmaps::get_base_layers(select.region = "ai",
                                            set.crs = "auto",
                                            design.year = 1991)
## Attach inpfc area to the stratum slot
ai_base_layers$survey.strata <- merge(x = ai_base_layers$survey.strata,
                                      y = inpfc,
                                      by = "STRATUM")

# inpfc_colors <- RColorBrewer::brewer.pal(n = 4, name = "Set1")
# names(x = inpfc_colors) <- c(299, 799, 3499, 5699)

# 1. Create a data frame with the longitudes and a range of latitudes
# This ensures the lines follow the "curve" of the projection
divider_lines <- data.frame(
  lon = c(-170, -177, -183),
  lower_lat = c(51, 51, 51),
  upper_lat = c(53.75, 52.75, 53.25),
  label = c("-170", "-177", "-183")
) |>
  # For each longitude, create a line from lat 50 to 56
  # This prevents the line from being a single dot
  (\(d) lapply(split(d, d$lon), function(row) {
    sf::st_linestring(x = matrix(data = c(rep(row$lon, 2),
                                          row$lower_lat, row$upper_lat),
                                 ncol = 2))
  }))() |>
  sf::st_sfc(crs = 4326) |>      # Define as WGS84 (Degrees)
  sf::st_transform(crs = 3338)

place_names <- data.frame(label = c("Western Aleutians", "Central Aleutians",
                                    "Eastern Aleutians", "Southern Bering Sea",
                                    "Stalemate Bank", "Dutch Harbor", "Adak",
                                    "Bering Sea", "Pacific Ocean"),
                          x = c(-1980000, -1700000,
                                -1300000, -890000,
                                -2140000, -850000,  -1540000,
                                -1300000, -2150000),
                          y = c(850000, 680000,
                                520000, 570000,
                                970000, 430000,  417000,
                                750000, 520000),
                          fontface = c("plain", "plain",
                                       "plain", "plain",
                                       "bold", "bold", "bold",
                                       "italic", "italic"),
                          size = c(2, 2, 2, 2, 1.75, 1.75, 1.75, 2, 2) )

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Create Station map
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ggplot() +
  ## Plot land
  geom_sf(data = ai_base_layers$akland, fill = "black", color = NA) +

  # {if (iplot == "strata")
  ## Plot strata color-filled by inpfc area
  # geom_sf(data = ai_base_layers$survey.strata,
  #         fill = inpfc_colors[paste(ai_base_layers$survey.strata$INPFC)],
  #         lwd = 0.1, color = "darkgrey", alpha = 0.25 ) else NULL} +

  ## Plot isobaths
  geom_sf(data = ai_base_layers$bathymetry,
          color = "grey", lwd = 0.1) +

  ## Plot Location of Stations by Vessel
  geom_sf(data = ai_allocation,
          aes(color = factor(VESSEL)),
          size = 0.15) +

  ## Plot INPFC longitudinal line boundaries
  geom_sf(data = divider_lines,
          aes(linetype = "INPFC Boundary"),
          color = "grey30",
          linewidth = 0.5,
          inherit.aes = FALSE) +

  ## Constrain plot to the AI survey area boundaries
  coord_sf(xlim = ai_base_layers$plot.boundary$x,
           ylim = ai_base_layers$plot.boundary$y) +

  ## Format Legend for Stations by Vessel
  scale_color_discrete(name = "",
                       labels = c("Ocean Explorer", "Alaska Provider")) +
  ## Set the line style of the boundaries
  scale_linetype_manual(
    name = NULL,
    values = c("INPFC Boundary" = "dashed")
  ) +
  ## Order the legend so the points are first and then the INPFC boundary
  ## is second (below). Increase the size of the points in the legend
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 1)),
    linetype = guide_legend(order = 2)
  ) +

  ## Plot inpfc area labels
  annotate("text", label = place_names$label,
           x = place_names$x, y = place_names$y,
           fontface = place_names$fontface, size = place_names$size) +

  ## Plot graticules and add x and y axes labels (empty)
  scale_x_continuous(name = "", breaks = c(-170, 180, 170)) +
  scale_y_continuous(name = "", breaks = c(50, 52, 54)) +
  ggspatial::annotation_scale(location = 'bl',
                              line_width = 0.5,
                              text_cex = 0.75) +

  ## Add compass direction
  ggspatial::annotation_north_arrow(
    location = "tr",            # "tl" = top left, "tr" = top right, etc.
    which_north = "true",       # Uses true north based on the projection
    pad_x = unit(0.25, "cm"),    # Distance from the left edge
    pad_y = unit(0.25, "cm"),    # Distance from the top edge
    style = ggspatial::north_arrow_fancy_orienteering(
      fill = c("black", "white"),
      line_col = "grey20"
    ),
    height = unit(0.75, "cm"),
    width = unit(0.75, "cm")
  ) +

  ## Additional plot settings
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.55, 0.925),
    legend.background = element_blank(),
    legend.title.align = 0.5,
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 5),
    legend.key = element_blank(),
    legend.key.width = unit(1, "cm"),
    legend.key.size = unit(0.2, "cm"),
    panel.grid.major = element_line(color = "grey", linewidth = 0.2,
                                    linetype = "dotted"),
    legend.spacing.y = unit(-0.25, "cm"),
    legend.title.align = 0.5,
    legend.box = "vertical"
  )

## Save plot
ggsave(path = "Y:/RACE_GF/ALEUTIAN/AI 2026/Station Allocation/",
       filename = "AI_station_map.png",
       width = 5.5, height = 2.5, units = "in", dpi = 750 )
## Open plot
shell.exec(paste0("Y:/RACE_GF/ALEUTIAN/AI 2026/Station Allocation/",
                  "AI_station_map.png"))
