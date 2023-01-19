library(terra)
library(RColorBrewer)

goa_strata <- terra::vect(x = "data/GOA/shapefiles_from_GDrive/goa_strata.shp")
goa_strata <- goa_strata[goa_strata$STRATUM != 0, ]
updated_goa_strata <- terra::vect(x = "data/GOA/processed_shapefiles/goa_strata_2023.shp")

nmfs_areas <- terra::vect(x = "data/GOA/shapefiles_from_GDrive/GOA_Shapes.shp")
nmfs_areas <- terra::project(x = nmfs_areas, goa_strata)

png(filename = "presentations/122022 SA/historical_subareas.png",
    width = 5, height = 3, units = "in", res = 500)

par(mfrow = c(1, 1))
plot(nmfs_areas[nmfs_areas$REP_AREA %in% c(610, 620, 630, 640, 650)],
     border = F, axes = F, mar = c(0,0,0,0))
plot(goa_strata,
     col = c(RColorBrewer::brewer.pal(n = 5, name = "Set1"))[
       with(as.data.frame(goa_strata),
            floor((STRATUM - floor(STRATUM/100)*100)/10))],
     add = TRUE, border = F)
plot(nmfs_areas[nmfs_areas$REP_AREA %in% c(610, 620, 630, 640, 650)], add = TRUE)
text(c(-728920.1, -155202.0,  223950.7,  603103.5, 1047111.4),
     c(296550.8, 406305.5, 635792.8, 930135.0, 730580.9),
     c(610, 620, 630, 640, 650))

legend("topleft", legend = c("Shumagin", "Chirikof", "Kodiak",
                             "Yakutat", "Southeastern"),
       fill = RColorBrewer::brewer.pal(n = 5, name = "Set1"),
       border = F,
       bty = "n")

dev.off()

png(filename = "presentations/122022 SA/updated_subareas.png",
    width = 5, height = 3, units = "in", res = 500)

plot(nmfs_areas[nmfs_areas$REP_AREA %in% c(610, 620, 630, 640, 650)],
     border = F, axes = F, mar = c(0,0,0,0))
plot(updated_goa_strata,
     col = c(RColorBrewer::brewer.pal(n = 5, name = "Set1"))[
       sapply(X = updated_goa_strata$INPFC_AREA,
              FUN = function(x) switch(x,
                                       "Shumagin" = 1, "Chirikof" = 2, "Kodiak" = 3,
                                       "Yakutat" = 4, "Southeastern" = 5))
       ], border = F, add = TRUE)

plot(nmfs_areas[nmfs_areas$REP_AREA %in% c(610, 620, 630, 640, 650)],
     add = TRUE)
text(c(-728920.1, -155202.0,  223950.7,  603103.5, 1047111.4),
     c(296550.8, 406305.5, 635792.8, 930135.0, 730580.9),
     c(610, 620, 630, 640, 650))
legend("topleft", legend = c(610, 620, 630, 640, 650),
       fill = RColorBrewer::brewer.pal(n = 5, name = "Set1"),
       border = F,
       bty = "n")
dev.off()
