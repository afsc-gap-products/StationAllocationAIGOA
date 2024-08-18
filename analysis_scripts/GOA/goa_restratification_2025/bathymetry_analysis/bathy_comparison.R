output_dir <- "//AKC0SS-n086/AKC_PubliC/Dropbox/Zimm/GEBCO/GOA/"

# terra::rast("C:/Users/zack.oyafuso/Desktop/goa_bathy/")

library(terra)
ARDEM <- terra::vect(x = paste0(output_dir,
                                "BTS_test_depths_ARDEM3/",
                                "June26_BTS_depths_ARDEM3_NO_NULLS.shp"))
ARDEM$RASTERVALU <- ARDEM$RASTERVALU * -1

GEBCO <- terra::vect(x = paste0(output_dir,
                                "BTS_test_depths_GEBCO3/",
                                "June26_BTS_depths_GEBCO3.shp"))
GEBCO$RASTERVALU <- GEBCO$RASTERVALU * -1

Mix <- terra::vect(x = paste0(output_dir,
                              "BTS_test_depths/",
                              "June26_BTS_depths_extract_goa_bathy.shp"))

bathys <- merge(x = subset(x = as.data.frame(ARDEM),
                           select = c("hauljoin", "c_start_lo", "c_start_la",
                                      "bottom_dep", "wire_lengt", "RASTERVALU")),
                y = subset(x = as.data.frame(GEBCO),
                           select = c("hauljoin", "RASTERVALU")) ,
                all = TRUE,
                by = "hauljoin",
                suffixes = c("_ARDEM", "_GEBCO"))
bathys <- merge(x = bathys,
                y = subset(x = as.data.frame(Mix),
                           select = c("hauljoin", "RASTERVALU")),
                by = "hauljoin")

names(x = bathys) <- c("hauljoin", "lon_dd", "lat_dd",
                       "bottom_depth_m", "wire_length_m",
                       "ARDEM_depth_m", "GEBCO_depth_m", "Mix_depth_m")

error <- sweep(x = bathys[, c("ARDEM_depth_m", "GEBCO_depth_m", "Mix_depth_m")],
               MARGIN = 1,
               STATS = bathys$bottom_depth_m,
               FUN = "-")
rel_error <- 100 * round(sweep(x = error,
                               MARGIN = 1,
                               STATS = bathys$bottom_depth_m,
                               FUN = "/"), 3)
rmse <- sqrt(apply(X = error^2,
                   MARGIN = 2,
                   FUN = mean,
                   na.rm = TRUE))

mae <- apply(X = abs(error),
             MARGIN = 2,
             FUN = mean,
             na.rm = TRUE)

png(filename = paste0("analysis_scripts/GOA/goa_restratification_2025/",
                      "bathymetry_analysis/FigXX_results_8_12_2024.png"),
    width = 190, height = 190, units = "mm", res = 500, family = "serif")

par(mfrow = c(2, 2), mar = c(4, 5, 2, 1), oma = c(0, 0, 0, 0))
hist(bathys$bottom_depth_m, ann = F, col = "darkgrey", nclass = 25, las = 1)
mtext(side = 1, "Observed Depth (m)", font = 2, line = 2.2)
mtext(side = 2, "Frequency", line = 3, font = 2)
mtext(side = 3,
      text = "A) Histogram of Observed Depths",
      line = 0.5, font = 2)
legend(x = 575, y = 750,
       legend = c(
         paste0("Min: ", min(x = bathys$bottom_depth_m, na.rm = T), " m"),
         paste0("Median: ", median(x = bathys$bottom_depth_m, na.rm = T), " m"),
         paste0("Mean: ", round(mean(x = bathys$bottom_depth_m, na.rm = T)), " m"),
         paste0("Max: ", max(x = bathys$bottom_depth_m, na.rm = T), " m")),
       bty = "n")
box()

for (isource in c("Mix", "GEBCO", "ARDEM")) {
  plot(bathys[, paste0(isource, "_depth_m")] ~ bathys$bottom_depth_m,
       xlim = c(0, 1000), ylim = c(0, 1000),
       pch = 16, asp = 1, type = "n", ann = F, las = 1)

  points(get(paste0(isource, "_depth_m")) ~ bottom_depth_m,
         data = bathys,
         cex = 0.5)

  mtext(side = 3,
        text = paste(c("Mix" = "B)", "GEBCO" = "C)", "ARDEM" = "D)")[isource],
                     isource),
        line = 0.5, font = 2)
  legend("bottomright",
         legend = c(paste0("Mean PB: ",
                           round(x = mean(error[, paste0(isource, "_depth_m")],
                                          na.rm = TRUE),
                                 digits = 1), '%'),
                    paste0("Median PB: ",
                           round(x = median(error[, paste0(isource, "_depth_m")],
                                            na.rm = TRUE),
                                 digits = 1), '%'),
                    paste0("MAE: ",
                           round(mae[paste0(isource, "_depth_m")], 1), " m"),
                    paste0("RMSE: ",
                           round(rmse[paste0(isource, "_depth_m")], 1), " m")),
         bty = "n")

  mtext(side = 1, "Observed Depth (m)", font = 2, line = 2.25)
  mtext(side = 2, "Raster Value (m)", line = 3, font = 2)
  abline(a = 0, b = 1, lwd = 2, col = "red")
}

dev.off()

summary(error[bathys$Mix_depth_m != bathys$GEBCO_depth_m, ])
