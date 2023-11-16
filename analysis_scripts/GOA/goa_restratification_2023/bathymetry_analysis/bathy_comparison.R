output_dir <- "//AKC0SS-n086/AKC_PubliC/Dropbox/Zimm/GEBCO/GOA/"

# terra::rast("C:/Users/zack.oyafuso/Desktop/goa_bathy/")

library(terra)
ARDEM <- terra::vect(x = paste0(output_dir,
                                "BTS_test_depths_ARDEM/",
                                "BTS_test_depths_ARDEM.shp"))
ARDEM$RASTERVALU <- ARDEM$RASTERVALU * -1

GEBCO <- terra::vect(x = paste0(output_dir,
                                "BTS_test_depths_GEBCO/",
                                "BTS_test_depths_GEBCO.shp"))
GEBCO$RASTERVALU <- GEBCO$RASTERVALU * -1

Mix <- terra::vect(x = paste0(output_dir,
                              "BTS_test_depths/",
                              "BTS_test_depths_extract_goa_bathy.shp"))

bathys <- merge(x = subset(x = as.data.frame(ARDEM),
                           subset = RASTERVALU != 9999,
                           select = -Diffs),
                y = subset(x = as.data.frame(GEBCO),
                           select = c("hauljoin", "RASTERVALU")) ,
                all = TRUE,
                by = "hauljoin",
                suffixes = c("_ARDEM", "_GEBCO"))
bathys <- merge(x = bathys,
                y = subset(x = as.data.frame(Mix),
                           select = c("hauljoin", "RASTERVALU")),
                by = "hauljoin")

names(x = bathys) <- c("hauljoin", "obs_depth",
                       "lon", "lat",
                       "ARDEM", "GEBCO", "Mix")

error <- sweep(x = bathys[, c("ARDEM", "GEBCO", "Mix")],
               MARGIN = 1,
               STATS = bathys$obs_depth,
               FUN = "-")
rel_error <- 100 * round(sweep(x = error,
                               MARGIN = 1,
                               STATS = bathys$obs_depth,
                               FUN = "/"), 3)
rmse <- sqrt(apply(X = error^2,
                   MARGIN = 2,
                   FUN = mean,
                   na.rm = TRUE))

mae <- apply(X = abs(error),
             MARGIN = 2,
             FUN = mean,
             na.rm = TRUE)

png(filename = paste0("analysis_scripts/GOA/goa_restratification_2023/",
                      "bathymetry_analysis/FigXX_results.png"),
    width = 190, height = 190, units = "mm", res = 500, family = "serif")

par(mfrow = c(2, 2), mar = c(4, 5, 2, 1), oma = c(0, 0, 0, 0))
hist(bathys$obs_depth, ann = F, col = "darkgrey", nclass = 25, las = 1)
mtext(side = 1, "Observed Depth (m)", font = 2, line = 2.2)
mtext(side = 2, "Frequency", line = 3, font = 2)
mtext(side = 3,
      text = "A) Histogram of Observed Depths",
      line = 0.5, font = 2)
legend(x = 575, y = 750,
       legend = c(
         paste0("Min: ", min(x = bathys$obs_depth, na.rm = T), " m"),
         paste0("Median: ", median(x = bathys$obs_depth, na.rm = T), " m"),
         paste0("Mean: ", round(mean(x = bathys$obs_depth, na.rm = T)), " m"),
         paste0("Max: ", max(x = bathys$obs_depth, na.rm = T), " m")),
       bty = "n")
box()

for (isource in c("Mix", "GEBCO", "ARDEM")) {
  plot(bathys[, isource] ~ bathys$obs_depth,
       xlim = c(0, 1000), ylim = c(0, 1000),
       pch = 16, asp = 1, type = "n", ann = F, las = 1)

  points(get(isource) ~ obs_depth,
         data = bathys,
         cex = 0.5)

  mtext(side = 3,
        text = paste(c("Mix" = "B)", "GEBCO" = "C)", "ARDEM" = "D)")[isource], isource),
        line = 0.5, font = 2)
  legend("bottomright",
         legend = c(paste0("Mean PB: ", round(x = mean(error[, isource],
                                                           na.rm = TRUE),
                                                  digits = 1)),
                    paste0("Median PB: ", round(x = median(error[, isource],
                                                               na.rm = TRUE),
                                                    digits = 1)),
                    paste0("MAE: ", round(mae[isource], 1)),
                    paste0("RMSE: ", round(rmse[isource], 1))),
         bty = "n")

  mtext(side = 1, "Observed Depth (m)", font = 2, line = 2.25)
  mtext(side = 2, "Raster Value (m)", line = 3, font = 2)
  abline(a = 0, b = 1, lwd = 2, col = "red")
}

dev.off()
