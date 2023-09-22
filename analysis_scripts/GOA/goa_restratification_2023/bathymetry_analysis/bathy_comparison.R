output_dir <- "//AKC0SS-n086/AKC_PubliC/Dropbox/Zimm/GEBCO/GOA/"

library(terra)
ardem <- terra::vect(x = paste0(output_dir,
                                "BTS_test_depths_ARDEM/",
                                "BTS_test_depths_ARDEM.shp"))
ardem$RASTERVALU <- ardem$RASTERVALU * -1

gebco <- terra::vect(x = paste0(output_dir,
                                "BTS_test_depths_GEBCO/",
                                "BTS_test_depths_GEBCO.shp"))
gebco$RASTERVALU <- gebco$RASTERVALU * -1

mix <- terra::vect(x = paste0(output_dir,
                              "BTS_test_depths/",
                              "BTS_test_depths_extract_goa_bathy.shp"))

bathys <- merge(x = subset(x = as.data.frame(ardem),
                           subset = RASTERVALU != 9999,
                           select = -Diffs),
                y = subset(x = as.data.frame(gebco),
                           select = c("hauljoin", "RASTERVALU")) ,
                all = TRUE,
                by = "hauljoin",
                suffixes = c("_ardem", "_gebco"))
bathys <- merge(x = bathys,
                y = subset(x = as.data.frame(mix),
                           select = c("hauljoin", "RASTERVALU")),
                by = "hauljoin")

names(x = bathys) <- c("hauljoin", "obs_depth",
                       "lon", "lat",
                       "ardem", "gebco", "mix")

error <- sweep(x = bathys[, c("ardem", "gebco", "mix")],
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

par(mfrow = c(2, 2), mar = c(3, 3, 2, 1), oma = c(3, 3, 0, 0))
for (isource in c("mix", "gebco", "ardem")) {

  lm_ <- lm(get(isource) ~ obs_depth,
            data = bathys)
  lm_predict <-
    predict(object = lm_,
            newdata = data.frame(obs_depth = 1:1000),
            interval = "predict")

  plot(bathys[, isource] ~ bathys$obs_depth,
       pch = 16, asp = 1, type = "n", ann = F, las = 1)
  polygon(x = c(1:1000, 1000:1),
          y = c(lm_predict[, "lwr"], rev(lm_predict[, "upr"])),
          col = "blue" )
  abline(a = 0, b = 1, lwd = 2, col = "red")

  points(get(isource) ~ obs_depth,
         data = bathys,
         cex = 0.5)

  mtext(side = 3, text = isource, line = 0.5, font = 2)
  text(x = 700,
       y = 200,
       labels = paste0("Mean % Bias: ", round(x = mean(error[, isource],
                                                       na.rm = TRUE),
                                              digits = 1),
                       "\nMedian % Bias: ", round(x = median(error[, isource],
                                                             na.rm = TRUE),
                                                  digits = 1),
                       "\nMAE: ", round(mae[isource], 1),
                       "\nRMSE: ", round(rmse[isource], 1)) )

}
mtext(side = 1, outer = TRUE, "Observed Depth (m)", font = 2)
mtext(side = 2, outer = TRUE, "Raster Value (m)", line = 1, font = 2)

plot(1, type = "n", axes = F)
legend("center", legend = c("identity line", "95% Prediction Interval") ,
       bty = "n", text.col = "white", cex = 1.5)
legend("center", legend = c("identity line", "95% Prediction Interval") ,
       lwd = c(2, NA), bty = "n", col = "red", text.col = "white", cex = 1.5)
legend("center", legend = c("Identity Line", "95% Prediction Interval") ,
       bty = "n", fill = c(NA,"blue"), border = NA, yjust = 0, cex = 1.5)

