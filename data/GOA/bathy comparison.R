dir_out <- ("//AKC0SS-n086/AKC_PubliC/Dropbox/Zimm/GEBCO/GOA/BTS_test_depths/")

library(terra)

goa <- terra::vect(paste0(dir_out, "BTS_test_depths_extract_goa_bathy.shp"))
goa$perc_diff <- goa$Diffs / goa$bottom_dep

hist(goa$perc_diff, nclass = 100); abline(v = 0, lwd = 2)
hist(goa$Diffs, nclass = 100); abline(v = 0, lwd = 2)

summary(goa$Diffs)

plot(goa, pch = 16, col = rev(terrain.colors(n = 100))[abs(100 * goa$perc_diff)], cex = 0.4)
points(goa[abs(goa$perc_diff) > .10, ], col = "red", pch = 16)

