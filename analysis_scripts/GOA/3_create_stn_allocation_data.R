##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Data synthesis for stratified survey optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Create dataset used for all optimization runs based on a
##                Gulf of Alaska groundfish VAST spatiotemporal
##                operating single-species models
##
##                Set up other constants used in downstream processes
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Library
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_data_geostat <-
  read.csv(file = "data/GOA/sdmtmb_data/goa_data_geostat.csv")
pred_grid <-
  read.csv(file = "data/GOA/sdmtmb_data/goa_2025_interpolation_grid.csv")

years <- paste0("year_", sort(x = unique(x = goa_data_geostat$year)))
n_years <- length(x = years)
n_cells <- nrow(x = pred_grid)
spp_names <- with(read.csv(file = "data/GOA/species_list/species_list.csv"),
                  unique(x = GROUP_CODE))
n_spp <- length(x = spp_names)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Collate predicted densities from univariate sdmTMB model funs
##  D_gct: Predicted density array, cell g, species c, and year t
##  Model results are run in 2_fit_sdmTMB_models.R and saved in temp/
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
D_gct <- array(dim = c(n_cells, n_spp, n_years),
               dimnames = list(NULL, spp_names, years))

for (ispp in spp_names) { ## Loop over species -- start

  ## Load data
  pred_dens_long <- readRDS(file = paste0("temp/", ispp, "/predictions.RDS"))

  ## Widen dataset
  pred_dens_wide <- reshape(pred_dens_long[, c("id", "year", "est")],
                            idvar = "id", timevar = "year",
                            direction = "wide", sep = "_", )
  names(x = pred_dens_wide) <- gsub(x = names(x = pred_dens_wide),
                                    pattern = "est_",
                                    replacement = "year_")
  ## Attach to predicted density array
  D_gct[, ispp, ] <- as.matrix(pred_dens_wide[order(pred_dens_wide$id), years])

}  ## Loop over species -- end

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Save data internally
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
usethis::use_data(D_gct, overwrite = TRUE)
usethis::use_data(pred_grid, overwrite = TRUE)
