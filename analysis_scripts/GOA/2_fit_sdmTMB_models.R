# Run all GOA sdmTMB spatiotemporal models ----
# Author: Lewis Barnett, Zack Oyafuso

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import packages
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(sdmTMB)
library(dplyr)
library(ggplot2)
library(here)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import data and interpolation grid
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
goa_data_geostat <-
  read.csv(file = "data/GOA/sdmtmb_data/goa_data_geostat.csv")
## Year is treated as a factor in sdmTMB
goa_data_geostat$year_f <- as.factor(x = goa_data_geostat$year)

species_list <-
  with(read.csv(file = "data/GOA/species_list/species_list.csv"),
       unique(x = GROUP_CODE))
year_set <- sort(x = unique(x = goa_data_geostat$year))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Format mesh and interpolation grid
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load grid and process prediction grid for all years desired
pred_grid <-
  read.csv(file = "data/GOA/sdmtmb_data/goa_2025_interpolation_grid.csv")
pred_grid <- cbind(id = 1:nrow(x = pred_grid),
                   pred_grid) |>
  sdmTMB::replicate_df(time_name = "year_f",
                       time_values = year_set)
pred_grid$year_f <- as.factor(x = pred_grid$year_f)
pred_grid$year <- as.integer(as.character(factor(pred_grid$year_f)))

mesh <- sdmTMB::make_mesh(subset(x = goa_data_geostat,
                                 ## there are 15 species, so just filter one
                                 subset = species == species_list[1],
                                 select = c("X", "Y")),
                          xy_cols = c("X", "Y"),
                          n_knots = 500)

for (species_name in species_list) {

  cat("Working on", species_name, "\n")

  ## Set up temporary directory where results will be saved
  result_dir <- paste0(getwd(), "/temp/", species_name)
  if (!dir.exists(paths = result_dir))
    dir.create(path = result_dir, recursive = TRUE)

  ## Fit model
  cat("Currently fitting model...")
  fit <- sdmTMB::sdmTMB(
    formula = cpue_kg_km2 ~ 0 + year_f,
    data = subset(x = goa_data_geostat,
                  subset = species == species_name),
    mesh = mesh,
    family = do.call(what = ifelse(
      test = species_name %in% c("northern rockfish", "dusky rockfish"),
      yes = "delta_lognormal",
      no = "delta_gamma"
    ),
    args = list(type = "poisson-link")
    ),
    time = "year",
    spatial = "on",
    spatiotemporal = "iid",
    anisotropy = TRUE,
    silent = TRUE
  )
  saveRDS(fit, file = paste0(result_dir, "/fit.RDS"))
  cat("Complete\n")

  ## Predict density and create simulated data on the interpolation grid
  cat("Predicting over interpolation grid...")
  p <- predict(fit, newdata = pred_grid, return_tmb_object = TRUE)
  saveRDS(predict(fit, newdata = pred_grid, type = "response"),
          file = paste0(result_dir, "/predictions.RDS"))
  cat("Complete\n")

  cat("Simulating densities over interpolation grid...")
  sim_data <- predict(fit, newdata = pred_grid, type = "response", nsim = 1000)
  saveRDS(sim_data, file = paste0(result_dir, "/sim_data.RDS"))
  rm(sim_data); gc()
  cat("Complete\n")

  ## Plot predicted density maps and fit diagnostics ----
  # q-q plot
  cat("Running model diagnostics...")
  pdf(file = paste0(result_dir, "/qq.pdf"), width = 5, height = 5)
  sims <- simulate(object = fit, nsim = 500, type = "mle-mvn")
  sims |> dharma_residuals(object = fit)
  dev.off()

  # residuals on map plot, by year
  resids <- sims |>
    dharma_residuals(fit, test_uniformity = FALSE, return_DHARMa = TRUE)
  fit$data$resids <- resids$scaledResiduals

  ggplot(subset(x = fit$data,
                subset = !is.na(x = resids) & is.finite(x = resids)),
         aes(X, Y, col = resids)) +
    scale_colour_gradient2(name = "residuals", midpoint = 0.5) +
    geom_point(size = 0.7) +
    facet_wrap(~year, ncol = 2) +
    coord_fixed() +
    theme_bw()
  ggsave(file = paste0(result_dir, "/residuals_map.pdf"),
         height = 9, width = 6.5, units = "in")

  ggplot(p$data, aes(X, Y, fill = exp(est1 + est2))) +
    geom_tile() +
    scale_fill_viridis_c(trans = "sqrt", name = "") +
    facet_wrap(~year, ncol = 2) +
    coord_fixed() +
    ggtitle("Predicted densities (kg / square km)") +
    theme_bw()
  ggsave(file = paste0(result_dir, "/predictions_map.pdf"),
         height = 9, width = 6.5, units = "in")

  cat("Complete\n")
  cat("Finished with", species_name, "\n\n")
}
