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
names(x = pred_grid)[names(x = pred_grid) == "DEPTH_M"] <- "depth_m"

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
  model_args <- list(
    formula = as.formula(cpue_kg_km2 ~ 0 + year_f + s(depth_m, k = 2)),
    data = subset(x = goa_data_geostat,
                  subset = species == species_name),
    mesh = mesh,
    family = do.call(what = ifelse(
      test = species_name %in% c("northern rockfish", "dusky rockfish"),
      yes = "delta_lognormal",
      no = "delta_gamma"),
      args = list(type = "poisson-link")
    ),
    time = "year",
    spatial = "on",
    spatiotemporal = "iid",
    anisotropy = TRUE,
    k_folds = 10,
    parallel = FALSE,
    use_initial_fit = TRUE,
    silent = TRUE
  )

  cat("Running model using depth as a covariate depth with 10-fold CV. ")
  fit_depth_cv <- do.call(what = "sdmTMB_cv", args = model_args)
  cat("NLL:", -fit_depth_cv$sum_loglik, "\n")
  saveRDS(object = fit_depth_cv, file = paste0(result_dir, "/fit_depth_cv.RDS"))

  cat("Running model with no covariates with 10-fold CV. ")
  model_args$formula <- as.formula(cpue_kg_km2 ~ 0 + year_f)
  fit_cv <- do.call(what = "sdmTMB_cv", args = model_args)
  cat("NLL:", -fit_cv$sum_loglik, "\n")
  saveRDS(object = fit_cv, file = paste0(result_dir, "/fit_cv.RDS"))

  better_model <- c("", "depth")[which.min(c(-fit_cv$sum_loglik,
                                             -fit_depth_cv$sum_loglik))]

  ## Refit the model with the lower -sum_loglik across folds but first,
  ## remove the arguments unique to the sdmTMB_cv() function because now we're
  ## using the original sdmTMB() function
  model_args$k_folds <- model_args$parallel <- model_args$use_initial_fit <- NULL
  if (better_model == "depth") model_args$formula <-
    as.formula(cpue_kg_km2 ~ 0 + year_f + s(depth_m, k = 2))
  fit <- do.call(what = "sdmTMB", args = model_args)
  saveRDS(object = fit, file = paste0(result_dir, "/fit.RDS"))

  ## Predict density and create simulated data on the interpolation grid
  cat("Predicting over interpolation grid\n")
  p <- predict(fit, newdata = pred_grid, return_tmb_object = TRUE)
  saveRDS(predict(fit, newdata = pred_grid, type = "response"),
          file = paste0(result_dir, "/predictions.RDS"))

  cat("Simulating densities over interpolation grid\n")
  sim_data <- predict(fit, newdata = pred_grid, type = "response", nsim = 1000)
  saveRDS(sim_data, file = paste0(result_dir, "/sim_data.RDS"))
  rm(sim_data); gc()

  ## Plot predicted density maps and fit diagnostics ----
  # q-q plot
  cat("Running model diagnostics\n")
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
  cat("Finished with", species_name, "\n\n")
}
