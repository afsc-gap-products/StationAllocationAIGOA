##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Synthesize the Cross Validation
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Loop through every CV fold and extract performance metrics
##                Maximum Gradient, RRMSE, Predictive Joint Negative LogLike.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Restart R Session before running
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Import Libraries ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyr)
library(VAST)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Set constants ----
##   Set VAST_dir to the external drive that the VAST runs are stored
##   Import the VAST data input
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
VAST_dir <- "G:/Oyafuso/VAST_Runs_GOA_2023/"
data_df <- read.csv("data/GOA/goa_vast_data_input.csv")

##################################################
####   Result object
##################################################
cross_val_results <- data.frame()

##################################################
####   Loop through each model (species, depth, CV fold) and calculate RRMSE
####   synthesize RRMSE, pred_jnll, maximum gradient
##################################################

spp_names <- sort(unique(data_df$COMMON_NAME))
for (ispp in spp_names){ ## Loop through species -- start

  for (depth_in_model in c("",
                           "_depth")) { ## Loop through depth models -- start
    for (icv in 1:10) { ## Loop through cv folds -- start

      ## Only go through CV folds where a model run was successful
      fit_file <- paste0(VAST_dir, ispp, depth_in_model,
                         "/CV_", icv, "/", "crossval_fit_performance.RDS")

      if(file.exists(fit_file)) {

        ## Load performance metrics
        cv_performance <- readRDS(paste0(VAST_dir, ispp, depth_in_model,
                    "/CV_", icv, "/", "crossval_fit_performance.RDS"))
        load(paste0(VAST_dir, ispp, depth_in_model,
                    "/CV_", icv, "/", "parameter_estimates.RData"))

        ## Record results
        cross_val_results <-
          rbind(cross_val_results,
                data.frame(spp_name = ispp,
                           depth_in_model = ifelse(depth_in_model == "", F, T),
                           cv_fold = icv,
                           max_gradient = parameter_estimates$max_gradient,
                           pred_jnll = cv_performance$prednll ))

      }  ## Loop through depth models -- start
    } ## Loop through cv folds -- end
  }

} ## Loop through species --end

##################################################
####   Synthesize Cross Validation Results
##################################################
converged <-
  tidyr::spread(data = aggregate(max_gradient ~ depth_in_model + spp_name,
                                 FUN = function(x) sum(x < 1e4),
                                 data = cross_val_results),
                key = depth_in_model,
                value = max_gradient)

pred_jnll <-
  tidyr::spread(data = aggregate(pred_jnll ~ depth_in_model + spp_name,
                                 FUN = function(x) round(mean(x, na.rm = T)),
                                 data = cross_val_results,
                                 subset = max_gradient < 1e-4),
                key = depth_in_model,
                value = pred_jnll)

##################################################
####   Synthesize the densities and abundance indices for the best models
##################################################
fit <- readRDS(paste0(VAST_dir, "arrowtooth flounder/fit.RDS"))
year_idx <- 1 + as.integer(names(table(fit$data_list$t_i)))
n_years <- length(year_idx)
n_cells <- dim(fit$Report$D_gct)[1]
n_spp <- nrow(pred_jnll)
spp_names <- pred_jnll$spp_name

D_gct <- array(dim = c(n_cells, n_spp, n_years),
                        dimnames = list(NULL, spp_names, NULL))

for(irow in 1:nrow(pred_jnll)) { ## Loop over species -- start

  ## Extract file name of best model
  ispp <- pred_jnll$spp_name[irow]
  depth_in_model <- c(FALSE, TRUE)[which.min(pred_jnll[irow, 2:3])]
  if (irow == 19) depth_in_model <- FALSE

  filename <- paste0(VAST_dir, ispp, ifelse(test = depth_in_model,
                                            yes = "_depth/",
                                            no = "/"), "fit.RDS")



  ## Load data, extract predicted density
  fit <- readRDS(filename)
  D_gct[, irow, ] <- fit$Report$D_gct[, 1, year_idx]

  print(paste(ispp, ifelse(depth_in_model, "with Depth", "without Depth")))
}  ## Loop over species -- end

##################################################
####   Save
##################################################
saveRDS(object = D_gct, file = "data/GOA/VAST_fit_D_gct.RDS")
save(list = c("converged", "cross_val_results", "pred_jnll"),
     file = "data/GOA/prednll_VAST_models.RData")
