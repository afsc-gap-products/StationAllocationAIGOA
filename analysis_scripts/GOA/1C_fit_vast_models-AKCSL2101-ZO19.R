##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Single-species VAST model runs
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Lewis Barnett (lewis.barnett@noaa.gov)
##                Jim Thorson"s VAST wiki example
##             (https://github.com/James-Thorson-NOAA/VAST/wiki/Crossvalidation)
## Description:   Run single-species VAST models with and without depth
##                as a smoother with 2 knots. Run 10-fold cross validation for
##                each model
##
##
##  NOTES:        Make sure R and package versions are consistent with those
##                versions speficied in the 2022 Terms of Reference (TOR)
##                used for providing model-based abundance indices for stock
##                assessment.
##                https://docs.google.com/document/d/1t-pIruLZ-F_iNzCysWLH8cdsM0gZpuUb/edit
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Set result directory ----
##   Model output is saved outside of repo due to size of the result outputs
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
VAST_dir <- c("C:/Users/zack.oyafuso/Desktop/VAST_Runs/")
if(!dir.exists(VAST_dir)) dir.create(VAST_dir, recursive = T)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Load packages ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(VAST)
library(effects)
library(splines)
library(units)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Package Version Preferences ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
R_version <- "R version 4.0.2 (2020-06-22)"
current_year <- 2022
VAST_cpp_version <- "VAST_v13_1_0"
pck_version <- c("VAST" = "3.9.0",
                 "FishStatsUtils" = "2.11.0",
                 "Matrix" = "1.4-0",
                 "TMB" = "1.7.22",
                 "DHARMa" = "0.4.5")

{
  if(sessionInfo()$R.version$version.string == R_version)
    message(paste0(sessionInfo()$R.version$version.string,
                   " is consistent with the ", current_year, " TOR."))

  if(!sessionInfo()$R.version$version.string == R_version)
    message(paste0("WARNING: ", sessionInfo()$R.version$version.string,
                   " is NOT consistent with the ", current_year, " TOR. ",
                   "Please update R version to ", R_version))

  for (pck in 1:length(pck_version)) {
    temp_version <- packageVersion(pkg = names(pck_version)[pck])

    if(temp_version == pck_version[pck])
      message(paste0("The version of the '", names(pck_version)[pck],
                     "' package (", temp_version, ") is consistent",
                     " with the ", current_year, " TOR."))

    if(!temp_version == pck_version[pck])
      message(paste0("WARNING: ",
                     "The version of the '", names(pck_version)[pck],
                     "' package (", temp_version, ") is NOT consistent",
                     " with the ", current_year, " TOR. Please update the '",
                     names(pck_version)[pck], "' package to ",
                     pck_version[pck]))
  }
  rm(pck, temp_version)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Data ----
##   Import CPUE data for each species
##   Import goa grid with depths
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
master_data <- read.csv(file = "data/GOA/goa_vast_data_input.csv" )
spp_names <- sort(unique(master_data$COMMON_NAME))

grid_goa <- read.csv("data/GOA/grid_goa.csv")
grid_goa$LOG10_DEPTH <- log10(grid_goa$DEPTH_EFH)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Fit models ----
##   Loop over species to fit models with and without depth covariates
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (ispp in spp_names[22]) { ## Loop over species -- start
  for (depth_in_model in c(F, T)[2]) { ## Loop over covariates -- start

    ## Create directory to store model results
    result_dir <- paste0(VAST_dir, ispp, ifelse(test = depth_in_model,
                                                yes = "_depth/",
                                                no = "/"))

    if (!dir.exists(result_dir)) dir.create(result_dir)

    ## Subset species
    data <- subset(master_data, COMMON_NAME == ispp)

    ## Prepare the dataframe for catch-rate data in the VAST format
    data_geostat <- data.frame(spp = data$COMMON_NAME,
                               HAULJOIN = data$HAULJOIN,
                               Year = data$YEAR,
                               Catch_KG = data$WEIGHT,
                               AreaSwept_km2 = data$EFFORT * 0.01,
                               Lat = data$LATITUDE,
                               Lon = data$LONGITUDE,
                               LOG10_DEPTH = log10(data$BOTTOM_DEPTH),
                               stringsAsFactors = T)

    ##  Assign 10 fold partitions of the data
    n_fold <- 10
    years <- paste0(unique(data_geostat$Year))
    NTime <- length(unique(data_geostat$Year))

    ## split data_geostat df by year, then on each split sub-df, randomly
    ## assign to a fold (1 of n_fold)
    set.seed(current_year)
    foldno <- lapply(
      # Split data_geostat by Year
      X = split.data.frame(data_geostat,
                           f = data_geostat$Year),

      # For each year split, randomly assign fold numbers so that each year is
      # equally split into n_fold folds
      FUN = function(x) {
        fold_no <- sample(x = 1:n_fold,
                          size = length(x$HAULJOIN),
                          replace = T)
        return(split(x$HAULJOIN, fold_no))
      })

    ## Attach fold number to the data_geostat
    for (iyear in years) { ## Loop over years -- start
      for (ifold in paste(1:n_fold)) { ## Loop over folds -- start
        data_geostat[data_geostat$HAULJOIN %in% foldno[[iyear]][[ifold]] ,
                     "fold"] = as.integer(ifold)
      } ## Loop over folds -- end
    } ## Loop over years -- end

    ## VAST Settings
    settings <- FishStatsUtils::make_settings(
      Version = VAST_cpp_version,
      n_x = 500,   # Number of knots
      Region = "User", #User inputted extrapolation grid
      purpose = "index2",
      bias.correct = FALSE,
      fine_scale = TRUE,
      ## Fields for spatial (Omega) or spatiotemporal (epsilon) variation
      ## on the 1st or 2nd predictor
      FieldConfig = c( "Omega1" = "IID", "Epsilon1" = "IID",
                       "Omega2" = "IID", "Epsilon2" = 0),
      "Options" = c("Calculate_Range" = F,
                    "Calculate_effective_area" = F),
      ObsModel = c(2, 1),
      max_cells = Inf,
      use_anisotropy = T)

    ####  Set general arguments to FishStatsUtils::fit_model()
    vast_arguments <- list(

      ## Input settings
      "working_dir" = result_dir,
      "settings" = settings,

      ## Interpolation grid locations and total areas
      "input_grid" = grid_goa[, c("Area_km2", "Lon", "Lat")],

      ## Data inputs
      "Lat_i" = data_geostat[, "Lat"],
      "Lon_i" = data_geostat[, "Lon"],
      "t_i" = data_geostat[, "Year"],
      "c_i" = rep(0, nrow(data_geostat)),
      "b_i" = units::as_units(data_geostat$Catch_KG, "kg"),
      "a_i" = units::as_units(data_geostat$AreaSwept_km2, "km2"),

      ## Output settings
      "getJointPrecision" = TRUE,

      ## Model tuning
      "newtonsteps" = 1,
      "test_fit" = F,

      ## Covariate data: bathy from grids come from the EFH bathymetry layer
      ## and bathy from data_geostat are those collected from the survey
      "covariate_data" = cbind(
        rbind(data_geostat[, c("Lat", "Lon", "LOG10_DEPTH")],
              grid_goa[, c("Lat", "Lon", "LOG10_DEPTH")]),
        Year = NA))

    if (depth_in_model)
      vast_arguments[["X1_formula"]] <- vast_arguments[["X2_formula"]] <-
      ~ splines::bs(LOG10_DEPTH,
                    degree = 2,
                    intercept = FALSE)

    ## Initial Fit
    fit <- do.call(what = FishStatsUtils::fit_model, args = vast_arguments)

    ## Plot depth effects
    if (depth_in_model) {

      # Must add data-frames to global environment (hope to fix in future)
      covariate_data_full = fit$effects$covariate_data_full
      catchability_data_full = fit$effects$catchability_data_full
      depth_effect_pred1 <-
        VAST::Effect.fit_model( fit,
                                focal.predictors = c("LOG10_DEPTH"),
                                which_formula = "X1",
                                xlevels = 100,
                                transformation = list(link=identity,
                                                      inverse=identity) )
      depth_effect_pred2 <-
        VAST::Effect.fit_model( fit,
                                focal.predictors = c("LOG10_DEPTH"),
                                which_formula = "X2",
                                xlevels = 100,
                                transformation = list(link=identity,
                                                      inverse=identity) )

      ## Save a reduced form of the predicted effects
      depth_effect_pred1 <-
        with(depth_effect_pred1,
             data.frame(x = x, pred = fit, lower = lower, upper = upper))

      depth_effect_pred2 <-
        with(depth_effect_pred2,
             data.frame(x = x, pred = fit, lower = lower, upper = upper))

      pdf(file = paste0(result_dir, "depth_effects.pdf"),
          width = 5, height = 5, onefile = TRUE )
      for (ipred in 1:2) {
        temp <- get(paste0("depth_effect_pred", ipred))
        plot(x = temp$LOG10_DEPTH, y = temp$pred, type = "n", xaxt = "n",
             ylim = range(temp[, c("lower", "upper")]),
             las = 1, ylab = "Marginal Effect", xlab = "Depth (m)",
             main = paste0("Depth Effect on 1st Predictor\n", ispp))
        with(temp, polygon(x = c(LOG10_DEPTH, rev(LOG10_DEPTH)),
                           y = c(lower, rev(upper)),
                           col = "cornflowerblue", border = F))
        with(temp, lines(x = LOG10_DEPTH, y = pred, lwd = 2))
        axis(side = 1, labels = NA, tck = -0.015,
             at = log10(c(seq(10, 90, 10), seq(100, 1000, 100))))
        axis(side = 1, labels = NA, tck = -0.03,
             at = log10(c(10,20,50,100,200,500,1000)))
        axis(side = 1, at = log10(c(10,20,50,100,200,500,1000)),
             labels = c(10,20,50,100,200,500,1000))
      }
      dev.off()

      save(list = c("depth_effect_pred1", "depth_effect_pred2"),
           file = paste0(result_dir, "depth_effect.RData"))
    }

    ##  Diagnostics plots
    if(!dir.exists(paste0(result_dir, "/diagnostics")))
      dir.create(paste0(result_dir, "/diagnostics"))
    plot(x = fit, working_dir = paste0(result_dir, "diagnostics/"))

    ## Save original model fit
    save(list = "fit", file = paste0(result_dir, "/fit.RDS"))

    ##################################################
    ####   10-fold Cross Validation
    ####   First, save MLE parameter values to the VAST arguments to be used
    ####   as starting values for the subsequent CV runs
    ##################################################
    vast_arguments$Parameters <- fit$ParHat

    ## Loop through partitions, refitting each time with a different PredTF_i
    for (fI in 1:n_fold) { ## Loop over folds -- start

      ## Create directory for CV run
      if (!dir.exists(paste0(result_dir, "CV_", fI))){
        dir.create(paste0(result_dir, "CV_", fI))
        file.copy(from = paste0(result_dir, "Kmeans_knots-",
                                settings$n_x, ".RData"),
                  to = paste0(result_dir, "CV_", fI, "/"))
      }

      ## Update which indices are withheld for prediction
      PredTF_i <- ifelse( test = data_geostat$fold == fI,
                          yes = TRUE,
                          no = FALSE )
      vast_arguments$PredTF_i <- PredTF_i
      vast_arguments$working_dir <- paste0(result_dir, "CV_", fI, "/")

      ## Fit CV run
      fit_CV <- do.call(what = FishStatsUtils::fit_model, args = vast_arguments)

      ## Save predicted and observed CPUEs
      obs_cpue <- with(data_geostat[PredTF_i, ], Catch_KG / AreaSwept_km2)
      pred_cpue <- fit_CV$Report$D_i[PredTF_i]

      cv_performance <-
        list(cpues = data.frame(cv_fold = fI,
                                obs_cpue,
                                pred_cpue),
             prednll = fit_CV$Report$pred_jnll,
             max_gradient = fit_CV$parameter_estimates$max_gradient)

      ## Save fit
      save(list = "fit_CV",
           file = paste0(result_dir, "CV_", fI, "/fit.RData"))
      save(cv_performance,
           file = paste0(result_dir, "CV_", fI,
                         "/crossval_fit_performance.RData"))

      if (depth_in_model) {
        # Must add data-frames to global environment (hope to fix in future)
        covariate_data_full = fit$effects$covariate_data_full
        catchability_data_full = fit$effects$catchability_data_full
        depth_effect_pred1 <-
          VAST::Effect.fit_model( fit,
                                  focal.predictors = c("LOG10_DEPTH"),
                                  which_formula = "X1",
                                  xlevels = 100,
                                  transformation = list(link=identity,
                                                        inverse=identity) )
        depth_effect_pred2 <-
          VAST::Effect.fit_model( fit,
                                  focal.predictors = c("LOG10_DEPTH"),
                                  which_formula = "X2",
                                  xlevels = 100,
                                  transformation = list(link=identity,
                                                        inverse=identity) )

        ## Save a reduced form of the predicted effects
        depth_effect_pred1 <-
          with(depth_effect_pred1,
               data.frame(x = x, pred = fit, lower = lower, upper = upper))

        depth_effect_pred2 <-
          with(depth_effect_pred2,
               data.frame(x = x, pred = fit, lower = lower, upper = upper))

        save(list = c("depth_effect_pred1", "depth_effect_pred2"),
             file = paste0(result_dir, "CV_", fI, "/depth_effect.RData"))
      }

    } ## Loop over folds -- end
  }  ## Loop over covariates -- end
} ## Loop over species -- end
