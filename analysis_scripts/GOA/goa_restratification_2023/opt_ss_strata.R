###############################################################################
## Project:       Spatiotemporal Survey Optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Conduct SamplingStrata R package multispecies stratified
##                survey optimization
###############################################################################
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Library
##  Install a forked version of the SamplingStrata Package from
##  zoyafuso-NOAA's Github page
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(devtools)
devtools::install_github(repo = "zoyafuso-NOAA/SamplingStrata")
library(SamplingStrata)
library(terra)

##################################################
####   Set up directories based on whether the optimization is being conducted
####        on a multi-species or single-species level
##################################################
github_dir <- "C:/Users/zack.oyafuso/Desktop/optim_res/"

##################################################
####   Load Data
####   Load Population CVs for use in the thresholds
####   Source plotting function
##################################################
source("analysis_scripts/GOA/goa_restratification_2023/plot_solution_results.R")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Load the true density, true index, and spatial domain dataset
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grid_goa <- read.csv(file = "data/GOA/vast_grid_goa.csv")
load("data/D_gct.rda")
load("data/optim_df.rda")

updated_goa_strata <-
  terra::vect(x = "data/GOA/processed_shapefiles/goa_strata_2023.shp")
depth_mods <- read.csv("data/GOA/strata_boundaries/depth_modifications_2023.csv")

nmfs <- terra::vect(x = "data/GOA/shapefiles_from_GDrive/GOA_Shapes.shp")
nmfs <- terra::project(x = nmfs, y = updated_goa_strata)
nmfs$area_name <- c("Southeastern", "Southeastern", "Shumagin", "Chirikof",
                    "Kodiak", "Yakutat", NA)
nmfs <- terra::aggregate(x = nmfs, by = "area_name")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Constants used throughout all scripts
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## crs used
lonlat_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
utm_crs <- "+proj=utm +zone=5N +units=km"
n_years <- dim(D_gct)[3]
n_spp <- dim(D_gct)[2]
n_cells <- dim(D_gct)[1]

##################################################
####   Constants to specify before doing optimization
##################################################

## Domain is the term used in the SamplingStrata package, is used to
## distinguish the five management districts (n_dom == n_districts)
n_dom <- 5

## depth input
# depth_input <- grid_goa$DEPTH_EFH[-removed_cells]
depth_input <- optim_df$DEPTH_EFH

## For the gulf-wide optimization, use 10 strata
## For the district-level optimization, use 5 strata per district
no_strata <-  rep(5, n_dom)

domain_input <-
  as.integer(factor(x = optim_df$INPFC_AREA,
                    levels = c("Shumagin", "Chirikof", "Kodiak",
                               "Yakutat", "Southeastern")))

## Subset strata variables
stratum_var_input <-  data.frame(X1 = depth_input)

spp_list <- c("walleye pollock", "Pacific cod", "arrowtooth flounder",
              "flathead sole", "rex sole", "northern rock sole",
              "southern rock sole", "Dover sole", "Pacific halibut",
              "Pacific ocean perch", "BS and RE rockfishes",
              "silvergray rockfish", "dusky rockfish", "northern rockfish",
              "shortspine thornyhead")

for (which_species in spp_list) { ## Start loop over species

  ## Which density values are we using?
  density_input <- D_gct[, which_species, ]

  ##################################################
  ####   Our df will have fields for:
  ####   domain: only one domain so the value is just 1
  ####   id: unique ID for each sampling cell
  ####   X1: strata variable 2: depth of cell (m)
  ####   X2: strata variable 1: longitude in eastings (km). Because the
  ####       optimization does not read in negative values, I shift the
  ####       values so that the lowest value is 0
  ####
  ####   Variables used to more efficiently calcualte stratum variance
  ####
  ####   WEIGHT: number of observed years
  ####   Y1, Y2, ... : density for a given cell summed across observed years
  ####   Y1_SQ_SUM, Y2_SQ_SUM, ... : density-squared for a given cell,
  ####           summed across observed years
  ##################################################
  frame <- cbind(
    data.frame(domainvalue = domain_input,
               id = (1:n_cells),
               stratum_var_input,
               WEIGHT = n_years),

    matrix(data = apply(X = density_input,
                        MARGIN = 1,
                        FUN = sum),
           ncol = 1,
           dimnames = list(NULL, paste0("Y1"))),

    matrix(data = apply(X = density_input,
                        MARGIN = 1,
                        FUN = function(x) sum(x^2)),
           ncol = 1,
           dimnames = list(NULL, paste0("Y", 1, "_SQ_SUM")))
  )

  ## Set the result directory to which optimization outputs will save
  result_dir = paste0(github_dir, "Single_Species_Optimization/",
                      which_species, "/")
  if(!dir.exists(result_dir)) dir.create(path = result_dir, recursive = T)

  ##################################################
  ####   Run optimization at SRS CV constraints
  ##################################################
  ## Initiate CVs to be those calculated under simple random sampling (SRS)
  srs_stats <- SamplingStrata::buildStrataDF(
    dataset = cbind( frame[, -grep(x = names(frame), pattern = "X")],
                     X1 = 1))

  srs_n <- as.numeric(550 * table(frame$domainvalue) / n_cells)

  ## SRS statistics
  srs_var <- srs_stats$S1^2 * (1 - srs_n / n_cells) / srs_n
  srs_cv <- sqrt(srs_var) / srs_stats$M1

  ## cv is a data input to the SamplingStrata package, assign the initial
  ## cv constraints
  cv <- list()
  cv[["CV1"]] <- srs_cv
  cv[["DOM"]] <- 1:n_dom
  cv[["domainvalue"]] <- 1:n_dom
  cv <- as.data.frame(cv)

  ## Set wd for output files, create a directory if it doesn"t exist yet
  setwd(result_dir)

  #Run optimization, set up a plot layout to show optimization updates
  solution <- optimStrata(method = "continuous",
                          errors = cv,
                          framesamp = frame,
                          iter = 200,
                          pops = 100,
                          elitism_rate = 0.1,
                          mut_chance = 1 / (no_strata[1] + 1),
                          nStrata = no_strata,
                          showPlot = T,
                          writeFiles = T)

  ## Organize result outputs
  solution$aggr_strata$STRATO <- as.integer(solution$aggr_strata$STRATO)
  solution$aggr_strata <-
    solution$aggr_strata[order(solution$aggr_strata$DOM1,
                               solution$aggr_strata$STRATO), ]

  sum_stats <- summaryStrata(solution$framenew,
                             solution$aggr_strata,
                             progress=FALSE)
  sum_stats$stratum_id <- 1:nrow(sum_stats)
  sum_stats$Population <- sum_stats$Population / n_years
  sum_stats$wh <- sum_stats$Allocation / sum_stats$Population
  sum_stats$Wh <- sum_stats$Population / n_cells
  sum_stats <- cbind(sum_stats,
                     subset(x = solution$aggr_strata,
                            select = -c(STRATO, N, COST, CENS, DOM1, X1)))

  plot_solution <- as.factor(paste0(
    "DOM", solution$framenew$DOMAINVALUE,
    " STR", solution$framenew$STRATO))


  plot_solution <- as.integer(plot_solution)

  ##################################################
  ####   Tune CV to hit 550 stations
  ##################################################
  temp_frame <- frame
  temp_frame$domainvalue <- 1
  srs_stats <- SamplingStrata::buildStrataDF(
    dataset = cbind( temp_frame[, -grep(x = names(temp_frame),
                                        pattern = "X")],
                     X1 = 1))
  srs_n <- as.numeric(550 * table(temp_frame$domainvalue) / n_cells)

  ## SRS statistics
  srs_var <- srs_stats$S1^2 * (1 - srs_n / n_cells) / srs_n
  srs_cv <- sqrt(srs_var) / srs_stats$M1

  error_df <- data.frame("DOM" = "DOM1",
                         "CV1" = srs_cv,
                         "domainvalue"  = 1)

  temp_stratif <- solution$aggr_strata
  temp_stratif$N <- temp_stratif$N / n_years
  temp_stratif$DOM1 <- 1

  temp_bethel <- SamplingStrata::bethel(
    errors = error_df,
    stratif = temp_stratif,
    realAllocation = T,
    printa = T)
  temp_n <- sum(ceiling(temp_bethel))

  while (temp_n != 550){
    over_under <- temp_n > 550
    CV_adj <- ifelse(over_under == TRUE,
                     yes = 1.001,
                     no = 0.999)

    error_df$CV1 <-
      as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "]) * CV_adj
    temp_bethel <- SamplingStrata::bethel(stratif = temp_stratif,
                                          errors = error_df,
                                          printa = TRUE)

    temp_n <- sum(as.numeric(temp_bethel))

    print(paste0("n = ", temp_n, ", CV = ",
                 as.numeric(attributes(temp_bethel)$outcv[, "ACTUAL CV"])) )
  }

  ##################################################
  ####   Update sample_allocations with optimal allocation
  ##################################################
  sample_allocations <- sum_stats$Allocation <- as.numeric(temp_bethel)

  cv_by_boat <- with(attributes(temp_bethel),
                     c("cv_constraint" = as.numeric(outcv[, "PLANNED CV "]),
                       "actual_cv" = as.numeric(outcv[, "ACTUAL CV"]),
                       "total_n" = temp_n))

  ##################################################
  ####   Save a plot of the solution
  ##################################################
  plot_solution_results(file_name = "solution.pdf",
                        grid_object =  optim_df,
                        districts_object = nmfs,
                        sol_by_cell = plot_solution,
                        strata_bounds = sum_stats)

  ##################################################
  ####   Save output
  ##################################################
  result_list <- list(solution = solution,
                      sum_stats = sum_stats,
                      cvs = cv_by_boat,
                      sample_allocations = sample_allocations,
                      sol_by_cell = plot_solution)
  saveRDS(obj = result_list, file = "result_list.RDS")

} ## End loop over species


