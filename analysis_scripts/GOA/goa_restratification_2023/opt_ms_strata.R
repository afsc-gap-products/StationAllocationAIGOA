###############################################################################
## Project:       Spatiotemporal Survey Optimization
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Conduct SamplingStrata R package multispecies stratified
##                survey optimization
###############################################################################
rm(list = ls())

##################################################
####  Install a forked version of the SamplingStrata Package from
####  zoyafuso-NOAA's Github page
####
####  Import other required packages
##################################################
library(devtools)
devtools::install_github(repo = "zoyafuso-NOAA/SamplingStrata")
library(SamplingStrata)
library(terra)
library(RColorBrewer)

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
# D_gct <- readRDS("data/GOA/VAST_fit_D_gct.RDS")

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

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Assign grid points to the new strata
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# grid_goa_sp <- terra::vect(x = cbind(ID = 1:nrow(grid_goa), grid_goa),
#                            geom = c("Lon", "Lat"),
#                            crs = lonlat_crs)
# grid_goa_sp <- terra::project(x = grid_goa_sp,
#                               y = updated_goa_strata)
#
# grid_goa_sp <-
#   terra::intersect(x = grid_goa_sp,
#                    y = nmfs)
# grid_goa_sp <- grid_goa_sp[!is.na(grid_goa_sp$area_name),]
# grid_goa_sp <- grid_goa_sp[grid_goa_sp$DEPTH_EFH <= 700, ]
#
# removed_cells <- (1:n_cells)[-grid_goa_sp$ID]
# D_gct <- D_gct[-removed_cells, , ]
# n_cells <- dim(D_gct)[1]

##################################################
####   Collect optimization results from each strata
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

density_input <- D_gct[, spp_list, ]

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
                      MARGIN = 1:2,
                      FUN = sum),
         ncol = length(spp_list),
         dimnames = list(NULL, paste0("Y", 1:length(spp_list)))),

  matrix(data = apply(X = density_input,
                      MARGIN = 1:2,
                      FUN = function(x) sum(x^2)),
         ncol = length(spp_list),
         dimnames = list(NULL, paste0("Y", 1:length(spp_list), "_SQ_SUM")))
)

## Compile Single-Species CVs to establish a lower limit for a given boat
ss_cvs <- vector(length = length(spp_list)); names(ss_cvs) <- spp_list
for (ispp in spp_list){
  result_list <- readRDS( paste0(github_dir, "/Single_Species_Optimization/",
                                 ispp, "/result_list.RDS"))
  ss_cvs[ispp] <- result_list$cvs["actual_cv"]
}

## Initiate CVs to be those calculated under SRS, assign to a variable
## named cv_constraints
## buildStrataDF calculates the stratum means and variances, X1 = 1
##     means to calculate those statics on the whole domain
srs_stats <- SamplingStrata::buildStrataDF(
  dataset = cbind( frame[, -grep(x = names(frame), pattern = "X")],
                   X1 = 1))

srs_n <- as.numeric(550 * table(frame$domainvalue) / nrow(frame))
srs_var <- as.matrix(srs_stats[, paste0("S", 1:length(spp_list))])^2

srs_var <- sweep(x = srs_var,
                 MARGIN = 1,
                 STATS = (1 - srs_n / n_cells) / srs_n,
                 FUN = "*")

srs_cv <- sqrt(srs_var) / srs_stats[, paste0("M", 1:length(spp_list))]

cv_constraints <- srs_cv

## Create CV constraint df
cv <- list()
for (spp in 1:length(spp_list))
  cv[[paste0("CV", spp)]] <- as.numeric(cv_constraints[, spp])
cv[["DOM"]] <- 1:n_dom
cv[["domainvalue"]] <- 1:n_dom
cv <- as.data.frame(cv)

## Set the result directory to which optimization outputs will save
result_dir <- paste0(github_dir, "/Multispecies_Optimization/")
if(!dir.exists(result_dir)) dir.create(path = result_dir, recursive = T)
setwd(result_dir)

#Run optimization
solution <- optimStrata(method = "continuous",
                        errors = cv,
                        framesamp = frame,
                        iter = 300,
                        pops = 100,
                        elitism_rate = 0.1,
                        mut_chance = 1 / (rep(5, n_dom) + 1),
                        nStrata = rep(5, n_dom),
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
####   Tune CV to hit 1, 2 or 3 boats (292, 550, 825 stations)
####   Assume CVs at the gulf-scale regardless on whether scale of the
####   optimization
##################################################
temp_frame <- frame
temp_frame$domainvalue <- 1
srs_stats <- SamplingStrata::buildStrataDF(
  dataset = cbind( temp_frame[, -grep(x = names(temp_frame),
                                      pattern = "X")],
                   X1 = 1))

srs_n <- sum(as.numeric(550 * table(frame$domainvalue) / nrow(frame)))
srs_var <- as.matrix(srs_stats[, paste0("S", 1:length(spp_list))])^2

## SRS statistics
srs_var <- sweep(x = srs_var,
                 MARGIN = 1,
                 STATS = (1 - srs_n / nrow(frame)) / srs_n,
                 FUN = "*")
srs_cv <- sqrt(srs_var) / srs_stats[, paste0("M", 1:length(spp_list))]

error_df <- data.frame("DOM" = "DOM1",
                       srs_cv,
                       "domainvalue"  = 1)
names(error_df)[2:(1 + length(spp_list))] <- paste0("CV", 1:length(spp_list))

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

  updated_cv_constraint <-
    as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "]) * (CV_adj) +
    ss_cvs * (1  - CV_adj)

  error_df[, paste0("CV", 1:length(spp_list))] <-
    as.numeric(updated_cv_constraint)

  temp_bethel <- SamplingStrata::bethel(stratif = temp_stratif,
                                        errors = error_df,
                                        printa = TRUE)

  temp_n <- sum(as.numeric(temp_bethel))

  print(paste0("n = ", temp_n) )
}

##################################################
####   Updated nh
##################################################
sample_allocations <- as.numeric(temp_bethel)
cv_by_boat <-
  data.frame(species = spp_list,
             cv_constraint = as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "]),
             actual_cv =     as.numeric(attributes(temp_bethel)$outcv[, "ACTUAL CV"]))

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
save(list = "result_list", file = "result_list.RData")




