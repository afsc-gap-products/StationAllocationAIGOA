##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Calculate MS allocation of STRS
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Gulf of Alaska 2023 bottom trawl survey
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import Library
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(devtools)
devtools::install_github(repo = "zoyafuso-NOAA/SamplingStrata")
library(SamplingStrata)
library(rgdal)
library(raster)
library(sp)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Load data
##   Only use predicted densities for the species included in optimization
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("data/GOA/optimization_data.RData")
grid_goa <- read.csv(file = "data/GOA/grid_goa.csv")
D_gct <- readRDS("data/GOA/VAST_fit_D_gct.RDS")[, spp_idx_opt, ]

## Think about ways to put this in package!!
updated_goa_strata <- rgdal::readOGR(dsn = "C:/Users/zack.oyafuso/Work/GitHub/Optimal_Allocation_GoA/products/updated_goa_strata/updated_goa_strata.shp")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Assign grid points to the new strata
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grid_goa_sp <- sp::SpatialPoints(coords = grid_goa[, c("Lon", "Lat")],
                                 proj4string = sp::CRS(lonlat_crs))
grid_goa_sp <- sp::spTransform(x = grid_goa_sp,
                               CRSobj = crs(updated_goa_strata))

grid_goa_sp <- raster::intersect(x = grid_goa_sp, y = updated_goa_strata)

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
  data.frame(domainvalue = 1,
             id = 1:n_cells,
             WEIGHT = n_years),

  matrix(data = apply(X = D_gct,
                      MARGIN = 1:2,
                      FUN = sum),
         ncol = ns_opt,
         dimnames = list(NULL, paste0("Y", 1:ns_opt))),

  matrix(data = apply(X = D_gct,
                      MARGIN = 1:2,
                      FUN = function(x) sum(x^2)),
         ncol = ns_opt,
         dimnames = list(NULL, paste0("Y", 1:ns_opt, "_SQ_SUM")))
)

## Initiate CVs to be those calculated under SRS, assign to a variable
## named cv_constraints
## buildStrataDF calculates the stratum means and variances, X1 = 1
##     means to calculate those statics on the whole domain
srs_stats <- SamplingStrata::buildStrataDF( dataset = cbind( frame, X1 = 1))
srs_n <- 550
srs_var <- as.matrix(srs_stats[, paste0("S", 1:ns_opt)])^2

srs_var <- sweep(x = srs_var,
                 MARGIN = 1,
                 STATS = (1 - srs_n / n_cells) / srs_n,
                 FUN = "*")

srs_cv <- sqrt(srs_var) / srs_stats[, paste0("M", 1:ns_opt)]
