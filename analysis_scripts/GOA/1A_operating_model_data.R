##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       GOA Groundfish CPUE Data Synthesis
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Contributors:  Lewis Barnett (lewis.barnett@noaa.gov)
## Description:   Create CPUE dataset used for VAST for species of interest
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Import libraries
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
library(terra)
library(sumfish)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Constants used
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## crs used
lonlat_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
utm_crs <- "+proj=utm +zone=5N +units=km"

## years used
start_year <- 1996
current_year <- 2021

##################################################
#### Import CPUE survey data
#### Import Haul-level data
#### Import Species codes
##################################################
sumfish::getSQL()
GOA <- sumfish::getRacebase(year = c(start_year, current_year),
                            survey = 'GOA')

data <- sumfish::sumHaul(GOA) %>% dplyr::mutate(REGION = "GOA")
species_codes <- GOA$species

##################################################
####   Merge together bathymetry rasters
##################################################
split_bathy <- list()
n_split_rasters <- length(dir("data/GOA/split_goa_bathy_ras/")) / 2
for (i in 1:n_split_rasters) {
  split_bathy[[i]] <- terra::rast(x = paste0("data/GOA/",
                                             "split_goa_bathy_ras",
                                             "/goa_bathy_processed_",
                                             i, ".grd"))
}

bathy <- do.call(what = terra::merge, args = split_bathy)
rm(split_bathy, i)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Join species names
##   Select and rename columns, dropping rows with mising depths
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
species_codes <- dplyr::select(species_codes, -YEAR_ADDED)
data <- dplyr::inner_join(data, species_codes)
data <- data %>% dplyr::select(YEAR, SURVEY = REGION, HAULJOIN = HAULJOIN,
                               SURFACE_TEMPERATURE, GEAR_TEMPERATURE,
                               BOTTOM_DEPTH,
                               EFFORT, WEIGHT,
                               LATITUDE = START_LATITUDE,
                               LONGITUDE = START_LONGITUDE,
                               SPECIES_NAME, COMMON_NAME) %>%
  tidyr::drop_na(BOTTOM_DEPTH, LATITUDE, LONGITUDE)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   Filter to GOA survey, remove tows with 0 bottom depth, and drop 2001,
##   the year when the survey was incomplete and years before 1996 when a
##   different net/soak time was used
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data <- data %>% filter(SURVEY == "GOA",
                        BOTTOM_DEPTH > 0,
                        YEAR != 2001 & YEAR >= 1996)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Sum catches of northern and southern rock sole with rock sole unid.
## (not distinguished until 1996), rename species complex
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rock_soles <- data %>%
  dplyr::filter(COMMON_NAME %in% c("rock sole unid.",
                                   "southern rock sole",
                                   "northern rock sole")) %>%
  dplyr::group_by_at(vars(-WEIGHT, -COMMON_NAME, -SPECIES_NAME)) %>%
  dplyr::summarise(WEIGHT = sum(WEIGHT)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(SPECIES_NAME = "Lepidopsetta spp.",
                COMMON_NAME = "rock soles")

data <- as.data.frame(rbind(data, rock_soles))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Sum catches of blackspooted and rougheye rocks with rougheye and
##  blackspotted rockfish unid.,  rename species complex
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
B_R_rockfishes <- data %>% dplyr::filter(
  COMMON_NAME %in% c("blackspotted rockfish",
                     "rougheye rockfish",
                     "rougheye and blackspotted rockfish unid.")) %>%
  dplyr::group_by_at(vars(-WEIGHT, -COMMON_NAME, -SPECIES_NAME)) %>%
  dplyr::summarise(WEIGHT = sum(WEIGHT)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(SPECIES_NAME = "Sebastes B_R",
                COMMON_NAME = "BS and RE rockfishes")
data <- as.data.frame(rbind(data, B_R_rockfishes))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Sum catches of sculpins
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sculpins <- data %>% dplyr::filter(
  COMMON_NAME %in% c("bigeye sculpin", "great sculpin",
                     "plain sculpin", "yellow Irish lord")) %>%
  dplyr::group_by_at(vars(-WEIGHT, -COMMON_NAME, -SPECIES_NAME)) %>%
  dplyr::summarise(WEIGHT = sum(WEIGHT)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(SPECIES_NAME = "sculpins",
                COMMON_NAME = "sculpins")
data <- as.data.frame(rbind(data, sculpins))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Change name of spiny dogfish to Pacific spiny dogifsh
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data$COMMON_NAME[data$COMMON_NAME == "spiny dogfish"] <-
  "Pacific spiny dogfish"

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Subset data to only the species of interest
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data <- subset(data,
               COMMON_NAME %in% c(
                 ## Species included in the survey optimization
                 "arrowtooth flounder", ## Atherestes stomias
                 "Pacific cod", ## Gadus macrocephalus
                 "walleye pollock", ## Gadus chalcogrammus
                 "rex sole", ## Glyptocephalus zachirus
                 "flathead sole", ## Hippoglossoides elassodon
                 "Pacific halibut", ## Hippoglossus stenolepis
                 "southern rock sole", ## Lepidopsetta bilineata
                 "northern rock sole", ## Lepidopsetta polyxystra
                 "Pacific ocean perch", ## Sebastes alutus
                 "silvergray rockfish", ## Sebastes brevispinis
                 "northern rockfish", ## Sebastes polyspinis
                 "dusky rockfish", ## Sebastes variabilis
                 "BS and RE rockfishes", ## Sebastes aleutianus and S. melanostictus
                 "Dover sole", ## Microstomus pacificus
                 "shortspine thornyhead", ## Sebastolobus alascanus

                 ## Species not included in the survey optimization, but
                 ## included when simulating surveys
                 "sablefish",
                 "Atka mackerel",
                 "shortraker rockfish",
                 "Pacific spiny dogfish",
                 "yelloweye rockfish",
                 "giant octopus",
                 "longnose skate",
                 "big skate",
                 "harlequin rockfish",
                 "giant grenadier",
                 "sculpins"
               ))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Assign station depths from EFH layer
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cpue_shape <- terra::vect(x = data[, c("LONGITUDE", "LATITUDE")],
                          geom = c("LONGITUDE", "LATITUDE"),
                          keepgeom = TRUE,
                          crs = lonlat_crs)
cpue_shape_aea <- terra::project(x = cpue_shape,
                                 y = terra::crs(bathy))
cpue_shape_aea$depth <- terra::extract(x = bathy,
                                       y = cpue_shape_aea)$dblbnd
cpue_shape_aea[, c("COMMON_NAME", "BOTTOM_DEPTH")] <-
  data[, c("COMMON_NAME", "BOTTOM_DEPTH")]

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Plot bathymetry and station locations along with stations without
##   assigned depths
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(bathy)
plot(cpue_shape_aea,
     add = T,
     pch = ".")
plot(cpue_shape_aea[is.na(cpue_shape_aea$depth),],
     add = T,
     pch = 16,
     col = 'red')

mismatched_idx = which(is.na(cpue_shape_aea$depth))
summary(data[mismatched_idx, "BOTTOM_DEPTH"])

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Plot correlation between EFH depths and reported depth from BTS
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(depth ~ BOTTOM_DEPTH,
     data = cpue_shape_aea[cpue_shape_aea$COMMON_NAME == "Pacific cod", ],
     xlab = "Depth recorded by the BTS",
     ylab = "Depth extracted from EFH layer")
cor(cpue_shape_aea$depth, cpue_shape_aea$BOTTOM_DEPTH, use = "complete.obs")
abline(a = 0, b = 1)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Attach depths to dataset, scaled
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data$DEPTH_EFH <- cpue_shape_aea$depth
data$DEPTH_EFH[is.na(data$DEPTH_EFH)] <-
  data$BOTTOM_DEPTH[is.na(data$DEPTH_EFH)]
data$LOG10_DEPTH_EFH <- log10(data$DEPTH_EFH)
data$LOG10_BOTTOM_DEPTH <- log10(data$BOTTOM_DEPTH)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Save
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data <- data[order(data$YEAR, data$SPECIES_NAME),]
write.csv(x = data,
          file = "data/GOA/goa_vast_data_input.csv",
          row.names = F)
