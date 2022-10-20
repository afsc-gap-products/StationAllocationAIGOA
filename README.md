# **Station Allocation Aleutian Islands and Gulf of Alaska**

This package is intended to be used to document and reproduce the station
allocation for the Aleutian Islands (even years) and Gulf of Alaska (odd years)
bottom trawl surveys. 

## Install
```
library(devtools)
devtools::install_github("zoyafuso-NOAA/SamplingStrata")
devtools::install_github("afsc-gap-products/StationAllocationAIGOA")
```

#### *Station Allocation: Gulf of Alaska 2023 Bottom Trawl Survey*
```
spp_list <- c("walleye pollock", "Pacific cod", "arrowtooth flounder", 
              "flathead sole", "rex sole", "northern rock sole", 
              "southern rock sole", "Dover sole", "Pacific halibut", 
              "Pacific ocean perch", "BS and RE rockfishes", 
              "silvergray rockfish", "dusky rockfish", "northern rockfish",
              "shortspine thornyhead")

GOA_allocation_2023 <- 
      StationAllocationAIGOA::goa_allocate_stations(
          species = spp_list, 
          n = 550, 
          min_n_per_stratum = 4, 
          year = 2023, 
          vessel_names = c("vessel_1",  "vessel_2"), 
          ## set a location to save output
          output_dir = NULL)
```


#### *Station Allocation: Aleutian Islands 2022 Bottom Trawl Survey* 
This package holds the functions used to randomly draw and allocate survey 
stations for the Aleutian Islands and Gulf of Alaska bottom trawl survey. The 
code was originally produced by Paul von Szalay and Ned Laman. This package is
maintained by Zack Oyafuso and is intended to formalize the set of functions 
saved on the G drive (G:/GOA/R/survey planning/.RData).


## Inputs 

goa.planning() is the main wrapper function that performs the allocations and consists of nested subfunctions. You should be able to just run goa.planning() and the function will ask for several user inputs before executing its main program, including:

* Survey name (AI/GOA/EBSSLOPE)

* Oracle password for the AIGOA_WORK_DATA schema

* Total number of stations to allocate

* Survey year

* Output directory

The (pseudocode) hierarchy of the goa.planning() function is as such:

```
goa.planning() { ## main wrapper function

get.ai.stations() ## query available Aleutian Islands stations

allocate.effort() { ## allocate available stations across strata
get.planning.data() ## caluclate stratum-level catch stats
plot.allocations.by.stratum() { ## output allocation plots by stratum 
plot.strata.stations() ##  plot stations across stratum
}
}

pick.gridpoints() ## randomly sample stations based on stratified random allocations

vessel.allocation() ## allocate stations across vessels

}

```


## Output

* AIallocations.xlsx: vessel allocation, station id, stratum id, and lat/lon of sampled stations.

* AIallocations.csv: same as AIallocations.xlsx but in csv form.

* AllocationPlotsByStratum.pdf: output plots 

## Legal disclaimer
This repository is a software product and is not official communication of the National Oceanic and Atmospheric Administration (NOAA), or the United States Department of Commerce (DOC). All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the DOC or DOC bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation, or favoring by the DOC. The DOC seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by the DOC or the United States Government.
