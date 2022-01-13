# AIGOASurveyPlanning

This package holds the functions used to allocate stations for the Aleutian Islands bottom trawl survey. The code was originally produced by Paul von Szalay and Ned Laman. This package is maintained by Zack Oyafuso and is intended to formalize the set of functions saved on the G drive (G:/GOA/R/survey planning/.RData).  

## Packages Required

```
library(devtools) ## initial installation of package
library(RODBC) ## connection to Oracle
library(XLConnect) ## reading and writing Excel workbooks and sheets
```

## Inputs 

goa.planning() is the main wrapper function that performs the allocations and consists of nested subfunctions. The function will ask for several user inputs before executing its main program:

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

