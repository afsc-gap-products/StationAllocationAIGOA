#' Wrapper function for plotting allocation plots by stratum
#'
#' @description
#' Creates AllocationPlotsByStratum.pdf, which includes output from
#' AIGOASurveyPlanning::plot.strata.stations as well as other output (finish!)
#'
#' Created 12/18/13
#'
#' Notes from Ned: added ", believeNRows = FALSE" to deal with version 3.0.1 of
#' R in the 64-bit environment
#'
#' @author Ned Laman \email{ned.laman@@noaa.gov}
#'
#' @param answer.years numeric array of dimensions (number of strata, number of
#' species, number of years) of the weighted optimal allocation (see reference...)
#' @param number.tows numeric. Vector output from
#' AIGOASurveyPlanning::allocate.effort
#' @param strata Dataframe. Output from AIGOASurveyPlanning::allocate.effort
#' @param channel open connection object. Created from RODBC::odbcConnect()
#' @param survey character string. Defaults to "AI" for Aleutian Islands
#' @param output.file character string. Output directory.
#'
#' @return AllocationPlotsByStratum.pdf saved in the output.file directory

plot.allocations.by.stratum <- function(answer.years,
                                        number.tows,
                                        strata,
                                        channel,
                                        survey,
                                        output.file){


  cat("\nNow retrieving species table from oracle table RACEBASE.species...\n")

  species <- sqlQuery(channel = channel, query = "select * from racebase.species", believeNRows = FALSE)
  names(species) <- casefold(names(species))

  cat("\nPlotting allocations by stratum and saving to PDF...\n")

  pdf(paste0(output.file, "\\AllocationPlotsByStratum.pdf"))

  plot.strata.stations(number.tows, strata, survey)

  species.summary <- apply(answer.years, 2., sum)
  all.sum <- sum(species.summary)

  for(sp in names(species.summary)){
    print(sp)
    common.name <- species$common_name[species$species_code == sp]
    tows <- t(answer.years[, sp,  ])
    percent.survey <- round((species.summary[sp]/all.sum) * 100, 1)

    if(max(tows) >= 3. | percent.survey > 1.){
      x <- barplot(tows, legend.text = names(
        answer.years[1, 1,]), axes = F, main = common.name, xlab = "Stratum", ylab = "% of survey", las = 2)
      #col = 1.:length(names(answer.years[1., 1.,  ])) + 1,
      box()
      text(0.85 * max(x), (1. - (length(names(answer.years[1., 1.,  ])) * 0.009)) * max(apply(tows, 2., sum)), paste(
        percent.survey, "% of survey"))
      # locator(n=1)
    }
  }

  stratum.summary <- apply(answer.years, 1, sum)

  for(st in names(stratum.summary)) {
    print(st)
    tows <- t(answer.years[st,  ,  ])
    percent.survey <- round((stratum.summary[st]/all.sum) * 100, 1)
    x <- barplot(tows, col = 1.:length(names(answer.years[1., 1.,  ])) + 1., legend.text = names(answer.years[1,1,]),
                 main = paste("Stratum", st), las = 2, xlab = "Species", ylab = "Number of Tows")
    box()
    text(0.85 * max(x), (1. - (length(names(answer.years[1., 1.,  ])) * 0.009)) * max(apply(tows, 2., sum)), paste(
      percent.survey, "% of survey"))
    # locator(n=1)
  }
  dev.off()
}
