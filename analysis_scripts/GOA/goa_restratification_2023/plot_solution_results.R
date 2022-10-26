###############################################################################
## Project:       Plot function
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Create function to output basic solution maps
###############################################################################


plot_solution_results <- function (file_name,
                                   grid_object,
                                   districts_object,
                                   sol_by_cell,
                                   strata_bounds){

  ## Import libraries
  grid_object$solution <- sol_by_cell

  ## Set up plot
  pdf(file = file_name, width = 6, height = 6)

  par(mfrow = c(3, 2))
  for (i in 1:5) {
    ireg <- c("Shumagin", "Chirikof", "Kodiak", "Yakutat", "Southeastern")[i]
    tmp <- terra::mask(x = grid_object, mask = nmfs[nmfs$area_name == ireg])
    tmp$solution <- as.integer(as.factor(tmp$solution))

    tmp_str_bound <- subset(strata_bounds, subset = Domain == i)

    plot(tmp, col = RColorBrewer::brewer.pal(n = max(tmp$solution),
                                             name = "Spectral")[tmp$solution],
         pch = 15,
         axes = F, main = ireg)
    legend(c("Shumagin" = "topleft", "Chirikof" = "topleft",
             "Kodiak" = "bottomright", "Yakutat" = "bottom",
             "Southeastern" = "bottomleft")[ireg],
           legend = paste0(round(tmp_str_bound$Lower_X1), "-",
                           round(tmp_str_bound$Upper_X1), " m, n = ",
                           tmp_str_bound$Allocation),
           fill = RColorBrewer::brewer.pal(n = max(tmp$solution),
                                           name = "Spectral"),
           cex = 0.7)
  }

  ## Close device
  dev.off()
}
