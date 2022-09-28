##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project:       Calculate MS allocation of STRS
## Author:        Zack Oyafuso (zack.oyafuso@noaa.gov)
## Description:   Gulf of Alaska 2023 bottom trawl survey
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())






##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   CVs under SRS: upper CV bounds for MS allocation
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   single species CV: lower CV bounds for MS allocation
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
frame$X1 <- as.integer(factor(grid_goa_sp@data$STRATUM))
strs_stats <- SamplingStrata::buildStrataDF(frame)
strs_stats <- strs_stats[order(as.integer(strs_stats$STRATO)), ]
strs_stats$N <- strs_stats$N / n_years

ss_cv <- data.frame()
ss_allocations <- matrix(nrow = nrow(strs_stats), ncol = ns_opt,
                         dimnames = list(levels(factor(grid_goa_sp@data$STRATUM)),
                                         common_names_opt))

for (ispp in 1:ns_opt) {
  error_df <- data.frame("DOM" = "DOM1",
                         "CV1" = as.numeric(srs_cv[ispp]),
                         "domainvalue"  = 1)
  temp_stratif <- cbind(DOM1 = 1,
                        CENS = 0,
                        strs_stats[, c(names(strs_stats)[1:2], paste0(c("M", "S"), ispp)) ])
  names(temp_stratif)[-c(1:4)] <- paste0(c("M", "S"), 1)

  temp_bethel <- SamplingStrata::bethel(
    errors = error_df,
    stratif = temp_stratif,
    realAllocation = T,
    printa = T,
    minnumstrat = 4)
  temp_n <- as.integer(sum(temp_bethel))
  temp_cv <- as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "])

  while (temp_n != 550){
    over_under <- temp_n > 550
    CV_adj <- ifelse(over_under == TRUE,
                     yes = 1.001,
                     no = 0.999)

    error_df$CV1 <- temp_cv* CV_adj
    temp_bethel <- SamplingStrata::bethel(stratif = temp_stratif,
                                          errors = error_df,
                                          printa = TRUE,
                                          minnumstrat = 4)

    temp_n <- sum(as.numeric(temp_bethel))
    temp_cv <- as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "])
    print(paste0("n = ", temp_n, ", CV = ", temp_cv) )
  }

  ss_cv <- rbind(ss_cv, data.frame(species = common_names_opt[ispp],
                                   ss_cv = temp_cv))
  ss_allocations[, ispp] <- temp_bethel
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##   multi species CV:
##   Start at SRS CV --> bethel algorithm --> optimal sampling size (n_srs)
##   If n_srs < 550 stations, reduce CV constraints across species by 0.1%
##   using the CVs optimized for each species (ss_cv) as a lower bound.
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
error_df <- cbind("DOM" = "DOM1",
                  srs_cv,
                  "domainvalue"  = 1)
names(error_df)[2:(1 + ns_opt)] <- paste0("CV", 1:ns_opt)

temp_stratif <- cbind(DOM1 = 1,
                      CENS = 0,
                      strs_stats[, c(names(strs_stats)[1:2],
                                     paste0("M", 1:ns_opt),
                                     paste0("S", 1:ns_opt)) ])

temp_bethel <- SamplingStrata::bethel(
  errors = error_df,
  stratif = temp_stratif,
  realAllocation = T,
  printa = T,
  minnumstrat = 4)
temp_n <- sum(ceiling(temp_bethel))


while (temp_n != 550){
  over_under <- temp_n > 550
  CV_adj <- ifelse(over_under == TRUE,
                   yes = 1.001,
                   no = 0.999)

  updated_cv_constraint <-
    as.numeric(attributes(temp_bethel)$outcv[, "PLANNED CV "]) * (CV_adj) +
    ss_cv$ss_cv * (1  - CV_adj)

  error_df[, paste0("CV", 1:ns_opt)] <- as.numeric(updated_cv_constraint)

  temp_bethel <- SamplingStrata::bethel(stratif = temp_stratif,
                                        errors = error_df,
                                        printa = TRUE,
                                        minnumstrat = 4)

  temp_n <- sum(as.numeric(temp_bethel))

  print(paste0("n = ", temp_n) )
}

ms_cv <- as.numeric(updated_cv_constraint)
ms_allocation <- as.integer(ceiling(temp_bethel))
names(ms_allocation) <- levels(factor(grid_goa_sp@data$STRATUM))
