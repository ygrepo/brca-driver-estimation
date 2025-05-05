# load libraries
library(BoutrosLab.plotting.general)
library(dplyr)

rm(list = ls())
date <- Sys.Date()
setwd("/Users/yvesgreatti/github/brca-driver-estimation")

source("code/tcga/helper_functions.R")

date <- "2025-05-05"
cancer <- "BRCA"
gene <- "BRCA1"
mutation <- "lohdeletionseg"
calculate_median_est_incidence_2(date, cancer, gene, mutation)


calculate_mutation_rate_ratio <- function(int, mut_rate, ddr, wt, anno, cancer,
                                          include_loh = TRUE,
                                          loh_column = "loh_rate",
                                          max_drivers = 15) {
  # 1. Sample ddr and matched wt
  ddr_sample <- sample(ddr, length(ddr), replace = TRUE)
  wt_sample  <- switch(cancer,
                       BRCA = standardize_clinical_characteristics_breast(anno, wt, ddr),
                       OV   = standardize_clinical_characteristics_ovarian(anno, wt, ddr),
                       stop("Please specify cancer = 'BRCA' or 'OV'")
  )
  
  # 2. If LOH/second-hits are in their own column, add them into the rate
  if (include_loh && loh_column %in% colnames(mut_rate)) {
    combined_rate <- mut_rate$rate + mut_rate[[loh_column]]
  } else {
    combined_rate <- mut_rate$rate
  }
  
  # 3. Compute medians and somatic‐to‐wt ratio
  ddr_median <- median(combined_rate[ddr_sample], na.rm=TRUE)
  wt_median  <- median(combined_rate[wt_sample], na.rm=TRUE)
  ratio      <- ddr_median / wt_median
  
  # 4. Build incidence_k = ratio^(k-1) for k = 1..max_drivers,
  #    since germline hit counts as the “first” event.
  ks        <- seq_len(max_drivers)
  incidences <- ratio^(ks - 1)
  names(incidences) <- paste0("incidence_", ks)
  
  # 5. Return
  data.frame(
    int        = int,
    num        = length(ddr),
    wt_median  = wt_median,
    ddr_median = ddr_median,
    ratio      = ratio,
    as.list(incidences),
    check.names = FALSE
  )
}
