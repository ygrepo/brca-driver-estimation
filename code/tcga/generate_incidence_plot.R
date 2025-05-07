# load libraries
library(BoutrosLab.plotting.general)
library(here)
library(dplyr)

rm(list = ls())
date <- Sys.Date()
#setwd("/Users/yvesgreatti/github/brca-driver-estimation")

source(here("code", "tcga", "helper_functions.R"))
#source("code/tcga/helper_functions.R")

date <- "2025-05-07"
#date <- "2020-05-01"
#date <- "2025-05-06"
cancer <- "BRCA"
gene <- "BRCA2"
#mutation <- "lohdeletionseg"
mutation <- "deletion"
calculate_median_est_incidence_detail(date, cancer, gene, mutation)

