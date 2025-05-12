# load libraries
library(BoutrosLab.plotting.general)
library(here)
library(dplyr)

rm(list = ls())
date <- Sys.Date()
#setwd("/Users/yvesgreatti/github/brca-driver-estimation")

source(here("code", "tcga", "helper_functions.R"))
#source("code/tcga/helper_functions.R")

date <- "2020-05-01"
#date <- "2020-04-22"
cancer <- "BRCA"
gene <- "BRCA1"
#mutation <- "lohdeletionseg"
mutation <- "cnaseg"
adj_flag <- TRUE
calculate_median_est_incidence_detail(date, 
                                      cancer, gene, mutation, 
                                      adj_flag = adj_flag)

