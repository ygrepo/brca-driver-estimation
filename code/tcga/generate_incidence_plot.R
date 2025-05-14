# load libraries
library(BoutrosLab.plotting.general)
library(argparse)
library(here)
library(dplyr)

rm(list = ls())
date <- Sys.Date()
#setwd("/Users/yvesgreatti/github/brca-driver-estimation")

source(here("code", "tcga", "helper_functions.R"))
if (interactive()) {
  # Mimic command-line input
  #argv <- c("-e", "BRCA1", "-c", "OV", "-m", "deletion")
  argv <- c("-e", "BRCA1", "-c", "BRCA", "-m", "cna")
} else {
  # Get real command-line arguments
  argv <- commandArgs(trailingOnly = TRUE)
}

parser <- ArgumentParser()
parser$add_argument("-e", "--gene", type = "character", help = "gene")
parser$add_argument("-c", "--cancer", type = "character", help = "cancer type")
parser$add_argument("-m", "--mutation", type = "character", help = "mutation type")
args <- parser$parse_args(argv)
print(args)


#source("code/tcga/helper_functions.R")

#date <- "2020-05-01"
#date <- "2020-04-22"
#cancer <- "BRCA"
#gene <- "BRCA1"
#mutation <- "lohdeletionseg"
#mutation <- "cnaseg"
#adj_flag <- TRUE
calculate_median_est_incidence_detail(date, 
                                      args$cancer, args$gene, args$mutation)
# , 
#                                       adj_flag = adj_flag)

