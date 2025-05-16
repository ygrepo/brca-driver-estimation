# load libraries
library(BoutrosLab.plotting.general)
library(argparse)
library(here)
library(dplyr)

rm(list = ls())

set.seed(42)
date <- Sys.Date()

source(here("code", "tcga", "helper_functions.R"))
if (interactive()) {
  # Mimic command-line input
  argv <- c("-e", "BRCA1", "-c", "BRCA", "-m", "cnaseg")
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

calculate_median_est_incidence_detail(
  date,
  args$cancer,
  args$gene,
  args$mutation
)
