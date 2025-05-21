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
  argv <- c("-e", "BRCA1", "-c", "BRCA", "-m", "cnaseg", "-a", "TRUE")
} else {
  # Get real command-line arguments
  argv <- commandArgs(trailingOnly = TRUE)
}

parser <- ArgumentParser()
parser$add_argument("-e", "--gene", type = "character", help = "gene")
parser$add_argument("-c", "--cancer", type = "character", help = "cancer type")
parser$add_argument("-m", "--mutation", type = "character", help = "mutation type")
parser$add_argument("-a", "--adj", type = "character", help = "Prop Adjustment (TRUE/FALSE or yes/no)")

args <- parser$parse_args(argv)
print(args)
prop_correction <- tolower(args$adj) %in% c("true", "t", "1", "yes", "y")
print(paste0("Proportion correction: ", prop_correction))

calculate_median_est_incidence_detail(
  date,
  args$cancer,
  args$gene,
  args$mutation,
  prop_correction
)
