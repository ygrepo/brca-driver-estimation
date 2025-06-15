# load libraries
library(BoutrosLab.plotting.general)
library(argparse)
library(here)
library(dplyr)

rm(list = ls())

set.seed(42)
date <- Sys.Date()
date <- "2025-06-04" # for testing purposes

source(here("code", "tcga", "helper_functions.R"))
if (interactive()) {
  # Mimic command-line input
  argv <- c("-e", "BRCA1", 
            "-c", "BRCA", 
            "-m", "amplification", 
            "-a", "TRUE",
            "-l", "FALSE",
            "-t", "TRUE")
  } else {
  # Get real command-line arguments
  argv <- commandArgs(trailingOnly = TRUE)
}

parser <- ArgumentParser()
parser$add_argument("-e", "--gene", type = "character", help = "gene")
parser$add_argument("-c", "--cancer", type = "character", help = "cancer type")
parser$add_argument("-m", "--mutation", type = "character", help = "mutation type")
parser$add_argument("-a", "--adj", type = "character", help = "Prop Correction (TRUE/FALSE or yes/no)")
parser$add_argument("-l", "--loh", type = "character", help = "LOH Correction (TRUE/FALSE or yes/no)")
parser$add_argument("-t", "--tp53", type = "character", help = "TP53 Correction (TRUE/FALSE or yes/no)")

args <- parser$parse_args(argv)
print(args)
prop_correction <- tolower(args$adj) %in% c("true", "t", "1", "yes", "y")
print(paste0("Proportion correction: ", prop_correction))
if (args$loh == "TRUE" || tolower(args$loh) %in% c("true", "t", "1", "yes", "y")) {
  # filter for LOH variants
  args$loh <- "TRUE"
  print("Filtering for LOH variants")
} else {
  args$loh <- "FALSE"
}

tp53_correction <- tolower(args$tp53) %in% c("true", "t", "1", "yes", "y")
if (tp53_correction) {
  tp53_correction <- "TRUE"
} else {
  tp53_correction <- "FALSE"
}
message(paste0("TP53 correction: ", tp53_correction))
results <- calculate_median_est_incidence_detail(
  date,
  args$cancer,
  args$gene,
  args$mutation,
  prop_correction,
  loh_correction = args$loh,
  tp53_correction = tp53_correction
)
results

filename <- paste(date, args$cancer, args$gene, args$mutation, 
                  "Prop", args$adj,
                  "loh", args$loh,
                  "tp53", tp53_correction,
                  "ks_test.tsv",
                  sep = "_"
)
filename <- here("output", "data", "TCGA/European", filename)

# write to file
write.table(
  results,
  file = filename,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
print(paste0(
  "Saving results to ", filename
))

