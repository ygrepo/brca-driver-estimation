### DESCRIPTION ###################################################################################
# Kuan's driver estimates

### USAGE #########################################################################################
# ensure in same directory as mutation and annotation files

### PREAMBLE ######################################################################################
# load libraries
library(argparse)
library(dplyr)
library(tidyr)
library(here)
library(readxl)
library(openxlsx)

rm(list = ls())
date <- Sys.Date()
#setwd("/Users/yvesgreatti/github/brca-driver-estimation")

source(here("code", "tcga", "helper_functions.R"))

# setwd('~/GermlineSomaticAssociations/genome-wide/output/somatic_gwas/pancan/driver_rate_comparison')
### OBTAIN COMMAND LINE ARGUMENTS #################################################################
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
### DATA PROCESSING ###############################################################################
# read in variant annotation
var_anno <- read.delim(
  #' data/TCGA_PanCanAtlas_2018/Germline_Huang_Cell2018/TCGA_ancestry_PC.txt',
  here("data", "TCGA", "PCA_pathVar_integrated_filtered_adjusted_ancestry.tsv"),
  as.is = TRUE
)
# find samples with variants in specified cancer
loh_flag <- FALSE
var_anno <- var_anno[var_anno$HUGO_Symbol == args$gene & var_anno$cancer == args$cancer, ]
if (loh_flag) {
  var_anno <- var_anno[var_anno$LOH_classification == "Significant", ]
}
ddr <- unique(var_anno$bcr_patient_barcode)
print(length(ddr))


# ddr <- var_anno[var_anno$HUGO_Symbol == args$gene & var_anno$cancer == args$cancer, 
#                 "bcr_patient_barcode"]

# ddr <- var_anno[var_anno$HUGO_Symbol == args$gene & var_anno$cancer == args$cancer,'bcr_patient_barcode']

# read in patient annotation
pat_anno <- read.delim(
  here("data", "TCGA", "clinical_PANCAN_patient_with_followup.tsv"),
  as.is = TRUE
)
# fn <- here("data", "TCGA", "TCGA-CDR-SupplementalTableS1.xlsx")
# pat_anno <- read_excel(fn, sheet = "TCGA-CDR")
#can_anno <- pat_anno[pat_anno$type == args$cancer, ]
can_anno <- pat_anno[pat_anno$acronym == args$cancer,]
# only consider individuals of European descent
can_anno <- can_anno[which(can_anno$race == "WHITE"), ]

# 
# # subset down to cancer type of interest
# can_anno <- pat_anno[pat_anno$acronym == args$cancer, ]
# # only consider individuals of European descent
# can_anno <- can_anno[which(can_anno$race == "WHITE"), ]

# get mutation rate
mut_rate <- get_mutation_rate(type = args$mutation, anno = can_anno)

# find samples without ddr
wt <- can_anno$bcr_patient_barcode[!can_anno$bcr_patient_barcode %in% var_anno$bcr_patient_barcode]

# only keep samples with mutation rate
ddr <- ddr[ddr %in% mut_rate$bcr_patient_barcode]
wt <- wt[wt %in% mut_rate$bcr_patient_barcode]

# print the number of ddr samples
print(paste0(
  "Number of ", args$cancer,
  " samples with variant in ", args$gene, ": ", length(ddr)
))

# run bootstrap
adj_flag <- FALSE
if (adj_flag) {
  adj_dir = "adjusted"
} else {
  adj_dir = "unadjusted"
}
results <- do.call(rbind, sapply(
  1:10,
  #1:10000,
  calculate_mutation_rate_ratio,
  date = date,
  mut_rate = mut_rate,
  ddr = ddr,
  wt = wt,
  anno = can_anno,
  cancer = args$cancer,
  gene = args$gene,
  mutation = args$mutation,
  adj_flag = adj_flag,
  simplify = FALSE
))

# generate file name
filename <- paste(date, args$cancer, args$gene, args$mutation, "incidence_estimates.tsv", sep = "_")
filename <- here("output", "data", "TCGA", adj_dir, filename)

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

