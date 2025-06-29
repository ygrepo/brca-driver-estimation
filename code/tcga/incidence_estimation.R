### DESCRIPTION ###################################################################################
# Kuan's driver estimates

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


source(here("code", "tcga", "helper_functions.R"))
# set seed for reproducibility
set.seed(42)
if (interactive()) {
  # Mimic command-line input
  argv <- c("-e", "BRCA1", 
            "-c", "BRCA", 
            "-m", "cnaseg", 
            "-a", "TRUE",
            "-l", "FALSE",
            "-t", "TRUE",
            "-n", "1",
            "--hrd", "yes")
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
parser$add_argument("-n", "--n_runs", type = "integer", default = 10000, help = "Number of bootstrap runs (default: 10000)")
parser$add_argument("--hrd", type = "character", default = 'yes', help = "Exclude HRD-gene variant")
args <- parser$parse_args(argv)
print(args)

### DATA PROCESSING ###############################################################################

var_anno <- read.delim(
  here("data", "TCGA", "PCA_pathVar_integrated_filtered_adjusted_ancestry.tsv"),
  as.is = TRUE
)
#var_anno |> glimpse() 

if (args$loh == "TRUE" || tolower(args$loh) %in% c("true", "t", "1", "yes", "y")) {
  # filter for LOH variants
  args$loh <- "TRUE"
  print("Filtering for LOH variants")
  var_anno <- var_anno[var_anno$LOH_classification == "Significant", ]
}

# find samples with variants in specified cancer
ddr <- var_anno[var_anno$HUGO_Symbol == args$gene & 
                  var_anno$cancer == args$cancer, "bcr_patient_barcode"]


# read in patient annotation
pat_anno <- read.delim(
  here("data", "TCGA", "clinical_PANCAN_patient_with_followup.tsv"),
  as.is = TRUE
)

# subset down to cancer type of interest
can_anno <- pat_anno[pat_anno$acronym == args$cancer, ]

# Select European only ----
target_ancestry <- readxl::read_xlsx(here("data", "TCGA", "mmc2.xlsx"),
  "S1 Calls per Patient",
  skip = 1
)
table(target_ancestry$consensus_ancestry)

target_ancestry <- target_ancestry |>
  select(patient, consensus_ancestry) |>
  filter(consensus_ancestry == "eur")

#length(intersect(target_ancestry$patient, can_anno$bcr_patient_barcode))
can_anno <- can_anno |>
  filter(bcr_patient_barcode %in% target_ancestry$patient)

# get mutation rate
mut_rate <- get_mutation_rate(type = args$mutation, anno = can_anno)


# find samples without ddr
wt <- can_anno$bcr_patient_barcode[!can_anno$bcr_patient_barcode %in% var_anno$bcr_patient_barcode]


# only keep samples with mutation rate
ddr <- ddr[ddr %in% mut_rate$bcr_patient_barcode]
wt <- wt[wt %in% mut_rate$bcr_patient_barcode]

if (args$hrd == "TRUE" || tolower(args$hrd) %in% c("true", "t", "1", "yes", "y")) {
  # filter for HRD-genes variants
  # filter for specified cancer and gene
  hrd_genes <- c("ATM", "CHEK2", "PALB2", "BARD1", "RAD51C", "RAD51D", "ATR")
  hrdbarccode <- var_anno |>
    dplyr::filter(HUGO_Symbol %in% hrd_genes) |>
    dplyr::filter(Overall_Classification %in% c("Pathogenic", "Likely Pathogenic")) |>
    dplyr::select(bcr_patient_barcode) |>
    dplyr::pull(bcr_patient_barcode) |>
    unique()
  
  # remove samples with HRD-genes variants
  print("Filtering for HRD-gene variants")
  cat("Number of HRD-gene variants: ", length(hrdbarccode), "\n")    
  l = length(intersect(wt, hrdbarccode))
  cat("Number of WT samples with HRD-gene variants: ", l, "\n")
  cat("Number of WT samples: ", length(wt), "\n")
  wt <- wt[!wt %in% hrdbarccode]
  cat("Number of WT samples after filtering: ", length(wt), "\n")
  
}  

# print the number of ddr samples
print(paste0("Number of ", args$cancer, " samples with variant in ", args$gene, ": ", length(ddr)))
n_runs <- 10000
#n_runs <- 1

print(paste0("Number of runs:", n_runs))

prop_correction <- tolower(args$adj) %in% c("true", "t", "1", "yes", "y")
print(paste0("Proportion correction: ", args$adj))

tp53_correction <- tolower(args$tp53) %in% c("true", "t", "1", "yes", "y")
message(paste0("TP53 correction: ", tp53_correction))

tp53_status_df <- NULL
if (tp53_correction) {
  message("Reading MC3 MAF for TP53 variants")
  mc3 <- read.delim(here("data", "TCGA", "mc3.v0.2.8.PUBLIC.maf.gene_vclass_HGVSp_sample_likelyDriverLoose_aggregated_matrix.tsv"), as.is = TRUE)
  df <- mc3 |>
    filter(Hugo_Symbol == "TP53") 
  all_samples <-  unique(colnames(df)[2:ncol(df)])
  all_samples <- unique(substr(all_samples, 1, 12))
  all_samples <- gsub("\\.", "-", all_samples)
  df <- df[, colSums(is.na(df)) < nrow(df)]
  tp53_del_samples<- unique(colnames(df)[2:ncol(df)])
  tp53_del_samples <- unique(substr(tp53_del_samples, 1, 12))
  tp53_del_samples <- gsub("\\.", "-", tp53_del_samples)
  tp53_status_df <- tibble::tibble(
    Sample.ID = all_samples,
    TP53_Status = ifelse(all_samples %in% tp53_del_samples, "TP53_DEL", "TP53_WT")
  )
  l = length(intersect(all_samples, tp53_del_samples))
  cat("Number of samples with TP53 deletion: ", l, "\n")
  l = length(intersect(ddr, all_samples))
  cat("Number of samples with DDR variant: ", l, "\n")
  l =  length(intersect(wt, all_samples))
  cat("Number of samples without DDR variant: ", l, "\n")
  l = length(intersect(ddr, tp53_del_samples))
  cat("Number of samples with DDR variant and TP53 deletion: ", l, "\n")
  l = length(intersect(wt, tp53_del_samples))
  cat("Number of samples without DDR variant and TP53 deletion: ", l, "\n")
}




# run bootstrap
results <- do.call(rbind, sapply(
  1:n_runs,
  calculate_mutation_rate_ratio,
  mut_rate = mut_rate,
  ddr = ddr,
  wt = wt,
  anno = can_anno,
  cancer = args$cancer,
  gene = args$gene,
  prop_correction = prop_correction,
  loh_correction = args$loh,
  tp53_status_df = tp53_status_df,
  simplify = FALSE
))

# generate file name
if (tp53_correction) {
  tp53_correction <- "TRUE"
} else {
  tp53_correction <- "FALSE"
}
cat("TP53 correction: ", tp53_correction, "\n")
filename <- paste(date, args$cancer, args$gene, args$mutation, 
                  "Prop", args$adj,
                  "loh", args$loh,
                  "tp53", tp53_correction,
  "incidence_estimates.tsv",
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
