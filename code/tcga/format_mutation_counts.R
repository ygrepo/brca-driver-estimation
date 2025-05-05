### DESCRIPTION ###################################################################################
# reformat cna calls; calculate number of deletions

### USAGE #########################################################################################
#

### PREAMBLE ######################################################################################
# load libraries
library(argparse)
library(dplyr)
rm(list = ls())
date <- Sys.Date()
setwd("/Users/yvesgreatti/github/brca-driver-estimation")
### DATA PROCESSING ###############################################################################
# read in cna data
segments <- read.delim("data/TCGA/TCGA_mastercalls.abs_segtabs.fixed.txt.gz")

# remove copy number neutral probes
cna <- segments[which(segments$Modal_Total_CN != 2), ]
# calculate the number of segments
cna_num <- as.data.frame(table(cna$Sample))
colnames(cna_num) <- c("bcr_patient_barcode", "count")
# add any samples without cnas
no_cna <- data.frame(
  bcr_patient_barcode = unique(segments$Sample[!segments$Sample %in% cna_num$bcr_patient_barcode]),
  count = 0
)
cna_num <- rbind(cna_num, no_cna)
# fin median for duplicated samples
cna_num$bcr_patient_barcode <- substr(cna_num$bcr_patient_barcode, 1, 12)
cna_num <- aggregate(cna_num$count, list(cna_num$bcr_patient_barcode), median, na.rm = TRUE)
colnames(cna_num) <- c("bcr_patient_barcode", "count")

# write to file
write.table(
  cna_num,
  file = "data/TCGA/mutations/TCGA_segments.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# calculate number of bps
cna_bp <- aggregate(cna$Length, list(cna$Sample), sum, na.rm = TRUE)
colnames(cna_bp) <- c("bcr_patient_barcode", "count")
# add samples without any cnas
cna_bp <- rbind(cna_bp, no_cna)
# fin median for duplicated samples
cna_bp$bcr_patient_barcode <- substr(cna_bp$bcr_patient_barcode, 1, 12)
cna_bp <- aggregate(cna_bp$count, list(cna_bp$bcr_patient_barcode), median, na.rm = TRUE)
colnames(cna_bp) <- c("bcr_patient_barcode", "count")

# write to file
write.table(
  cna_bp,
  file = "data/TCGA/mutations/TCGA_total_cna_bp.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

### CALCULATE NUMBER OF DELETIONS #################################################################

# count the number of deletions per patient


del <- aggregate(
  segments$Modal_Total_CN,
  list(segments$Sample),
  function(x) {
    sum(x < 2, na.rm = TRUE)
  }
)
colnames(del) <- c("bcr_patient_barcode", "count")
del$bcr_patient_barcode <- substr(del$bcr_patient_barcode, 1, 12)

# remove duplicated lines and find median value for duplicated samples
del <- del[!duplicated(del), ]
del_dedup <- aggregate(del$count, list(del$bcr_patient_barcode), median, na.rm = TRUE)
colnames(del_dedup) <- c("bcr_patient_barcode", "count")

# write to file
write.table(
  del_dedup,
  file = "data/TCGA/mutations/TCGA_segments_deletions_only.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# count the number of deletions per patient including LOH
segments <- read.delim("data/TCGA/TCGA_mastercalls.abs_segtabs.fixed.txt.gz")

del_loh <- aggregate(
  (segments$Modal_Total_CN < 2 & segments$LOH == 1),
  list(segments$Sample),
  sum
)
colnames(del_loh) <- c("bcr_patient_barcode", "count")
del_loh$bcr_patient_barcode <- substr(del_loh$bcr_patient_barcode, 1, 12)
del_loh <- del_loh[!duplicated(del_loh), ]
del_loh_dedup <- aggregate(del_loh$count,
  list(del_loh$bcr_patient_barcode), median,
  na.rm = TRUE
)
colnames(del_loh_dedup) <- c("bcr_patient_barcode", "count")

# write to file
write.table(
  del_loh_dedup,
  file = "data/TCGA/mutations/TCGA_LOH_deletions.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


### CALCULATE PGA FOR EACH PATIENT ################################################################
# read in cna data
segments <- read.delim("data/TCGA/TCGA_mastercalls.abs_segtabs.fixed.txt.gz")
segments <- segments[!is.na(segments$Length), ]

# only keep deletions
segments_del <- segments[segments$Modal_Total_CN < 2, ]

# count the number of abberated bases per patient
ab_bp <- aggregate(
  segments_del$Length,
  list(segments_del$Sample),
  sum
)
colnames(ab_bp) <- c("bcr_patient_barcode", "count")

# add samples without deletion
missing <- segments$Sample[!segments$Sample %in% ab_bp$bcr_patient_barcode]
nodel <- data.frame(
  bcr_patient_barcode = unique(missing),
  count = 0
)
ab_bp <- rbind(ab_bp, nodel)
# reformat id
ab_bp$bcr_patient_barcode <- substr(ab_bp$bcr_patient_barcode, 1, 12)

# remove duplicated lines and find median value for duplicated samples
ab_bp <- ab_bp[!duplicated(ab_bp), ]
ab_bp_dedup <- aggregate(ab_bp$count, list(ab_bp$bcr_patient_barcode), median, na.rm = TRUE)
colnames(ab_bp_dedup) <- c("bcr_patient_barcode", "count")

# write to file
write.table(
  ab_bp_dedup,
  file = "data/TCGA/mutations/TCGA_cna_deletion_bp.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

### CALCULATE PGA FOR EACH PATIENT ################################################################
# read in cna data
segments <- read.delim("data/TCGA/TCGA_mastercalls.abs_segtabs.fixed.txt.gz")
segments <- segments[!is.na(segments$Length), ]

# only keep deletions with LOH
segments_del <- segments[segments$Modal_Total_CN < 2 & segments$LOH == 1, ]

# count the number of abberated bases per patient
ab_bp <- aggregate(
  segments_del$Length,
  list(segments_del$Sample),
  sum
)
colnames(ab_bp) <- c("bcr_patient_barcode", "count")

# add samples without deletion
missing <- segments$Sample[!segments$Sample %in% ab_bp$bcr_patient_barcode]
nodel <- data.frame(
  bcr_patient_barcode = unique(missing),
  count = 0
)
ab_bp <- rbind(ab_bp, nodel)
# reformat id
ab_bp$bcr_patient_barcode <- substr(ab_bp$bcr_patient_barcode, 1, 12)

# remove duplicated lines and find median value for duplicated samples
ab_bp <- ab_bp[!duplicated(ab_bp), ]
ab_bp_dedup <- aggregate(ab_bp$count, list(ab_bp$bcr_patient_barcode), median, na.rm = TRUE)
colnames(ab_bp_dedup) <- c("bcr_patient_barcode", "count")

# write to file
write.table(
  ab_bp_dedup,
  file = "data/TCGA/mutations/TCGA_LOH_deletions_bp.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

### FORMAT INDELS #################################################################################
# read in mc3
rm(list = ls())

# mc3 <- read.delim('data/MC3_Ellrott_CellSys2018/mc3.v0.2.8.CONTROLLED.maf.gz', as.is = TRUE)
mc3 <- read.delim("data/TCGA/mc3.v0.2.8.PUBLIC.maf.gz", as.is = TRUE)
# find all samples in mc3
all_samples <- substr(mc3$Tumor_Sample_Barcode, 1, 12)
# subset down to deletions
mc3_del <- mc3[mc3$Variant_Type == "DEL", ]
# calculate length of deletion
mc3_del$length <- nchar(mc3_del$Tumor_Seq_Allele1) - nchar(mc3_del$Tumor_Seq_Allele2)
mc3_del[mc3_del$Tumor_Seq_Allele2 == "-", "length"] <- mc3_del[mc3_del$Tumor_Seq_Allele2 == "-", "length"] + 1
# reformat
indel_df <- mc3_del[, c("Tumor_Sample_Barcode", "length")]
# count up indels
indel_df <- aggregate(indel_df$length, list(indel_df$Tumor_Sample_Barcode), sum)
colnames(indel_df) <- c("bcr_patient_barcode", "count")
indel_df$bcr_patient_barcode <- substr(indel_df$bcr_patient_barcode, 1, 12)

# include other samples in mc3 without indels
no_indel <- data.frame(
  bcr_patient_barcode = unique(all_samples[!all_samples %in% indel_df$bcr_patient_barcode]),
  count = 0
)
indel_df <- rbind(indel_df, no_indel)


# write to table
write.table(
  indel_df,
  file = "data/TCGA/mutations/TCGA_indels_deletions_bp_only.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# subset down to insertions
mc3_ins <- mc3[mc3$Variant_Type == "INS", ]
# calculate length of insertion
mc3_ins$length <- nchar(mc3_ins$Tumor_Seq_Allele2) - nchar(mc3_ins$Tumor_Seq_Allele1)
mc3_ins[mc3_ins$Tumor_Seq_Allele1 == "-", "length"] <- mc3_ins[mc3_ins$Tumor_Seq_Allele1 == "-", "length"] + 1
# reformat
ins_df <- mc3_ins[, c("Tumor_Sample_Barcode", "length")]
# count up indels
ins_df <- aggregate(ins_df$length, list(ins_df$Tumor_Sample_Barcode), sum)
colnames(ins_df) <- c("bcr_patient_barcode", "count")
ins_df$bcr_patient_barcode <- substr(ins_df$bcr_patient_barcode, 1, 12)

# include other samples in mc3 without inss
no_ins <- data.frame(
  bcr_patient_barcode = unique(all_samples[!all_samples %in% ins_df$bcr_patient_barcode]),
  count = 0
)
ins_df <- rbind(ins_df, no_ins)


# write to table
write.table(
  ins_df,
  file = "data/TCGA/mutations/TCGA_indels_insertions_bp_only.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
