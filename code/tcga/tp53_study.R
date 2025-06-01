library(argparse)
library(dplyr)
library(here)


rm(list = ls())
date <- Sys.Date()

# mc3 <- read.delim('data/MC3_Ellrott_CellSys2018/mc3.v0.2.8.CONTROLLED.maf.gz', as.is = TRUE)
mc3 <- read.delim(here("data", "TCGA", "mc3.v0.2.8.PUBLIC.maf.gz"), as.is = TRUE)

df <- mc3 |>
  filter(Hugo_Symbol == "TP53") 
# 
#   select(Tumor_Sample_Barcode, Variant_Type) 
# 
# # find all samples in mc3
# all_samples <- substr(mc3$Tumor_Sample_Barcode, 1, 12)
# subset down to deletions
df_del <- df[df$Variant_Type == "DEL", ]
# 
# # reformat
# indel_df <- mc3_del[, c("Tumor_Sample_Barcode", "length")]
# # count up indels
# indel_df <- aggregate(indel_df$length, list(indel_df$Tumor_Sample_Barcode), sum)
# colnames(indel_df) <- c("bcr_patient_barcode", "count")
# indel_df$bcr_patient_barcode <- substr(indel_df$bcr_patient_barcode, 1, 12)