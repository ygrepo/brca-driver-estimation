### DESCRIPTION ###################################################################################
# create densityplots for confounding factors

### USAGE #########################################################################################
#

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(dplyr)
rm(list = ls())

setwd("/Users/yvesgreatti/github/brca-driver-estimation")

date <- Sys.Date()
### DATA PROCESSING ###############################################################################
# read in aneuploidy file
aneu <- read.csv("data/TCGA/mutations/TCGA_aneuploidy.csv", as.is = TRUE, skip = 1)
# only keep ovarian and breast
aneu <- aneu[aneu$Type %in% c("BRCA", "OV"), ]
# subset id
aneu$Sample <- substr(aneu$Sample, 1, 12)

# read in variant annotation
var_anno <- read.delim(
  #' data/TCGA_PanCanAtlas_2018/Germline_Huang_Cell2018/TCGA_ancestry_PC.txt',
  "data/TCGA/PCA_pathVar_integrated_filtered_adjusted_ancestry.tsv",
  as.is = TRUE
)
# generate tissue variant pairs
pairs <- expand.grid(cancer = c("BRCA", "OV"), gene = c("BRCA1", "BRCA2"))
# find samples with variants in specified cancer
ddr <- do.call(rbind, apply(
  pairs,
  1,
  function(x) {
    tmp <- unlist(var_anno[var_anno$HUGO_Symbol == x["gene"] & var_anno$cancer == x["cancer"], "bcr_patient_barcode"])
    data.frame(Type = x["cancer"], Variant = x["gene"], Sample = tmp)
  }
))
# create densityplot of genomic doublings
plot_data <- aggregate(aneu$Genome_doublings, list(aneu$Sample), median)
colnames(plot_data) <- c("Sample", "GD")
plot_data <- merge(plot_data, unique(aneu[, c("Sample", "Type")]), by = "Sample")
# add variant
plot_data <- merge(plot_data, ddr, by = c("Sample", "Type"), all.x = TRUE)
plot_data$Variant <- as.character(plot_data$Variant)
plot_data[is.na(plot_data$Variant), "Variant"] <- "WT"

create.densityplot(
  list(
    ov_wt = plot_data[plot_data$Type == "OV" & plot_data$Variant == "WT", "GD"],
    ov_brca1 = plot_data[plot_data$Type == "OV" & plot_data$Variant == "BRCA1", "GD"],
    ov_brca2 = plot_data[plot_data$Type == "OV" & plot_data$Variant == "BRCA2", "GD"],
    brca_wt = plot_data[plot_data$Type == "BRCA" & plot_data$Variant == "WT", "GD"],
    brca_brca1 = plot_data[plot_data$Type == "BRCA" & plot_data$Variant == "BRCA1", "GD"],
    brca_brca2 = plot_data[plot_data$Type == "BRCA" & plot_data$Variant == "BRCA2", "GD"]
  ),
  filename = paste0("./output/figures/tcga/", date, "_genome_doubling_densityplot.tiff"),
  resolution = 300,
  col = default.colours(6),
  legend = list(
    inside = list(
      fun = draw.key,
      args = list(
        key = list(
          points = list(
            col = default.colours(6),
            pch = 21,
            cex = 1.5,
            fill = default.colours(6)
          ),
          text = list(
            lab = c("OV|WT", "OV|BRCA1", "OV|BRCA2", "BRCA|WT", "BRCA|BRCA1", "BRCA|BRCA2")
          ),
          title = "Cancer|Variant",
          padding = 2,
          cex = 1
        )
      ),
      x = 0.65,
      y = 0.97,
      draw = FALSE
    )
  ),
)
