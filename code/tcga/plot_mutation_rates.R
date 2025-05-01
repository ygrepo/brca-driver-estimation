### DESCRIPTION ###################################################################################
# create somatic mutation rate plots

### USAGE #########################################################################################
# ensure in same directory as mutation and annotation files

### PREAMBLE ######################################################################################
# load libraries
library(BoutrosLab.plotting.general)
library(dplyr)

rm(list = ls())
date <- Sys.Date()
setwd("/Users/yvesgreatti/github/brca-driver-estimation")

source("code/tcga/helper_functions.R")

# date <- Sys.Date();
#
# source('~/svn/Collaborators/KuanHuang/helper_functions.R')
#
# setwd('~/GermlineSomaticAssociations/genome-wide/output/somatic_gwas/pancan/driver_rate_comparison')
### PLOT MUTATION RATES ###########################################################################
# read in patient annotation
pat_anno <- read.delim(
  "data/TCGA/clinical_PANCAN_patient_with_followup.tsv",
  as.is = TRUE
)
# only consider individuals of European ancestry
pat_anno <- pat_anno[pat_anno$race == "WHITE", ]
# subset down to cancer type of interest
brca_anno <- pat_anno[pat_anno$acronym == "BRCA", ]
ov_anno <- pat_anno[pat_anno$acronym == "OV", ]

# calculate mutation rate in breast cancer samples
mut_rate_brca <- calculate_mut_rate_cna_snv_meth_clones(brca_anno)
mut_rate_ov <- calculate_mut_rate_cna_snv_meth_clones(ov_anno)

# add subtype and stage
mut_rate_brca <- add_subtype_and_stage(mut_rate_brca, pat_anno,
  grade = FALSE
)
mut_rate_ov <- add_subtype_and_stage(mut_rate_ov, pat_anno,
  subtype = FALSE, stage.type = "clinical"
)

# read in variant annotation
var_anno <- read.delim(
  #' data/TCGA_PanCanAtlas_2018/Germline_Huang_Cell2018/TCGA_ancestry_PC.txt',
  "data/TCGA/PCA_pathVar_integrated_filtered_adjusted_ancestry.tsv",
  as.is = TRUE
)

# var_anno <- read.delim(
# 	'PCA_pathVar_integrated_filtered_adjusted.tsv',
# 	as.is = TRUE
# 	)
# find samples with variants in BRCA2 and BRCA2 in BRCA
ddr_brca <- find_brca1_brca2_samples("BRCA")
# find samples with variants in BRCA1 and BRCA2 in OV
ddr_ov <- find_brca1_brca2_samples("OV")

# create plot data
mut_rate_brca$cancer <- "BRCA"
mut_rate_ov$cancer <- "OV"
plot_data <- rbind(mut_rate_brca, mut_rate_ov)

# convert stages to numeric
plot_data <- reformat_stage(plot_data)
# convert grade to numeric
plot_data$grade <- as.numeric(gsub("G", "", plot_data$grade))

# add variant status
plot_data$variant <- "WT"
plot_data[plot_data$bcr_patient_barcode %in% c(ddr_brca$brca1, ddr_ov$brca1), "variant"] <- "BRCA1"
plot_data[plot_data$bcr_patient_barcode %in% c(ddr_brca$brca2, ddr_ov$brca2), "variant"] <- "BRCA2"
plot_data[plot_data$bcr_patient_barcode %in% c(ddr_brca$both, ddr_ov$both), "variant"] <- "BOTH"

# remove samples with both BRCA1 and BRCA2
plot_data <- plot_data[plot_data$variant != "BOTH", ]

# merge mutation and cancer into a group
plot_data$group <- paste(plot_data$cancer, plot_data$mutation)

# conver rate to log scale
plot_data$log10count <- log10(plot_data$count)
plot_data$log10rate <- log10(plot_data$rate)
# plot_data[plot_data$mutation != 'METH', 'log10count'] <- log10(plot_data[plot_data$mutation != 'METH', 'log10count'])

groups <- c(
  "BRCA SNV", "BRCA CNA", "BRCA CNA DEL", "BRCA INDEL DEL",
  "BRCA CNA GAIN", "BRCA INSERT", "OV SNV", "OV CNA", "OV CNA DEL",
  "OV INDEL DEL", "OV CNA GAIN", "OV INSERT"
)

# run linear regression
params <- expand.grid(group = groups, variant = c("BRCA1", "BRCA2"), value = c("rate", "count"))
statsdf <- as.data.frame(t(apply(
  params,
  1,
  function(x) {
    fit_linear_reg(
      plot_data = plot_data, group = x["group"],
      value = x["value"], variant = x["variant"]
    )
  }
)))
statsdf <- cbind(params, statsdf)
statsdf$fdr <- p.adjust(statsdf$pvalue, method = "fdr")

boxplots <- list()

for (group in groups) {
  print(paste("Plotting", group))

  rate_label <- paste(group, "rate", sep = "_")
  count_label <- paste(group, "count", sep = "_")

  # Set Y-axis labels
  if (grepl("CNA$|DEL|INSERT", group)) {
    ylab_label_rate <- ""
    ylab_label_count <- ""
  } else {
    ylab_label_rate <- "Bps Mutated/Age at diagnosis\n(bps mutated/year)"
    ylab_label_count <- "Bps Mutated"
  }

  # Define axis parameters per group pattern
  yparams <- switch(TRUE,
    grepl("METH", group) ~ list(rate = c(0, 0.005), rate_at = seq(0, 0.005, 0.001), count = c(0, 0.3), count_at = seq(0, 0.3, 0.05)),
    grepl("CLONES", group) ~ list(rate = c(0, 0.35), rate_at = seq(0, 0.35, 0.1), count = c(0, 16), count_at = seq(0, 16, 8)),
    grepl("INDEL", group) ~ list(rate = c(-2.4, 4), rate_at = seq(-2, 3, 1), count = c(0, 5.5), count_at = 0:5),
    grepl("DEL", group) ~ list(rate = c(0, 4.5e7), rate_at = seq(0, 4e8, 1e7), count = c(0, 2.5e9), count_at = seq(0, 2e9, 1e9)),
    grepl("BRCA SNV", group) ~ list(rate = c(0, 100), rate_at = seq(0, 100, 20), count = c(0, 5500), count_at = seq(0, 5000, 1000)),
    grepl("OV SNV", group) ~ list(rate = c(0, 40), rate_at = seq(0, 40, 10), count = c(0, 2000), count_at = seq(0, 2000, 1000)),
    grepl("INSERT", group) ~ list(rate = c(0, 1.1), rate_at = seq(0, 1, 0.2), count = c(0, 80), count_at = seq(0, 80, 20)),
    grepl("CNA GAIN", group) ~ list(rate = c(0, 1.5e8), rate_at = seq(0, 1.5e8, 4e7), count = c(0, 4.8e9), count_at = seq(0, 4e9, 1e9)),
    grepl("CNA$", group) ~ list(rate = c(0, 1.5e8), rate_at = seq(0, 1.5e8, 4e7), count = c(0, 4.8e9), count_at = seq(0, 4e9, 1e9)),
    TRUE ~ list(rate = c(0, 1.5e8), rate_at = seq(0, 1.5e8, 4e7), count = c(0, 3.1e9), count_at = seq(0, 3e9, 1e9))
  )

  # Choose variable names based on INDEL logic
  is_indel <- grepl("INDEL", group)
  rate_var <- if (is_indel) "log10rate" else "rate"
  count_var <- if (is_indel) "log10count" else "count"

  # Subset data
  group_data <- plot_data[plot_data$group == group, ]

  # Plot rate boxplot
  boxplots[[rate_label]] <- create.boxplot(
    as.formula(paste(rate_var, "~ variant")),
    data = group_data,
    add.stripplot = TRUE,
    ylimits = yparams$rate,
    yat = yparams$rate_at,
    xlab.label = "",
    ylab.label = ylab_label_rate,
    xaxis.cex = 1.5,
    main = group,
    main.cex = 1.5,
    yaxis.cex = 1.5,
    ylab.cex = 1.5,
    xaxis.lab = rep("", 3),
    resolution = 300,
    key = add_linear_reg_key(statsdf, group, "rate", y = 0.98)
  )

  # Plot count boxplot
  boxplots[[count_label]] <- create.boxplot(
    as.formula(paste(count_var, "~ variant")),
    data = group_data,
    add.stripplot = TRUE,
    ylimits = yparams$count,
    yat = yparams$count_at,
    xlab.label = "",
    ylab.label = ylab_label_count,
    xaxis.cex = 1.5,
    main = "",
    main.cex = 1.5,
    yaxis.cex = 1.5,
    ylab.cex = 1.5,
    xaxis.lab = c("BRCA1", "BRCA2", "WT"),
    resolution = 300,
    key = add_linear_reg_key(statsdf, group, "count", y = 0.98)
  )
}

create.multipanelplot(
  boxplots[c(1, 3, 5, 7, 2, 4, 6, 8)],
  filename = paste0("./output/figures/tcga/", date, "_BRCA_mutation_rate_boxplot.tiff"),
  plot.objects.width = c(1.1, 1, 1, 1),
  layout.width = 4,
  layout.height = 2,
  resolution = 300,
  width = 21
)


create.multipanelplot(
  boxplots[c(9, 11, 10, 12)],
  filename = paste0("./output/figures/tcga/", date, "_BRCA_gains_mutation_rate_boxplot.tiff"),
  plot.objects.widths = c(0.55, 0.45),
  layout.width = 2,
  layout.height = 2,
  resolution = 300,
  width = 12
)

create.multipanelplot(
  boxplots[c(13, 15, 17, 19, 14, 16, 18, 20)],
  filename = paste0("./output/figures/tcga/", date, "_OV_mutation_rate_boxplot.tiff"),
  plot.objects.width = c(1.1, 1, 1, 1),
  layout.width = 4,
  layout.height = 2,
  resolution = 300,
  width = 21
)

create.multipanelplot(
  boxplots[c(21, 23, 22, 24)],
  filename = paste0("./output/figures/tcga/", date, "_OV_gains_mutation_rate_boxplot.tiff"),
  plot.objects.widths = c(0.55, 0.45),
  layout.width = 2,
  layout.height = 2,
  resolution = 300,
  width = 12
)
