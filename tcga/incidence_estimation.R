### DESCRIPTION ###################################################################################
# Kuan's driver estimates

### USAGE #########################################################################################
# ensure in same directory as mutation and annotation files 

### PREAMBLE ######################################################################################
# load libraries
library(argparse);
library(dplyr);

date <- Sys.Date();

source('~/svn/Collaborators/KuanHuang/helper_functions.R')

setwd('~/GermlineSomaticAssociations/genome-wide/output/somatic_gwas/pancan/driver_rate_comparison')
### OBTAIN COMMAND LINE ARGUMENTS #################################################################
parser <- ArgumentParser();

parser$add_argument('-e', '--gene', type = 'character', help = 'gene');
parser$add_argument('-c', '--cancer', type = 'character', help = 'cancer type');
parser$add_argument('-m', '--mutation', type = 'character', help = 'mutation type');

args <- parser$parse_args();

### DATA PROCESSING ###############################################################################
# read in variant annotation 
var_anno <- read.delim(
	'PCA_pathVar_integrated_filtered_adjusted.tsv',
	as.is = TRUE
	)
# find samples with variants in specified cancer
ddr <- var_anno[var_anno$HUGO_Symbol == args$gene & var_anno$cancer == args$cancer,'bcr_patient_barcode']

# read in patient annotation 
pat_anno <- read.delim(
	'clinical_PANCAN_patient_with_followup.tsv',
	as.is = TRUE
	)
# subset down to cancer type of interest
can_anno <- pat_anno[pat_anno$acronym == args$cancer,]
# only consider individuals of European descent 
can_anno <- can_anno[which(can_anno$race == 'WHITE'),]

# get mutation rate 
mut_rate <- get_mutation_rate(type = args$mutation, anno = can_anno)

# find samples without ddr 
wt <- can_anno$bcr_patient_barcode[!can_anno$bcr_patient_barcode %in% var_anno$bcr_patient_barcode]

# only keep samples with mutation rate 
ddr <- ddr[ddr %in% mut_rate$bcr_patient_barcode]
wt 	<- wt[wt %in% mut_rate$bcr_patient_barcode]

# print the number of ddr samples 
print(paste0("Number of ", args$cancer, " samples with variant in ", args$gene, ": ", length(ddr)))

# run bootstrap 
results <- do.call(rbind, sapply(
	1:10000,
	calculate_mutation_rate_ratio,
	mut_rate = mut_rate,
	ddr = ddr,
	wt = wt,
	anno = can_anno,
	cancer = args$cancer,
	simplify = FALSE
	))

# generate file name 
filename <- paste(date, args$cancer, args$gene, args$mutation, 'incidence_estimates.tsv', sep = '_')

# write to file 
write.table(
	results,
	file = filename,
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)
