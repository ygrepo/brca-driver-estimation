### DESCRIPTION ###################################################################################
# Kuan's driver estimates

### USAGE #########################################################################################
# ensure in same directory as mutation and annotation files 

### PREAMBLE ######################################################################################
# load libraries
library(argparse);
library(dplyr);

date <- Sys.Date();

source('helper_functions.R')

setwd('/N/projects/curtis/METABRIC')
outdir <- "/oak/stanford/groups/ccurtis2/users/khoulaha/BRCA_drivers"
### OBTAIN COMMAND LINE ARGUMENTS #################################################################
parser <- ArgumentParser();

parser$add_argument('-e', '--gene', type = 'character', help = 'gene');
parser$add_argument('-m', '--mutation', type = 'character', help = 'mutation type');

args <- parser$parse_args();

### DATA PROCESSING ###############################################################################
# read in variant annotation 
var_anno <- read.delim(
	'targetSeq/allVariants_forChristina.txt',
	as.is = TRUE
	)
# find samples with variants in specified cancer
ddr <- var_anno[var_anno$gene == args$gene & var_anno$type == 'PATH', 'sample']

# read in samples to consider
samples <- read.delim(
	file.path(outdir, 'sample_list.txt'),
	header = FALSE,
 	as.is = TRUE
 	)
wt <- samples[!samples$V1 %in% ddr,1]

# read in patient annotation 
pat_anno <- read.delim(
	'cbioportal/brca_metabric/data_clinical_patient.txt',
	skip = 4,
	as.is = TRUE
	)
pat_anno <- pat_anno[pat_anno$PATIENT_ID %in% samples$V1,]

# add PAM50 
add_anno <- read.delim(
	'clinical/Clinical1980V2.txt',
	as.is = TRUE
	)
add_anno <- add_anno[,c('METABRIC_ID','pam50_official','stage')]
colnames(add_anno) <- c('PATIENT_ID','pam50','stage')
pat_anno <- merge(pat_anno, add_anno, by = 'PATIENT_ID')

# get mutation rate 
mut_rate <- get_mutation_rate(type = args$mutation, anno = pat_anno)

# only keep samples with mutation rate 
ddr <- ddr[ddr %in% mut_rate$sample]
wt 	<- wt[wt %in% mut_rate$sample]

# print the number of ddr samples 
print(paste0("Number of samples with variant in ", args$gene, ": ", length(ddr)))

# run bootstrap 
results <- do.call(rbind, sapply(
	1:10000,
	calculate_mutation_rate_ratio,
	mut_rate = mut_rate,
	ddr = ddr,
	wt = wt,
	anno = pat_anno,
	simplify = FALSE
	))

# generate file name 
filename <- file.path(outdir, paste(date, args$gene, args$mutation, 'incidence_estimates.tsv', sep = '_'))

# write to file 
write.table(
	results,
	file = filename,
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)
