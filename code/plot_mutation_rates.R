### DESCRIPTION ###################################################################################
# create somatic mutation rate plots

### USAGE #########################################################################################
# ensure in same directory as mutation and annotation files 

### PREAMBLE ######################################################################################
# load libraries
library(BoutrosLab.plotting.general)
library(dplyr)

date <- Sys.Date();

source('~/svn/Collaborators/KuanHuang/helper_functions.R')

setwd('~/GermlineSomaticAssociations/genome-wide/output/somatic_gwas/pancan/driver_rate_comparison')
### MUTATION RATE CNA AND SNV #####################################################################
calculate_mut_rate_cna_snv_meth_clones <- function(anno) {
	mut_rate_snv 			<- get_mutation_rate(type = 'snv', anno = anno)
	mut_rate_snv$mutation 	<- 'SNV'
	mut_rate_cna 			<- get_mutation_rate(type = 'cna', anno = anno)
	mut_rate_cna$mutation 	<- 'CNA'
	mut_rate_del			<- get_mutation_rate(type = 'deletion', anno = anno)
	mut_rate_del$mutation	<- 'CNA DEL'
	mut_rate_indel 			<- get_mutation_rate(type = 'indel', anno = anno)
	mut_rate_indel$mutation <- 'INDEL DEL'
	mut_rate_amp			<- get_mutation_rate(type = 'amplification', anno = anno)
	mut_rate_amp$mutation	<- 'CNA AMP'
	# mut_rate_meth 			<- get_mutation_rate(type = 'meth', anno = anno)
	# mut_rate_meth$mutation 	<- 'METH'
	# mut_rate_clones 		<- get_mutation_rate(type = 'clones', anno = anno)
	# mut_rate_clones$mutation <- 'CLONES'
	#mut_rate 				<- rbind(mut_rate_snv, mut_rate_cna, mut_rate_meth, mut_rate_clones)
	mut_rate 				<- rbind(mut_rate_snv, mut_rate_cna, mut_rate_del, mut_rate_indel, mut_rate_amp)
	return(mut_rate)
}

### FIND BRCA1 and BRCA2 MUTATED SAMPLES ##########################################################
find_brca1_brca2_samples <- function(cancer) {
	ddr <- list()
	ddr[['brca1']] 	<- var_anno[var_anno$HUGO_Symbol == 'BRCA1' & var_anno$cancer == cancer,'bcr_patient_barcode']
	ddr[['brca2']] 	<- var_anno[var_anno$HUGO_Symbol == 'BRCA2' & var_anno$cancer == cancer,'bcr_patient_barcode']	
	ddr[['both']] 	<- intersect(ddr$brca1, ddr$brca2)
	return(ddr)
}

### GET KRUSKAL P-VALUE ###########################################################################
get_kruskal_pvalue <- function(plot_data, group, value) {
	pvalue <- kruskal.test(
		plot_data[plot_data$group == group,value], 
		factor(plot_data[plot_data$group == group,'variant']))$p.value
	# convert to scientific notation 
	sci <- scientific.notation(pvalue, digits = 2)
	return(sci)
	}

### PLOT MUTATION RATES ###########################################################################
# read in patient annotation 
pat_anno <- read.delim(
	'clinical_PANCAN_patient_with_followup.tsv',
	as.is = TRUE
	)
# subset down to cancer type of interest
brca_anno 	<- pat_anno[pat_anno$acronym == 'BRCA',]
ov_anno		<- pat_anno[pat_anno$acronym == 'OV',]

# calculate mutation rate in breast cancer samples
mut_rate_brca 	<- calculate_mut_rate_cna_snv_meth_clones(brca_anno)
mut_rate_ov 	<- calculate_mut_rate_cna_snv_meth_clones(ov_anno)

# read in variant annotation 
var_anno <- read.delim(
	'PCA_pathVar_integrated_filtered_adjusted.tsv',
	as.is = TRUE
	)
# find samples with variants in BRCA2 and BRCA2 in BRCA
ddr_brca <- find_brca1_brca2_samples('BRCA')
# find samples with variants in BRCA1 and BRCA2 in OV
ddr_ov <- find_brca1_brca2_samples('OV')

# create plot data 
mut_rate_brca$cancer 	<- 'BRCA'
mut_rate_ov$cancer 		<- 'OV'
plot_data <- rbind(mut_rate_brca, mut_rate_ov)

# add variant status 
plot_data$variant <- 'WT'
plot_data[plot_data$bcr_patient_barcode %in% c(ddr_brca$brca1, ddr_ov$brca1),'variant'] <- 'BRCA1'
plot_data[plot_data$bcr_patient_barcode %in% c(ddr_brca$brca2, ddr_ov$brca2),'variant'] <- 'BRCA2'
plot_data[plot_data$bcr_patient_barcode %in% c(ddr_brca$both, ddr_ov$both),'variant'] <- 'BOTH'

# remove samples with both BRCA1 and BRCA2 
plot_data <- plot_data[plot_data$variant != 'BOTH',]

# merge mutation and cancer into a group 
plot_data$group <- paste(plot_data$cancer, plot_data$mutation)

# conver rate to log scale 
plot_data$log10count <- plot_data$count
plot_data[plot_data$mutation != 'METH', 'log10count'] <- log10(plot_data[plot_data$mutation != 'METH', 'log10count'])

boxplots <- list()
sci <- list()
for (group in c('BRCA SNV','BRCA CNA','BRCA CNA DEL','BRCA INDEL DEL','OV SNV','OV CNA','OV CNA DEL','OV INDEL DEL')) {
#for (group in c('BRCA SNV','BRCA CNA','BRCA METH','BRCA CLONES','OV SNV','OV CNA','OV METH','OV CLONES')) {
	rate_label <- paste(group, 'rate', sep = '_')
	count_label <- paste(group, 'count', sep = '_')
	# calculate p-value
	sci[[rate_label]] <- get_kruskal_pvalue(plot_data, group, value = 'rate')
	sci[[count_label]] <- get_kruskal_pvalue(plot_data, group, value = 'count')
	
	# skip ylab label for some plots
	if (grepl('CNA|DEL', group)) {
	#if (grepl('CNA|METH|CLONES', group)) {
		ylab_label_rate <- ''
		ylab_label_count <- ''
	} else {
		ylab_label_rate <- 'Mutation Rate'
		#ylab_label_count <- bquote('log'[10]*' Count')
		ylab_label_count <- 'Count'
		#ylab.label <- bquote('log'[10]*' mutation rate')
	}

	if (grepl('METH', group)) {
		ylimits_rate <- c(0, 0.005)
		yat_rate <- seq(0, 0.005, 0.001)
		ylimits_count <- c(0, 0.3)
		yat_count <- seq(0, 0.3, 0.05)
	} else if (grepl('CLONES', group)) {
		ylimits_rate <- c(0,0.35)
		yat_rate <- seq(0,0.35,0.1)
		ylimits_count <- c(0,16)
		yat_count <- seq(0,16,8)
	} else if (grepl('DEL', group)) {
		ylimits_rate <- c(0,4)
		yat_rate <- 0:4
		ylimits_count <- c(0,250)
		yat_count <- seq(0,250,50)
	} else {
		ylimits_rate <- c(0,10)
		yat_rate <- seq(0,10,2)
		ylimits_count <- c(0,400)
		yat_count <- seq(0,400,200)
	}

	# create boxplot 
	boxplots[[rate_label]] <- create.boxplot(
		rate ~ variant,
		add.stripplot = TRUE,
		data = plot_data[plot_data$group == group,], 
		#xaxis.rot = 90,
		ylimits = ylimits_rate,
		yat = yat_rate,
		xlab.label = '',
		xaxis.cex = 1.5,
		main = group,
		main.cex = 1.5,
		yaxis.cex = 1.5,
		ylab.cex = 1.5,
		ylab.label = ylab_label_rate,
		xaxis.lab = rep('',3),
		resolution = 300,
		key = list(
			text = list(
				lab = sci[[rate_label]]
				),
			x = 0.02,
			y = 0.98,
			cex = 1.5,
			padding.text = 3,
			corner = c(0,1)
			)
		)
	# create boxplot of counts 
	#ymax <- max(plot_data[plot_data$group == group,'log10count'])
	boxplots[[count_label]] <- create.boxplot(
		count ~ variant,
		add.stripplot = TRUE,
		data = plot_data[plot_data$group == group,], 
		#xaxis.rot = 90,
		ylimits = ylimits_count,
		yat = yat_count,
		#ylimits = c(0,ymax+ymax*0.2),
		#yat = seq(0,ymax+ymax*0.2,1),
		xlab.label = '',
		xaxis.cex = 1.5,
		main = '',
		main.cex = 1.5,
		yaxis.cex = 1.5,
		ylab.cex = 1.5,
		ylab.label = ylab_label_count,
		xaxis.lab = c('BRCA1','BRCA2','WT'),
		resolution = 300,
		key = list(
			text = list(
				lab = sci[[count_label]]
				),
			x = 0.02,
			y = 0.98,
			cex = 1.5,
			padding.text = 3,
			corner = c(0,1)
			)
		)

}

create.multipanelplot(
	boxplots[c(1,3,5,7,2,4,6,8)],
	filename = paste0(date, '_BRCA_mutation_rate_boxplot.tiff'),
	layout.width = 4,
	layout.height = 2,
	resolution = 300,
	width = 18
	)

create.multipanelplot(
	boxplots[c(9,11,13,15,10,12,14,16)],
	filename = paste0(date, '_OV_mutation_rate_boxplot.tiff'),
	layout.width = 4,
	layout.height = 2,
	resolution = 300,
	width = 18
	)


