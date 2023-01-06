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
### PLOT MUTATION RATES ###########################################################################
# read in patient annotation 
pat_anno2 <- read.delim(
	'clinical_PANCAN_patient_with_followup.tsv',
	as.is = TRUE
	)
# only consider individuals of European ancestry 
pat_anno <- pat_anno[pat_anno$race == 'WHITE',]
# subset down to cancer type of interest
brca_anno 	<- pat_anno[pat_anno$acronym == 'BRCA',]
ov_anno		<- pat_anno[pat_anno$acronym == 'OV',]

# calculate mutation rate in breast cancer samples
mut_rate_brca 	<- calculate_mut_rate_cna_snv_meth_clones(brca_anno)
mut_rate_ov 	<- calculate_mut_rate_cna_snv_meth_clones(ov_anno)

# add subtype and stage 
mut_rate_brca 	<- add_subtype_and_stage(mut_rate_brca, pat_anno, 
	grade = FALSE)
mut_rate_ov 	<- add_subtype_and_stage(mut_rate_ov, pat_anno, 
	subtype = FALSE, stage.type = 'clinical') 

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

# convert stages to numeric 
plot_data <- reformat_stage(plot_data)
# convert grade to numeric 
plot_data$grade <- as.numeric(gsub('G', '', plot_data$grade))

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
plot_data$log10count 	<- log10(plot_data$count)
plot_data$log10rate 	<- log10(plot_data$rate)
#plot_data[plot_data$mutation != 'METH', 'log10count'] <- log10(plot_data[plot_data$mutation != 'METH', 'log10count'])

groups <- c('BRCA SNV','BRCA CNA','BRCA CNA DEL','BRCA INDEL DEL',
	'BRCA CNA GAIN','BRCA INSERT','OV SNV','OV CNA','OV CNA DEL',
	'OV INDEL DEL','OV CNA GAIN', 'OV INSERT')

# run linear regression 
params <- expand.grid(group = groups, variant = c('BRCA1','BRCA2'), value = c('rate','count'))
statsdf <- as.data.frame(t(apply(
	params,
	1,
	function(x) {
		fit_linear_reg(plot_data = plot_data, group = x['group'],
			value = x['value'], variant = x['variant'])	
		}
	)))
statsdf <- cbind(params, statsdf)
statsdf$fdr <- p.adjust(statsdf$pvalue, method = 'fdr')

boxplots <- list()
#sci <- list()
for (group in groups) {
	print(paste("Plotting", group))
#for (group in c('BRCA SNV','BRCA CNA','BRCA METH','BRCA CLONES','OV SNV','OV CNA','OV METH','OV CLONES')) {
	rate_label <- paste(group, 'rate', sep = '_')
	count_label <- paste(group, 'count', sep = '_')
	# # calculate p-value
	# sci[[rate_label]] <- get_kruskal_pvalue(plot_data, group, value = 'rate')
	# sci[[count_label]] <- get_kruskal_pvalue(plot_data, group, value = 'count')
	
	# skip ylab label for some plots
	if (grepl('CNA$|DEL|INSERT', group)) {
	#if (grepl('CNA|METH|CLONES', group)) {
		ylab_label_rate <- ''
		ylab_label_count <- ''
	} else {
		ylab_label_rate <- 'Bps Mutated/Age at diagnosis\n(bps mutated/year)'
		#ylab_label_count <- bquote('log'[10]*' Count')
		ylab_label_count <- 'Bps Mutated'
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
	} else if (grepl('INDEL', group)) {
		ylimits_rate <- c(-2.4,4)
		yat_rate <- seq(-2,3,1)
		yaxis_rate <- c(expression(1%*%10^-2), expression(1%*%10^-1), expression(1), expression(10), expression(1%*%10^2),expression(1%*%10^3))
		ylimits_count <- c(0,5.5)
		yat_count <- 0:5
		yaxis_count <- c(expression(0), expression(1%*%10^1), expression(1%*%10^2), expression(1%*%10^3), expression(1%*%10^4), expression(1%*%10^5))
	} else if (grepl('DEL', group)) {
		ylimits_rate <- c(0,4.5e7)
		yat_rate <- seq(0,4e8,1e7)
		yaxis_rate <- c(expression(0), expression(1%*%10^7), expression(2%*%10^7), expression(3%*%10^7), expression(4%*%10^7))
		ylimits_count <- c(0,2.5e9)
		yat_count <- seq(0,2e9,1e9)
		yaxis_count <- c(expression(0), expression(1%*%10^9), expression(2%*%10^9))
	} else if (grepl('BRCA SNV', group)) {
		ylimits_rate <- c(0,100)
		yat_rate <- seq(0,100,20)
		yaxis_rate <- c(expression(0), expression(20), expression(40), expression(60), expression(80), expression(100))
		ylimits_count <- c(0,5500)
		yat_count <- seq(0,5000,1000)
		yaxis_count <- c(expression(0), expression(1%*%10^3), expression(2%*%10^3), expression(3%*%10^3), expression(4%*%10^3), expression(5%*%10^3))
	} else if (grepl('OV SNV', group)) {
		ylimits_rate <- c(0,40)
		yat_rate <- seq(0,40,10)
		yaxis_rate <- c(expression(0), expression(10), expression(20), expression(30), expression(40))
		ylimits_count <- c(0,2000)
		yat_count <- seq(0,2000,1000)
		yaxis_count <- c(expression(0), expression(1%*%10^3), expression(2%*%10^3))
	} else if (grepl('INSERT', group)) {
		ylimits_rate <- c(0,1.1)
		yat_rate <- seq(0,1,0.2)
		yaxis_rate <- seq(0,1,0.2)
		ylimits_count <- c(0,80)
		yat_count <- seq(0,80,20) 
		yaxis_count <- seq(0,80,20)
	} else if (grepl('CNA GAIN', group)) {
		ylimits_rate <- c(0,1.5e8)
		yat_rate <- seq(0,1.5e8,4e7)
		yaxis_rate <- c(expression(0.0), expression(2.0%*%10^7), expression(8.0%*%10^7), expression(1.2%*%10^8))
		ylimits_count <- c(0,4.8e9)
		yat_count <- seq(0,4e9,1e9)
		yaxis_count <- c(expression(0), expression(1%*%10^9), expression(2%*%10^9), expression(3%*%10^9), expression(4%*%10^9))
	} else if (grepl('CNA$', group)) {
		ylimits_rate <- c(0,1.5e8)
		yat_rate <- seq(0,1.5e8,4e7)
		yaxis_rate <- c(expression(0.0), expression(2.0%*%10^7), expression(8.0%*%10^7), expression(1.2%*%10^8))
		ylimits_count <- c(0,4.8e9)
		yat_count <- seq(0,4e9,1e9)
		yaxis_count <- c(expression(0), expression(1%*%10^9), expression(2%*%10^9), expression(3%*%10^9), expression(4%*%10^9))
	} else {
		ylimits_rate <- c(0,1.5e8)
		yat_rate <- seq(0,1.5e8,4e7)
		yaxis_rate <- c(expression(0.0), expression(2.0%*%10^7), expression(8.0%*%10^7), expression(1.2%*%10^8))
		ylimits_count <- c(0,3.1e9)
		yat_count <- seq(0,3e9,1e9)
		yaxis_count <- c(expression(0), expression(1%*%10^9), expression(2%*%10^9), expression(3%*%10^9))
	}


	if (!grepl('INDEL', group)) {
		# create boxplot 
		boxplots[[rate_label]] <- create.boxplot(
			rate ~ variant,
			add.stripplot = TRUE,
			data = plot_data[plot_data$group == group,], 
			#xaxis.rot = 90,
			ylimits = ylimits_rate,
			yat = yat_rate,
			yaxis.lab = yaxis_rate,
			xlab.label = '',
			xaxis.cex = 1.5,
			main = group,
			main.cex = 1.5,
			yaxis.cex = 1.5,
			ylab.cex = 1.5,
			ylab.label = ylab_label_rate,
			xaxis.lab = rep('',3),
			resolution = 300,
			key = add_linear_reg_key(
				statsdf = statsdf, 
				group = group,
				value = 'rate',
				y = 0.98
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
			yaxis.lab = yaxis_count,
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
			key = add_linear_reg_key(
				statsdf = statsdf, 
				group = group,
				value = 'count',
				y = 0.98
				)
			)
	} else {
		# create boxplot 
		boxplots[[rate_label]] <- create.boxplot(
			log10rate ~ variant,
			add.stripplot = TRUE,
			data = plot_data[plot_data$group == group,], 
			#xaxis.rot = 90,
			ylimits = ylimits_rate,
			yat = yat_rate,
			yaxis.lab = yaxis_rate,
			xlab.label = '',
			xaxis.cex = 1.5,
			main = group,
			main.cex = 1.5,
			yaxis.cex = 1.5,
			ylab.cex = 1.5,
			ylab.label = ylab_label_rate,
			xaxis.lab = rep('',3),
			resolution = 300,
			key = add_linear_reg_key(
				statsdf = statsdf, 
				group = group,
				value = 'rate',
				y = 0.98
				)
			)
		# create boxplot of counts 
		#ymax <- max(plot_data[plot_data$group == group,'log10count'])
		boxplots[[count_label]] <- create.boxplot(
			log10count ~ variant,
			add.stripplot = TRUE,
			data = plot_data[plot_data$group == group,], 
			#xaxis.rot = 90,
			ylimits = ylimits_count,
			yat = yat_count,
			yaxis.lab = yaxis_count,
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
			key = add_linear_reg_key(
				statsdf = statsdf, 
				group = group,
				value = 'count',
				y = 0.98
				)
			)

		}

}

create.multipanelplot(
	boxplots[c(1,3,5,7,2,4,6,8)],
	filename = paste0(date, '_BRCA_mutation_rate_boxplot.tiff'),
	plot.objects.width = c(1.1,1,1,1),
	layout.width = 4,
	layout.height = 2,
	resolution = 300,
	width = 21
	)


create.multipanelplot(
	boxplots[c(9,11,10,12)],
	filename = paste0(date, '_BRCA_gains_mutation_rate_boxplot.tiff'),
	plot.objects.widths = c(0.55,0.45),
	layout.width = 2,
	layout.height = 2,
	resolution = 300,
	width = 12
	)

create.multipanelplot(
	boxplots[c(13,15,17,19,14,16,18,20)],
	filename = paste0(date, '_OV_mutation_rate_boxplot.tiff'),
	plot.objects.width = c(1.1,1,1,1),
	layout.width = 4,
	layout.height = 2,
	resolution = 300,
	width = 21
	)

create.multipanelplot(
	boxplots[c(21,23,22,24)],
	filename = paste0(date, '_OV_gains_mutation_rate_boxplot.tiff'),
	plot.objects.widths = c(0.55,0.45),
	layout.width = 2,
	layout.height = 2,
	resolution = 300,
	width = 12
	)
