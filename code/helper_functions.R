### GET MUTATION RATE #############################################################################
get_mutation_rate <- function(type, anno) {
	if (type == 'snv') {
	# read in somatic mutation rate
		# mut <- read.delim(
		# 	'TCGA_TMB.tsv',
		# 	as.is = TRUE
		# 	)
		mut <- read.delim(
			'TCGA_Tumor_Sample_patient_uniq_somatic_mutation_burden.tsv',
			as.is = TRUE
			)
		colnames(mut) <- gsub('nonsynonymous_count','count', colnames(mut))
	} else if (type == 'cna') {
		# read in cna segments
		# mut <- read.delim(
		# 	#'TCGA_cna_count_per_gene.tsv',
		# 	'TCGA_mastercalls.abs_segtabs.fixed.CNVperSample.tsv',
		# 	header = FALSE,
		# 	skip = 1,
		# 	col.names = c('count', 'bcr_patient_barcode'),
		# 	as.is = TRUE
		# 	)
		# mut$bcr_patient_barcode <- substr(mut$bcr_patient_barcode, 1, 12)
		mut <- read.delim(
			'TCGA_cna_deletion_bp.tsv',
			#'TCGA_segments_deletions_only.tsv',
			as.is = TRUE
			)
	} else if (type == 'deletion') {
		# read in cna segments, deletions only 
		mut <- read.delim(
			'TCGA_segments_deletions_only.tsv',
			as.is = TRUE
			)
	} else if (type == 'indel') {
		# read in cna segments 
		mut <- read.delim(
			'TCGA_indels_deletions_only.tsv',
			as.is = TRUE
			)
	} else if (type == 'amplification') {
		# read in cna segments 
		mut <- read.delim(
			'TCGA_cna_amplification_bp.tsv',
			as.is = TRUE
			)
	} else if (type == 'insertion') {
		# read in cna segments 
		mut <- read.delim(
			'TCGA_indels_insertions_only.tsv',
			as.is = TRUE
			)
	} else if (type == 'meth') {
		# read in methylation values 
		mut <- read.csv(
			'ciriello_methylation_rates.csv',
			as.is = TRUE,
			header = TRUE
			);
		mut <- mut[,c('SampleID','DMI_score')]
		mut$SampleID <- gsub('.','-',substr(mut$SampleID, 1, 12), fixed = TRUE)
		colnames(mut) <- c('bcr_patient_barcode','count')
	} else if (type == 'clones') {
		mut <- read.csv('ciriello_number_clones.csv', as.is = TRUE)
		mut <- mut[,c('sample_name','number.of.clones')]
		mut$sample_name <- substr(mut$sample_name, 1, 12)
		colnames(mut) <- c('bcr_patient_barcode','count')
	} else {
	 	stop("Invalid mutation type specified. Please specifiy either snv, meth or cna ...")
	}
	# if multiple samples per patient, find median mutation count 
	mut <- aggregate(mut$count, by = list(mut$bcr_patient_barcode), median, na.rm = TRUE)
	colnames(mut) <- c('bcr_patient_barcode','count')
	# convert number of mutations to mutation rate
	# divide by age at diagnosis 
	mut_rate <- merge(
		mut[,c('count','bcr_patient_barcode')], 
		anno[,c('bcr_patient_barcode','age_at_initial_pathologic_diagnosis')],
		by = 'bcr_patient_barcode'
		)
	mut_rate$rate <- mut_rate$count/as.numeric(mut_rate$age_at_initial_pathologic_diagnosis)
	rownames(mut_rate) <- mut_rate$bcr_patient_barcode
	return(mut_rate)
}

### CALCULATE MUTATION RATE RATIO #################################################################
calculate_mutation_rate_ratio <- function(int, mut_rate, ddr, wt) {
	# sample ddr and non ddr samples 
	ddr_sample <- sample(ddr, length(ddr), replace = TRUE)
	wt_sample <- sample(wt, length(ddr), replace = TRUE)

	# calculate median mutation rate of ddr and non ddr samples 
	ddr_median <- median(mut_rate[ddr_sample,'rate'])
	wt_median <- median(mut_rate[wt_sample,'rate'])

	# calculate ratio 
	ddr_ratio <- ddr_median/wt_median

	# calculate ratios 
	data.frame(
		int = int,
		num = length(ddr),
		wt_median = wt_median,
		ddr_median = ddr_median,
		ratio = ddr_ratio,
		incidence_two = ddr_ratio*ddr_ratio,
		incidence_three = ddr_ratio*ddr_ratio*(ddr_ratio^(1/2)),
		incidence_four = ddr_ratio*ddr_ratio*(ddr_ratio^(1/2))*(ddr_ratio^(1/3)),
		incidence_five = ddr_ratio*ddr_ratio*(ddr_ratio^(1/2))*(ddr_ratio^(1/3))*(ddr_ratio^(1/4)),
		incidence_six = ddr_ratio*ddr_ratio*(ddr_ratio^(1/2))*(ddr_ratio^(1/3))*(ddr_ratio^(1/4))*(ddr_ratio^(1/5)),
		incidence_seven = ddr_ratio*ddr_ratio*(ddr_ratio^(1/2))*(ddr_ratio^(1/3))*(ddr_ratio^(1/4))*(ddr_ratio^(1/5))*(ddr_ratio^(1/6)),
		incidence_eight = ddr_ratio*ddr_ratio*(ddr_ratio^(1/2))*(ddr_ratio^(1/3))*(ddr_ratio^(1/4))*(ddr_ratio^(1/5))*(ddr_ratio^(1/6))*(ddr_ratio^(1/7)),
		incidence_nine = ddr_ratio*ddr_ratio*(ddr_ratio^(1/2))*(ddr_ratio^(1/3))*(ddr_ratio^(1/4))*(ddr_ratio^(1/5))*(ddr_ratio^(1/6))*(ddr_ratio^(1/7))*(ddr_ratio^(1/8)),
		incidence_ten = ddr_ratio*ddr_ratio*(ddr_ratio^(1/2))*(ddr_ratio^(1/3))*(ddr_ratio^(1/4))*(ddr_ratio^(1/5))*(ddr_ratio^(1/6))*(ddr_ratio^(1/7))*(ddr_ratio^(1/8))*(ddr_ratio^(1/9)),
		incidence_eleven = ddr_ratio*ddr_ratio*(ddr_ratio^(1/2))*(ddr_ratio^(1/3))*(ddr_ratio^(1/4))*(ddr_ratio^(1/5))*(ddr_ratio^(1/6))*(ddr_ratio^(1/7))*(ddr_ratio^(1/8))*(ddr_ratio^(1/9))*(ddr_ratio^(1/10)),
		incidence_twelve = ddr_ratio*ddr_ratio*(ddr_ratio^(1/2))*(ddr_ratio^(1/3))*(ddr_ratio^(1/4))*(ddr_ratio^(1/5))*(ddr_ratio^(1/6))*(ddr_ratio^(1/7))*(ddr_ratio^(1/8))*(ddr_ratio^(1/9))*(ddr_ratio^(1/10))*(ddr_ratio^(1/11)),
		incidence_thirteen = ddr_ratio*ddr_ratio*(ddr_ratio^(1/2))*(ddr_ratio^(1/3))*(ddr_ratio^(1/4))*(ddr_ratio^(1/5))*(ddr_ratio^(1/6))*(ddr_ratio^(1/7))*(ddr_ratio^(1/8))*(ddr_ratio^(1/9))*(ddr_ratio^(1/10))*(ddr_ratio^(1/11))*(ddr_ratio^(1/12)),
		incidence_fourteen = ddr_ratio*ddr_ratio*(ddr_ratio^(1/2))*(ddr_ratio^(1/3))*(ddr_ratio^(1/4))*(ddr_ratio^(1/5))*(ddr_ratio^(1/6))*(ddr_ratio^(1/7))*(ddr_ratio^(1/8))*(ddr_ratio^(1/9))*(ddr_ratio^(1/10))*(ddr_ratio^(1/11))*(ddr_ratio^(1/12))*(ddr_ratio^(1/13)),
		incidence_fifteen = ddr_ratio*ddr_ratio*(ddr_ratio^(1/2))*(ddr_ratio^(1/3))*(ddr_ratio^(1/4))*(ddr_ratio^(1/5))*(ddr_ratio^(1/6))*(ddr_ratio^(1/7))*(ddr_ratio^(1/8))*(ddr_ratio^(1/9))*(ddr_ratio^(1/10))*(ddr_ratio^(1/11))*(ddr_ratio^(1/12))*(ddr_ratio^(1/13))*(ddr_ratio^(1/14)),
		incidence_sixteen = ddr_ratio*ddr_ratio*(ddr_ratio^(1/2))*(ddr_ratio^(1/3))*(ddr_ratio^(1/4))*(ddr_ratio^(1/5))*(ddr_ratio^(1/6))*(ddr_ratio^(1/7))*(ddr_ratio^(1/8))*(ddr_ratio^(1/9))*(ddr_ratio^(1/10))*(ddr_ratio^(1/11))*(ddr_ratio^(1/12))*(ddr_ratio^(1/13))*(ddr_ratio^(1/14))*(ddr_ratio^(1/15))
		)
}

### CALCULATE MEDIAN ESTIMATED INCIDENCES #########################################################
calculate_median_est_incidence <- function() {
	# create grid of cancer, gene and mutation type 
	par <- expand.grid(cancer = c('BRCA','OV'), gene = c('BRCA1','BRCA2'), mutation = c('snv','cna','meth', 'clones', 'indels','amplification','insertion'))
	# loop over grid and calculate medians 
	incidence <- apply(
		par[par$mutation == 'amplification',], 
		1,
		function(x) {
			filename <- paste(date, x['cancer'], x['gene'], x['mutation'], 'incidence_estimates.tsv', sep = '_') 
			print(filename)
			#filename <- list.files(pattern = fileflag)
			tmp <- read.delim(filename, as.is = TRUE)
			# read in observed 
			# # read in observed values 
			observed <- read.delim('observed_incidence_rates.tsv', as.is = TRUE)
			# extract median and CIs
			ob_median <- observed[observed$Cancer == x['cancer'] & observed$Gene == x['gene'],'Median']
			ob_L95 	<- observed[observed$Cancer == x['cancer'] & observed$Gene == x['gene'],'L95']
			ob_U95 	<- observed[observed$Cancer == x['cancer'] & observed$Gene == x['gene'],'U95']
			ob_n 	<- observed[observed$Cancer == x['cancer'] & observed$Gene == x['gene'],'Num']
			# calculate standard deviation 
			ob_sd <- sqrt(ob_n)*(ob_U95 - ob_L95)/3.92
			ob_dist <- rnorm(10000, mean = ob_median, sd = ob_sd)
			plots <- create_incidence_segplot(tmp, ob_median = ob_median, ob_L95, ob_U95, filename = paste(date, x['cancer'], x['gene'], x['mutation'], 'segplot.tiff', sep = '_'), 
				main = paste(x['cancer'], x['gene']))
			# calculate median 
			apply(tmp[,grep('incidence', colnames(tmp))], 2, mean)
			})
	colnames(incidence) <- apply(par, 1, paste, collapse = '_')
	rownames(incidence) <- c('two','three','four','five','six','seven','eight','nine','ten','eleven','twelve','thirteen')
	return(incidence)
}

create_incidence_segplot <- function(tmp, ob_median, ob_L95, ob_U95, filename, main, driver_max = NULL) {
	# create CI for bootstrap 
	mean_inc <- apply(tmp[,grep('incidence', colnames(tmp))], 2, mean)
	sd_inc <-  apply(tmp[,grep('incidence', colnames(tmp))], 2, sd)
	CI_inc <- 1.96*(sd_inc/sqrt(nrow(tmp)))

	plot_data <- data.frame(num = 2:(length(mean_inc)+1), mean = mean_inc, CI = CI_inc)
	plot_data$L95 <- plot_data$mean-plot_data$CI
	plot_data$U95 <- plot_data$mean+plot_data$CI
	# set maximum xlimits 
	if (!is.null(driver_max)) {
		plot_data <- plot_data[plot_data$num <= driver_max,]
	}
	plot_data$num <- factor(plot_data$num)
	max_y <- max(plot_data$mean+plot_data$CI)
	ylimits <- max(max_y+0.2*max_y, ob_U95+ob_U95*0.1)
	# create segplot
	create.scatterplot(
		mean ~ num,
		data = plot_data, 
		filename = filename,
		main = main,
		ylab.label = 'Incidence',
		xlab.label = 'Number Drivers',
		xaxis.lab = 2:(length(mean_inc)+1),
		ylimits = c(0, ylimits),
		yat = seq(0, ylimits, 10),
		#centers = plot_data$mean,
		y.error.up = plot_data$CI,
		resolution = 300,
		add.rectangle = TRUE,
		alpha.rectangle = 0.5,
		col.rectangle = 'grey',
		xleft.rectangle = 0,
		error.bar.lwd = 0.5,
		xright.rectangle = length(mean_inc)+1.5,
		ybottom.rectangle = ob_L95,
		ytop.rectangle = ob_U95,
		width = 9,
		#plot.horizontal = FALSE,
		abline.h = ob_median
		)
	}


#### CREATE BARPLOT OF MEDIANS ####################################################################
create_barplot_of_medians <- function(incidence) {
	# create barplot of medians
	incidence <- as.data.frame(incidence) 
	incidence$driver_number <- rownames(incidence)
	plot_data <- gather(incidence, 'group', 'incidence', -driver_number)

	observed <- read.delim('observed_incidence_rates.tsv', as.is = TRUE)
	observed$group <- paste(observed$Cancer, observed$Gene, sep = '_')
	# reformat observed data 
	ob_df <- data.frame(driver_number = 'observed', group = c(paste(observed$group, 'snv', sep = '_'), paste(observed$group, 'cna', sep = '_')),
		incidence = rep(observed$Median, 2))
	plot_data <- rbind(plot_data, ob_df)
	# make sure observed is last 
	plot_data$index <- NA
	plot_data[plot_data$driver_number == 'two','index'] <- 'a'
	plot_data[plot_data$driver_number == 'three','index'] <- 'b'
	plot_data[plot_data$driver_number == 'four','index'] <- 'c'
	plot_data[plot_data$driver_number == 'five','index'] <- 'd'
	plot_data[plot_data$driver_number == 'observed','index'] <- 'e'

	create.barplot(
		incidence ~ group, 
		data = plot_data, 
		filename = paste0(date, '_incidence_means_barplot.tiff'),
		groups = plot_data$index,
		resolution = 300,
		col = c(default.colours(4),'black'),
		#xaxis.rot = 90,
		ylab.label = 'Incidence',
		xaxis.lab = c('BRCA\nBRCA1\nCNA', 'BRCA\nBRCA1\nSNV', 'BRCA\nBRCA2\nCNA', 'BRCA\nBRCA2\nSNV',
			'OV\nBRCA1\nCNA', 'OV\nBRCA1\nSNV', 'OV\nBRCA2\nCNA', 'OV\nBRCA1\nSNV'),
		xaxis.cex = 1,
		xlab.label = '',
		width = 8,
		legend = list(
			inside = list(
				fun = draw.key,
				args = list(
					key = list(
						points = list(
							col = c(default.colours(4), 'black'),
							pch = 21,
							cex = 1.5,
							fill = c(default.colours(4), 'black')
							),
						text = list(
							lab = c('Two','Three','Four','Five', 'Observed')
							),
						cex = 1
						)
					),
				x = 0.75,
				y = 0.97,
				draw = FALSE
				)
			)
		)
}

### CALCULATE CONFIDENCE INTERVALS ################################################################
calculate_CI <- function(x, upper = TRUE) {
	stand_dev <- sd(x)
	med <- median(x)
	if (upper) {
		CI <- med + 1.96*(stand_dev/sqrt(length(x)))
	} else {
		CI <- med - 1.96*(stand_dev/sqrt(length(x)))
		}
	return(CI)
}