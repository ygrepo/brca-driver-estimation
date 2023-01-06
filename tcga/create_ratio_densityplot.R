### DESCRIPTION ###################################################################################
# create density plot of ratio bootstrapped distributions

### USAGE #########################################################################################
# 

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)

setwd("/.mounts/labs/cpcgene/private/projects/GermlineSomaticAssociations/genome-wide/output/somatic_gwas/pancan/driver_rate_comparison/manuscript_figures/")

date <- Sys.Date()
### READ IN BOOTSTRAP RESULTS #####################################################################
read_in_bootstrap <- function(mutation, variant, cancer) {
	# if (mutation == 'snv') {
	# 	date <- '2020-01-02'
	# } else {
	# 	date <- '2019-12-24'
	# }
	date <- '2020-02-10'
	filename <- paste(date, cancer, variant, mutation, 'incidence_estimates.tsv', sep = '_')
	# read in 
	tmp <- read.delim(filename, as.is = TRUE)
	return(tmp)
}

### DATA PROCESSING ###############################################################################
# read in boostrap results for each mutation variant cancer triplet 
triplets <- expand.grid(mutation = c('snv','cna', 'deletion', 'amplification','indel'), variant = c('BRCA1','BRCA2'), cancer = c('BRCA','OV'))
bootstrap <- apply(
	triplets, 
	1,
	function(x) {
		read_in_bootstrap(mutation = x['mutation'], variant = x['variant'], cancer = x['cancer'])
		})

# reformat plot data
plot_data <- lapply(bootstrap, '[[', 5)
names(plot_data) <- paste(triplets$mutation, triplets$variant, triplets$cancer, sep = '|')

for (mut in c('snv','cna', 'deletion', 'amplification','indel')) {
	# create densityplot
	create.densityplot(
		plot_data[grep(mut, names(plot_data))],
		# list(brca1_brca = log10(plot_data[[3]]), 
		# 	brca2_brca = log10(plot_data[[7]]), 
		# 	brca1_ov = log10(plot_data[[11]]), 
		# 	brca2_ov = log10(plot_data[[15]])),
		filename = paste0(date, '_mutation_ratio_densityplot_', mut, '.tiff'),
		col = default.colours(4),
		legend = list(
	        inside = list(
	            fun = draw.key,
	            args = list(
	                key = list(
	                    points = list(
	                        col = default.colours(4),
	                        pch = 21,
	                        cex = 1.5,
	                        fill = default.colours(4)
	                        ),
	                    text = list(
	                        lab = gsub(paste0(mut,'|'), '', grep(mut, names(plot_data), value = TRUE), fixed = TRUE)
	                        ),
	                    title = 'Gene|Cancer',
	                    padding = 2,
	                    cex = 1
	                    )
	                ),
	            x = 0.65,
	            y = 0.97,
	            draw = FALSE
	            )
	        ),
		xlab.label = 'Mutation Rate Ratio',
		#xlab.label = expression('log'[10]*' Mutation Rate Ratio'),
		ylab.label = 'Density',
		resolution = 300
		)
}

### CALCULATE MEDIAN AND 95 CI ####################################################################
cis <- list()
for (i in 1:length(plot_data)) {
	tmp <- plot_data[[i]]
	tmp <- tmp[order(tmp)]
	l95 <- tmp[length(tmp)*0.025]
	u95 <- tmp[length(tmp)*0.975]
	cis[[i]] <- data.frame(
		id = names(plot_data)[i],
		median = median(tmp),
		l95 = l95,
		u95 = u95
		)
}
cis <- do.call(rbind, cis)
