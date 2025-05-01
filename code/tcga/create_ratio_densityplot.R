### DESCRIPTION ###################################################################################
# create density plot of ratio bootstrapped distributions

### USAGE #########################################################################################
# 

### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(dplyr)
rm(list = ls())

setwd('/Users/yvesgreatti/github/brca-driver-estimation')
# Plot one file ----
#filename <- paste0('data/TCGA/output/2025-04-30_BRCA_BRCA1_snv_incidence_estimates.tsv')
filename <- paste0('data/TCGA/incidence_estimates/2019-06-24_BRCA_BRCA1_snv_incidence_estimates.tsv')
filename

data <- read.delim(filename, as.is = TRUE)

filename <- paste0('data/TCGA/incidence_estimates/2019-06-24_mutation_ratio_densityplot_BRCA_BRCA1_snv.tsv')
data <- read.delim(filename, as.is = TRUE)
# Assuming the variable of interest is in the 5th column (adjust if different)
# Define label and extract the numeric 'ratio' column
label <- "BRCA1|BRCA"
data_list <- list()
data_list[[label]] <- data$incidence_two[!is.na(data$incidence_two)]  # Ensure no NA values

# Create the density plot
BoutrosLab.plotting.general::create.densityplot(
  data_list,
  filename = "output/figures/tcga/test_single_densityplot.png",
  col = default.colours(1),
  xlab.label = 'Mutation Rate Ratio',
  ylab.label = 'Density',
  legend = list(
    inside = list(
      fun = draw.key,
      args = list(
        key = list(
          points = list(
            col = default.colours(1),
            pch = 21,
            cex = 1.5,
            fill = default.colours(1)
          ),
          text = list(lab = label),
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
  resolution = 300
)

date <- Sys.Date()
### READ IN BOOTSTRAP RESULTS #####################################################################
# Define the function to map dates by mutation
get_date <- function(mutation, cancer, variant) {
  if (mutation %in% c('snv', 'cna')) {
    if  (variant == 'BRCA1') {
    #if (cancer == 'BRCA' && variant == 'BRCA1') {
        return('2019-06-24')
    } else {
      return('2019-06-25')
    }
  } else if (mutation == 'bp_deletions') {
    return('2019-10-29')
  } else if (mutation == 'bp_amplification') {
    return('2019-10-30')
  } else {
    stop("Unknown mutation type")
  }
}

# Read-in function
read_in_bootstrap <- function(mutation, variant, cancer) {
  date <- get_date(mutation, cancer, variant)
  filename <- paste(date, cancer, variant, mutation, 'incidence_estimates.tsv', sep = '_')
  filename <- paste0("data/TCGA/incidence_estimates/", filename)
  #filename <- file.path('incidence_estimates', filename)
  cat('Reading in:', filename, '\n')
  tmp <- read.delim(filename, as.is = TRUE)
  return(tmp)
}

# Generate triplets for all combinations
triplets <- expand.grid(
  mutation = c('snv', 'cna', 'bp_deletions', 'bp_amplification'),
  variant = c('BRCA1', 'BRCA2'),
  cancer = c('BRCA', 'OV'),
  stringsAsFactors = FALSE
)

# Read all bootstrap data
bootstrap <- apply(
  triplets,
  1,
  function(x) {
    read_in_bootstrap(mutation = x['mutation'], variant = x['variant'], cancer = x['cancer'])
  }
)

# 
# read_in_bootstrap <- function(mutation, variant, cancer) {
# 	# if (mutation == 'snv') {
# 	# 	date <- '2020-01-02'
# 	# } else {
# 	# 	date <- '2019-12-24'
# 	# }
#   if (mutation == 'bp_deletions' || mutation == 'bp_amplifications') {
#   	date <- '2019-10-29'
#   } else {
#   	date <- '2019-06-24'
#   }
#   #date <- '2020-02-10'
# 	filename <- paste(date, cancer, variant, mutation, 'incidence_estimates.tsv', sep = '_')
# 	filename <- paste0('incidence_estimates/', filename)
# 	cat('reading in ', filename, '\n')
# 	# read in 
# 	tmp <- read.delim(filename, as.is = TRUE)
# 	return(tmp)
# }
# 
# ### DATA PROCESSING ###############################################################################
# # read in boostrap results for each mutation variant cancer triplet 
# triplets <- expand.grid(mutation = c('snv','cna', 'bp_deletions', 'bp_amplifications'), variant = c('BRCA1','BRCA2'), cancer = c('BRCA','OV'))
# #triplets <- expand.grid(mutation = c('snv','cna', 'deletion', 'amplification','indel'), variant = c('BRCA1','BRCA2'), cancer = c('BRCA','OV'))
# bootstrap <- apply(
# 	triplets, 
# 	1,
# 	function(x) {
# 		read_in_bootstrap(mutation = x['mutation'], variant = x['variant'], cancer = x['cancer'])
# 		})

# reformat plot data
plot_data <- lapply(bootstrap, '[[', 5)
names(plot_data) <- paste(triplets$mutation, triplets$variant, triplets$cancer, sep = '|')

for (mut in c('snv','cna', 'bp_deletions', 'bp_amplification')) {
  fn = paste0(date, '_mutation_ratio_densityplot_', mut, '.tiff')
  fn = paste0('./output/figures/tcga/', fn, '.png')
#for (mut in c('snv','cna', 'deletion', 'amplification','indel')) {
    # create densityplot
  BoutrosLab.plotting.general::create.densityplot(
		plot_data[grep(mut, names(plot_data))],
		# list(brca1_brca = log10(plot_data[[3]]), 
		# 	brca2_brca = log10(plot_data[[7]]), 
		# 	brca1_ov = log10(plot_data[[11]]), 
		# 	brca2_ov = log10(plot_data[[15]])),
		#filename = paste0(date, '_mutation_ratio_densityplot_', mut, '.tiff'),
		filename = fn,
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
