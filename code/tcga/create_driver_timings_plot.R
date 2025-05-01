### DESCRIPTION ###################################################################################
# create early drivers plot


### USAGE #########################################################################################


### PREAMBLE ######################################################################################
library(BoutrosLab.plotting.general)
library(dplyr)
rm(list = ls())

setwd('/Users/yvesgreatti/github/brca-driver-estimation')

date <- Sys.Date()
### COUNT NUMBER OF EARLY CLONAL ##################################################################
count_number_of_early_clonal <- function(df, ids) {
	# subset pcawg data to specified samples
	df <- df[df$sampleNames %in% ids,]
	# subset down to early clonal mutations
	df_early <- df[which(df$CLS == "clonal [early]"),]
	# categorize as snv, deletion and insertion 
	df_early$type <- NA
	df_early$type <- apply(
		df_early, 
		1,x
		function(x) {
			ref <- nchar(x['ref'])
			alt <- nchar(x['alt'])
			if (is.na(ref) | is.na(alt)) {
				return(NA)
			} else if (ref == alt) {
				return('SNV')
			} else if (ref > alt) {
				return('DEL')
			} else if (ref < alt) {
				return('INS')
			} else {
				return(NA)
				}
			})
	# filter out deplicate SNVs 
	df_early <- unique(df_early[,c('seqnames','start','type','ref','alt')])
	# create plot data 
	plot_data <- as.data.frame(table(df_early$type))
	colnames(plot_data) <- c('Mut','Count')
	return(plot_data)
	}

### SNVS AND INDELS ###############################################################################
# read in evolutionary timings data from pcawg 
pcawg <- read.delim(
	'2018-07-25-driversTiming.icgc.controlled.txt',
	as.is = TRUE
	)

# read in annotation 
anno <- read.delim(
	'~/isilon/private/Resources/data/ICGC_PanCancer/PCAWG/data/annotation/2019-11-26_DCC/sample.all_projects.tsv.gz',
	as.is = TRUE
	)
# find all breast samples 
breast_ids <- anno[grep('BRCA', anno$project_code),'icgc_donor_id']
ovarian_ids <- anno[grep('OV', anno$project_code),'icgc_donor_id']

# count number of early clonal mut types in breast cancer 
breast <- count_number_of_early_clonal(pcawg, breast_ids)
breast$cancer <- 'BRCA'

# count number of early clonal mut types in ovarian cancer 
ovarian <- count_number_of_early_clonal(pcawg, ovarian_ids)
ovarian$cancer <- 'OV'

plot_data <- rbind(breast, ovarian)

# create barplot 
create.barplot(
	Count ~ cancer,
	data = plot_data,
	filename = paste0(date, '_early_clonal_driver_count.tiff'),
	groups = plot_data$Mut,
	col = default.colours(3),
	ylab.label = 'Number Mutations',
	xlab.label = 'Cancer',
	ylimits = c(0,50),
	yat = seq(0,50,10),
	legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = list(
                    points = list(
                        col = 'black',
                        pch = 22,
                        cex = 3,
                        fill = default.colours(3)
                        ),
                    text = list(
                        lab = c('Deletion','Insertion','SNV')
                        ),
                    padding.text = 5,
                    cex = 1
                    )
                ),
            	title = 'Mutation',
                # Positioning legend on plot
                x = 0.01,
                y = 0.99
            )
        ),
	resolution = 300
	)

### SEGMENTS ######################################################################################
# read in segments 
seg <- read.delim(
	'2018-07-25-allSegmentsTime.icgc.controlled.txt', 
	as.is = TRUE
	)

# read in annotation 
anno <- read.delim(
	'icgc_sample_annotations_summary_table.txt',
	as.is = TRUE
	)
# find breast cancer ids 
breast_ids <- anno[anno$histology_abbreviation == 'Breast-AdenoCA', 'tumour_aliquot_id']
ovarian_ids <- anno[anno$histology_abbreviation == 'Ovary-AdenoCA', 'tumour_aliquot_id']

# subset down to just breast samples 
breast_seg <- seg[seg$samplename %in% breast_ids,]


tmp <- breast_seg[breast_seg$samplename == 'f393baf9-2710-9203-e040-11ac0d484504',]
sgain <- tmp[tmp$type == 'SingleGain',]
cnloh <- tmp[tmp$type == 'CNLOH',]
dgain <- tmp[tmp$type == 'DoubleGain',]

create.densityplot(
	list(sgain = sgain$time, cnloh = cnloh$time, dgain = dgain$time),
	filename = paste0(date, '_f393baf9-2710-9203-e040-11ac0d484504_segments_timing_densityplot.tiff'),
	col = default.colours(3),
	resolution = 300,
	legend = list(
		inside = list(
			fun = draw.key,
			args = list(
                key = list(
                    points = list(
                        col = default.colours(3),
                        pch = 21,
                        cex = 1.5,
                        fill = default.colours(3)
                        ),
                    text = list(
                        lab = c('CNLOH','Double Gain','Single Gain')
                        ),
                    cex = 1
                    )
                ),
            x = 0.01,
            y = 0.99,
            draw = FALSE
            )
        )
)

### CREATE LEAGUE PLOTTING DATA ###################################################################
create_league_plotting_data <- function(df) {
	# subset down to clonal 
	df <- df[which(df$assignment == 'clonal'),]
	# remove WGD 
	df <- df[df$event != 'WGD',]
	# categorize mutations 
	df$cat <- NA
	df[grep('loss|homdel', df$event),'cat'] <- 'Del'
	df[grep('gain', df$event),'cat'] <- 'Gain'
	df[is.na(df$cat),'cat'] <- "SNV"
	# create plot data 
	plot_data <- as.data.frame(table(df$cat))
	colnames(plot_data) <- c('Mut','Count')
	return(plot_data)
	}
 
### LEAGUE RESULTS ################################################################################
# read in league results 
league <- read.delim(
	'timelines_leagueModelEvents.txt', 
	as.is = TRUE
	)

# subset down to breast cancer
breast <- league[which(league$cohort == 'Breast-AdenoCA'),]
# create breast cancer plot data 
bplot_data <- create_league_plotting_data(breast)
bplot_data$Cancer <- 'BRCA'

# subset down to ovarian cancer 
ovarian <- league[which(league$cohort == 'Ovary-AdenoCA'),]
# create ovarian cancer plot data
oplot_data <- create_league_plotting_data(ovarian)
oplot_data$Cancer <- 'OV'

plot_data <- rbind(bplot_data, oplot_data)

# create barplot 
create.barplot(
	Count ~ Cancer,
	data = plot_data,
	filename = paste0(date, '_clonal_driver_count.tiff'),
	groups = plot_data$Mut,
	col = default.colours(3),
	ylab.label = 'Number Mutations',
	xlab.label = 'Cancer',
	ylimits = c(0,20),
	yat = seq(0,20,5),
	legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = list(
                    points = list(
                        col = 'black',
                        pch = 22,
                        cex = 3,
                        fill = default.colours(3)
                        ),
                    text = list(
                        lab = c('Deletion','Gain','SNV')
                        ),
                    padding.text = 5,
                    cex = 1
                    )
                ),
            	title = 'Mutation',
                # Positioning legend on plot
                x = 0.01,
                y = 0.99
            )
        ),
	resolution = 300
	)

# calculate hypergeometric stats 
calculate_hypergeo_stats <- function(df) {
	df <- df[!is.na(df$assignment),]
	q <- sum(df$type == 'deletion' & df$assignment == 'clonal')
	m <- sum(df$type == 'deletion')
	n <- nrow(df)-m
	k <- sum(df$assignment == 'clonal')
	1-phyper(q=q-1, m=m, n=n, k=k)
	}

### CREATE TIMING VS RECURRENCE PLOT ##############################################################
# calculate frequency 
# setting number of breast cancer samples from supplementary table 1
breast_num <- 198
breast$freq <- rowSums(breast[,c('num_clonal','num_subclonal')])/breast_num

# set mutation type 
breast$type <- NA
breast[grep('loss|homdel', breast$event),'type'] <- 'deletion'
breast[grep('gain', breast$event),'type'] <- 'gain'
breast[grep('WGD', breast$event),'type'] <- 'wgd'
breast[is.na(breast$type),'type'] <- 'snv'

# assignment ranks based on supplementary figure from pcawg evolutionary work 
breast_ranks <- data.frame(
	event = c('loss_17p','PIK3CA','CBFB','homdel_10q23.31','TP53','loss_19p13.3',
		'homdel_9p21.3','loss_3p','loss_13q','loss_16q','GATA3','gain_8q','gain_1q',
		'loss_8p','loss_11q','gain_11q13.3','gain_8q24.21','MAP3K1', 'loss_22q',
		'loss_11p15.5', 'loss_4p','gain_16p','gain_17q','loss_15q','loss_6q',
		'loss_2q37.3','WGD', 'gain_10p'),
	rank = 1:28
	)
breast <- merge(breast, breast_ranks, by = 'event')
breast <- breast[order(breast$rank),]


create.scatterplot(
	freq ~ rank,
	data = breast,
	groups = breast$type,
	col = c('dodgerblue','firebrick','orange','darkorchid'),
	filename = paste0(date, '_breast_driver_timing_scatterplot.tiff'),
	resolution = 300,
	ylab.label = 'Recurrence',
	xlab.label = 'Timing',
	ylimits = c(0,1.02),
	yat = seq(0,1,0.2),
	xlimits = c(0, 28.5),
	xat = c(0,28.5),
	xaxis.lab = c('Early','Late'),
	main = 'BRCA',
	add.text = TRUE,
	text.labels = gsub('loss_|homdel_','',c(breast$event[1:5], 'loss_8p','loss_13q','loss_16q','loss_11q')),
	#text.x = breast$rank[1:5]+0.5,
	text.x = c(1.5, 2.7, 3.5, 4.2, 5.5, 14.5, 9.5, 10.5, 15.5), 
	text.y = c(0.80, 0.29, 0.02, 0.13, 0.52, 0.72, 0.66, 0.57, 0.58),
	#text.y = breast$freq[1:5]+0.05,
	text.cex = 1.3,
	legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = list(
                    points = list(
                        col = 'black',
                        pch = 21,
                        cex = 2,
                        fill = c('dodgerblue','firebrick','orange','darkorchid')
                        ),
                    text = list(
                        lab = c('Deletion','Gain','SNV','WGD')
                        ),
                    padding.text = 5,
                    cex = 1
                    )
                ),
            	title = 'Mutation',
                # Positioning legend on plot
                x = 0.7,
                y = 0.99
            )
        )
	)


### OVARIAN ###
# calculate frequency 
# setting number of ovarian cancer samples from supplementary table 1
ovarian_num <- 113
ovarian$freq <- rowSums(ovarian[,c('num_clonal','num_subclonal')])/ovarian_num

# set mutation type 
ovarian$type <- NA
ovarian[grep('loss|homdel', ovarian$event),'type'] <- 'deletion'
ovarian[grep('gain', ovarian$event),'type'] <- 'gain'
ovarian[grep('WGD', ovarian$event),'type'] <- 'wgd'
ovarian[is.na(ovarian$type),'type'] <- 'snv'

# order by recurrence to plot labels 
ovarian_ranks <- data.frame(
	event = c("loss_17", "TP53", "loss_19p13.3", "loss_22q", "loss_13q",
		'LATS1','loss_4q', "gain_7q", "gain_3q26.2", "loss_11p15.5", 
		"gain_3q", "gain_20q", "gain_5p", "gain_8q24.21", "loss_16q",
		"loss_8p", "loss_18q", "loss_5q","CDK12", "WGD","loss_6q",
		"RB1","gain_17q","loss_14q", "POLR3E","gain_5p15.33"),
	rank = 1:26
	)
ovarian <- merge(ovarian, ovarian_ranks, by = 'event')
ovarian <- ovarian[order(ovarian$rank),]

create.scatterplot(
	freq ~ rank,
	data = ovarian,
	groups = ovarian$type,
	col = c('dodgerblue','firebrick','orange','darkorchid'),
	filename = paste0(date, '_ovarian_driver_timing_scatterplot.tiff'),
	resolution = 300,
	ylab.label = 'Recurrence',
	xlab.label = 'Timing',
	ylimits = c(0,1.02),
	yat = seq(0,1,0.2),
	xlimits = c(0, 28.5),
	xat = c(0,28.5),
	xaxis.lab = c('Early','Late'),
	main = 'OV',
	add.text = TRUE,
	text.labels = gsub('loss_','',c(ovarian$event[1:5], 'loss_4q')),
	text.x = c(1, 4.3, 3.5, 4.5, 5.1, 8.2),
	text.y = c(0.98, 0.94, 0.65, 0.75, 0.83, 0.82), #ovarian$freq[1:5]+0.03,
	text.cex = 1.3,
	legend = list(
        inside = list(
            fun = draw.key,
            args = list(
                key = list(
                    points = list(
                        col = 'black',
                        pch = 21,
                        cex = 2,
                        fill = c('dodgerblue','firebrick','orange','darkorchid')
                        ),
                    text = list(
                        lab = c('Deletion','Gain','SNV','WGD')
                        ),
                    padding.text = 5,
                    cex = 1
                    )
                ),
            	title = 'Mutation',
                # Positioning legend on plot
                x = 0.7,
                y = 0.99
            )
        )
	)
