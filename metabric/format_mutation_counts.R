### DESCRIPTION ###################################################################################
# format mutations calls 

setwd('/N/projects/curtis/METABRIC')
### FIND SAMPLES TO USE ###########################################################################
# read in samples to consider
samples <- read.delim(
	'/oak/stanford/groups/ccurtis2/users/khoulaha/BRCA_drivers/sample_list.txt',
	header = FALSE,
 	as.is = TRUE
 	)

# set output dir 
outdir <- "/oak/stanford/groups/ccurtis2/users/khoulaha/BRCA_drivers/mutations"
### CNA COUNTS ####################################################################################
# read in CNA data 
cnas <- read.delim(
	'cbioportal/brca_metabric/data_CNA.txt',
	check.names = FALSE,
	as.is = TRUE
	)
cnas <- cnas[,samples$V1]

# create dataframe of number of cnas
cnas_count <- data.frame(
	sample = colnames(cnas),
	count = colSums(cnas != 0, na.rm = TRUE)
	)
# write to file 
write.table(
	cnas_count,
	file = file.path(outdir, 'Metabric_cna_count.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

# create data frame of number of gains 
gains_count <- data.frame(
	sample = colnames(cnas),
	count = colSums(cnas > 0, na.rm = TRUE)
	)
# write to file 
write.table(
	gains_count,
	file = file.path(outdir, 'Metabric_gains_count.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

# create data frame of number of losses
losses_count <- data.frame(
	sample = colnames(cnas),
	count = colSums(cnas < 0, na.rm = TRUE)
	)
# write to file 
write.table(
	losses_count,
	file = file.path(outdir, 'Metabric_losses_count.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

### CNA SEGS ######################################################################################
# read in segments 
seg <- read.delim(
	'copynumber/CNAall1980_CBS.seg',
	as.is = TRUE
	)
# threshold at |0.1|
seg <- seg[abs(seg$seg.mean) > 0.3,]
# only keep tn pairs 
seg <- seg[seg$SampleID %in% samples$V1,]
# calculate seg length 
seg$length <- seg$loc.end-seg$loc.start
# count segments 
seg_count <- as.data.frame(table(seg$SampleID))
colnames(seg_count) <- c('sample','count')
# write to file 
write.table(
	seg_count,
	file = file.path(outdir, 'Metabric_cna_seg_count.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)
# count number of bps 
cna_bps <- aggregate(seg$length, list(seg$SampleID), sum)
colnames(cna_bps) <- c('sample','count')
# write to file 
write.table(
	cna_bps,
	file = file.path(outdir, 'Metabric_cna_bps_count.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

# count segments gains
seg_gain <- seg[seg$seg.mean > 0,]
seg_gains_count <- as.data.frame(table(seg_gain$SampleID))
colnames(seg_gains_count) <- c('sample','count')
# write to file 
write.table(
	seg_gains_count,
	file = file.path(outdir, 'Metabric_gains_seg_count.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)
# count number of bps in gains
gains_bps <- aggregate(seg_gain$length, list(seg_gain$SampleID), sum)
colnames(gains_bps) <- c('sample','count')
# write to file 
write.table(
	gains_bps,
	file = file.path(outdir, 'Metabric_gains_bps_count.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

# count segments losses
seg_loss <- seg[seg$seg.mean < 0,]
seg_losses_count <- as.data.frame(table(seg_loss$SampleID))
colnames(seg_losses_count) <- c('sample','count')
# write to file 
write.table(
	seg_losses_count,
	file = file.path(outdir, 'Metabric_losses_seg_count.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)
# count number of bps in gains
losses_bps <- aggregate(seg_loss$length, list(seg_loss$SampleID), sum)
colnames(losses_bps) <- c('sample','count')
# write to file 
write.table(
	losses_bps,
	file = file.path(outdir, 'Metabric_losses_bps_count.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)

### SNV ###########################################################################################
# read in snvs 
snv <- read.delim(
	'cbioportal/brca_metabric/data_mutations_extended.txt',
	skip = 1,
	as.is = TRUE
	)
# read in ids profiled
ids <- read.delim(
	'cbioportal/brca_metabric/data_mutations_extended.txt',
	as.is = TRUE,
	header = FALSE
	)
ids <- unlist(strsplit(ids[1,1], split = ' '))[-1]
ids <- ids[ids %in% samples$V1]

# reformat 
snv_count <- as.data.frame(table(snv$Tumor_Sample_Barcode))
# find samples without snvs 
nosnvs <- data.frame(
	sample = ids[which(!ids %in% snv_count$Var1)],
	count = 0
	)
colnames(snv_count) <- c('sample','count')
snv_count <- rbind(snv_count, nosnvs)
snv_count <- snv_count[snv_count$sample %in% samples$V1,]
# write to file 
write.table(
	snv_count,
	file = file.path(outdir, 'Metabric_snv_count.txt'),
	sep = '\t',
	row.names = FALSE,
	quote = FALSE
	)
