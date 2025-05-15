### GET MUTATION RATE #############################################################################
# get_mutation_rate <- function(type, anno) {
# 	if (type == 'snv') {
# 	# read in somatic mutation rate
# 		# mut <- read.delim(
# 		# 	'TCGA_TMB.tsv',
# 		# 	as.is = TRUE
# 		# 	)
# 		mut <- read.delim(
# 			'TCGA_Tumor_Sample_patient_uniq_somatic_mutation_burden.tsv',
# 			as.is = TRUE
# 			)
# 		colnames(mut) <- gsub('nonsynonymous_count','count', colnames(mut))
# 	} else if (type == 'cna') {
# 		# read in cna segments
# 		mut <- read.delim(
# 			"TCGA_total_cna_bp.tsv",
# 			header = TRUE,
# 			as.is = TRUE
# 			)
# 	} else if (type == 'cnaseg') {
# 		# read in cna segments
# 		mut <- read.delim(
# 			"TCGA_segments.tsv",
# 			header = TRUE,
# 			as.is = TRUE
# 			)
# 	} else if (type == 'deletion') {
# 		# read in cna segments, deletions only
# 		mut <- read.delim(
# 			'TCGA_cna_deletion_bp.tsv',
# 			as.is = TRUE
# 			)
# 	} else if (type == 'deletionseg') {
# 		# read in cna segments, deletions only
# 		mut <- read.delim(
# 			'TCGA_indels_deletions_only.tsv',
# 			as.is = TRUE
# 			)
# 	} else if (type == 'indel') {
# 		# read in cna segments
# 		mut <- read.delim(
# 			'TCGA_indels_deletions_bp_only.tsv',
# 			as.is = TRUE
# 			)
# 	} else if (type == 'indelseg') {
# 		# read in cna segments
# 		mut <- read.delim(
# 			'TCGA_indels_deletions_only.tsv',
# 			as.is = TRUE
# 			)
# 	} else if (type == 'amplification') {
# 		# read in cna segments
# 		mut <- read.delim(
# 			'TCGA_cna_amplification_bp.tsv',
# 			as.is = TRUE
# 			)
# 	} else if (type == 'amplificationseg') {
# 		# read in cna segments
# 		mut <- read.delim(
# 			'TCGA_segments_amplifications_only.tsv',
# 			as.is = TRUE
# 			)
# 	} else if (type == 'insertion') {
# 		# read in cna segments
# 		mut <- read.delim(
# 			'TCGA_indels_insertions_bp_only.tsv',
# 			as.is = TRUE
# 			)
# 	} else if (type == 'meth') {
# 		# read in methylation values
# 		mut <- read.csv(
# 			'ciriello_methylation_rates.csv',
# 			as.is = TRUE,
# 			header = TRUE
# 			);
# 		mut <- mut[,c('SampleID','DMI_score')]
# 		mut$SampleID <- gsub('.','-',substr(mut$SampleID, 1, 12), fixed = TRUE)
# 		colnames(mut) <- c('bcr_patient_barcode','count')
# 	} else if (type == 'clones') {
# 		mut <- read.csv('ciriello_number_clones.csv', as.is = TRUE)
# 		mut <- mut[,c('sample_name','number.of.clones')]
# 		mut$sample_name <- substr(mut$sample_name, 1, 12)
# 		colnames(mut) <- c('bcr_patient_barcode','count')
# 	} else {
# 	 	stop("Invalid mutation type specified. Please specifiy either snv, meth or cna ...")
# 	}
# 	# if multiple samples per patient, find median mutation count
# 	mut <- aggregate(mut$count, by = list(mut$bcr_patient_barcode), median, na.rm = TRUE)
# 	colnames(mut) <- c('bcr_patient_barcode','count')
# 	# convert number of mutations to mutation rate
# 	# divide by age at diagnosis
# 	mut_rate <- merge(
# 		mut[,c('count','bcr_patient_barcode')],
# 		anno[,c('bcr_patient_barcode','age_at_initial_pathologic_diagnosis')],
# 		by = 'bcr_patient_barcode'
# 		)
# 	mut_rate$rate <- mut_rate$count/as.numeric(mut_rate$age_at_initial_pathologic_diagnosis)
# 	rownames(mut_rate) <- mut_rate$bcr_patient_barcode
# 	return(mut_rate)
# }

get_mutation_rate <- function(type, anno) {
  # Mapping of mutation types to file paths
  # file_map <- list(
  #   snv = "data/TCGA/TCGA_Tumor_Sample_patient_uniq_somatic_mutation_burden.tsv",
  #   cna = "data/TCGA/mutations/TCGA_total_cna_bp.tsv",
  #   cnaseg = "data/TCGA/mutations/TCGA_segments.tsv",
  #   deletion = here("data", "TCGA", "mutations", "TCGA_cna_deletion_bp.tsv"),
  #   deletionseg = "data/TCGA/mutations/TCGA_indels_deletions_only.tsv",
  #   indel = "data/TCGA/mutations/TCGA_indels_deletions_bp_only.tsv",
  #   indelseg = "data/TCGA/mutations/TCGA_indels_deletions_only.tsv",
  #   amplification = "data/TCGA/mutations/TCGA_cna_amplification_bp.tsv",
  #   amplificationseg = "data/TCGA/mutations/TCGA_segments_amplifications_only.tsv",
  #   insertion = "data/TCGA/mutations/TCGA_indels_insertions_bp_only.tsv"
  # )
  file_map <- list(
    snv = here("data", "TCGA", "mutations", "TCGA_Tumor_Sample_patient_uniq_somatic_mutation_burden.tsv"),
    cna = here("data", "TCGA", "mutations", "TCGA_total_cna_bp.tsv"),
    cnaseg = here("data", "TCGA", "mutations", "TCGA_segments.tsv"),
    deletion = here("data", "TCGA", "mutations", "TCGA_cna_deletion_bp.tsv"),
    deletionseg = here("data", "TCGA", "mutations", "TCGA_segments_deletions_only.tsv"),
    indel_insertion = here("data", "TCGA", "mutations", "TCGA_indels_insertions_bp_only.tsv"),
    indelseg_insertion = here("data", "TCGA", "mutations", "TCGA_indels_insertions_only.tsv"),
    indel_deletion = here("data", "TCGA", "mutations", "TCGA_indels_deletions_bp_only.tsv"),
    indelseg_deletion = here("data", "TCGA", "mutations", "TCGA_indels_deletions_only.tsv"),
    amplification = here("data", "TCGA", "mutations", "TCGA_cna_amplification_bp.tsv"),
    amplificationseg = here("data", "TCGA", "mutations", "TCGA_segments_amplifications_only.tsv")
  )

  cat("Reading mutation data for type:", type, "\n")
  if (type %in% names(file_map)) {
    fn <- file_map[[type]]
    cat("Reading file:", fn, "\n")
    mut <- read.delim(fn, as.is = TRUE, header = TRUE)
  } else if (type == "meth") {
    mut <- read.csv("ciriello_methylation_rates.csv", as.is = TRUE)
    mut <- mut[, c("SampleID", "DMI_score")]
    mut$SampleID <- gsub(".", "-", substr(mut$SampleID, 1, 12), fixed = TRUE)
    colnames(mut) <- c("bcr_patient_barcode", "count")
  } else if (type == "clones") {
    mut <- read.csv("ciriello_number_clones.csv", as.is = TRUE)
    mut <- mut[, c("sample_name", "number.of.clones")]
    mut$sample_name <- substr(mut$sample_name, 1, 12)
    colnames(mut) <- c("bcr_patient_barcode", "count")
  } else {
    stop("Invalid mutation type specified. Please specify one of: snv, cna, meth, clones, etc.")
  }

  # Special renaming for SNV data
  if (type == "snv") {
    colnames(mut) <- gsub("nonsynonymous_count", "count", colnames(mut))
  }

  # Aggregate by patient barcode using median
  mut <- aggregate(mut$count, by = list(mut$bcr_patient_barcode), FUN = median, na.rm = TRUE)
  colnames(mut) <- c("bcr_patient_barcode", "count")

  mut_rate <- merge(
    mut[, c("count", "bcr_patient_barcode")],
    anno[, c("bcr_patient_barcode", "age_at_initial_pathologic_diagnosis")],
    by = "bcr_patient_barcode"
  )
  mut_rate$rate <- mut_rate$count / as.numeric(mut_rate$age_at_initial_pathologic_diagnosis)
  rownames(mut_rate) <- mut_rate$bcr_patient_barcode
  return(mut_rate)
  #
  # # Merge with clinical annotation and compute mutation rate
  # mut_rate <- merge(
  #   mut,
  #   anno[, c("bcr_patient_barcode", "age_at_initial_pathologic_diagnosis")],
  #   by = "bcr_patient_barcode"
  # )
  # mut_rate$rate <- mut_rate$count / as.numeric(mut_rate$age_at_initial_pathologic_diagnosis)
  # rownames(mut_rate) <- mut_rate$bcr_patient_barcode
  #
  # return(mut_rate)
}

### CALCULATE MUTATION RATE RATIO #################################################################
calculate_mutation_rate_ratio <- function(int, mut_rate, ddr, wt, anno, cancer) {
  # sample ddr and non ddr samples
  #set.seed(42)  # You can choose any integer seed value
  ddr_sample <- sample(ddr, length(ddr), replace = TRUE)
  if (cancer == 'BRCA') {
    wt_sample <- standardize_clinical_characteristics_breast(
      anno = anno,
      wt = wt,
      ddr = ddr
    )
  } else if (cancer == 'OV') {
    wt_sample <- standardize_clinical_characteristics_ovarian(
      anno = anno,
      wt = wt,
      ddr = ddr
    )
  } else {
    stop("Please specify a valid cancer type. Options are BRCA or OV ...")
  }

 # wt_sample <- sample(wt, length(ddr), replace = TRUE)

  # calculate median mutation rate of ddr and non ddr samples
  ddr_median <- median(mut_rate[ddr_sample, "rate"])
  wt_median <- median(mut_rate[wt_sample, "rate"])

  # calculate ratio
  ddr_ratio <- ddr_median / wt_median

  # calculate ratios
  data.frame(
    int = int,
    num = length(ddr),
    wt_median = wt_median,
    ddr_median = ddr_median,
    ratio = ddr_ratio,
    incidence_two = ddr_ratio * ddr_ratio,
    incidence_three = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)),
    incidence_four = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)),
    incidence_five = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)),
    incidence_six = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)) * (ddr_ratio^(1 / 5)),
    incidence_seven = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)) * (ddr_ratio^(1 / 5)) * (ddr_ratio^(1 / 6)),
    incidence_eight = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)) * (ddr_ratio^(1 / 5)) * (ddr_ratio^(1 / 6)) * (ddr_ratio^(1 / 7)),
    incidence_nine = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)) * (ddr_ratio^(1 / 5)) * (ddr_ratio^(1 / 6)) * (ddr_ratio^(1 / 7)) * (ddr_ratio^(1 / 8)),
    incidence_ten = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)) * (ddr_ratio^(1 / 5)) * (ddr_ratio^(1 / 6)) * (ddr_ratio^(1 / 7)) * (ddr_ratio^(1 / 8)) * (ddr_ratio^(1 / 9)),
    incidence_eleven = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)) * (ddr_ratio^(1 / 5)) * (ddr_ratio^(1 / 6)) * (ddr_ratio^(1 / 7)) * (ddr_ratio^(1 / 8)) * (ddr_ratio^(1 / 9)) * (ddr_ratio^(1 / 10)),
    incidence_twelve = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)) * (ddr_ratio^(1 / 5)) * (ddr_ratio^(1 / 6)) * (ddr_ratio^(1 / 7)) * (ddr_ratio^(1 / 8)) * (ddr_ratio^(1 / 9)) * (ddr_ratio^(1 / 10)) * (ddr_ratio^(1 / 11)),
    incidence_thirteen = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)) * (ddr_ratio^(1 / 5)) * (ddr_ratio^(1 / 6)) * (ddr_ratio^(1 / 7)) * (ddr_ratio^(1 / 8)) * (ddr_ratio^(1 / 9)) * (ddr_ratio^(1 / 10)) * (ddr_ratio^(1 / 11)) * (ddr_ratio^(1 / 12)),
    incidence_fourteen = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)) * (ddr_ratio^(1 / 5)) * (ddr_ratio^(1 / 6)) * (ddr_ratio^(1 / 7)) * (ddr_ratio^(1 / 8)) * (ddr_ratio^(1 / 9)) * (ddr_ratio^(1 / 10)) * (ddr_ratio^(1 / 11)) * (ddr_ratio^(1 / 12)) * (ddr_ratio^(1 / 13)),
    incidence_fifteen = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)) * (ddr_ratio^(1 / 5)) * (ddr_ratio^(1 / 6)) * (ddr_ratio^(1 / 7)) * (ddr_ratio^(1 / 8)) * (ddr_ratio^(1 / 9)) * (ddr_ratio^(1 / 10)) * (ddr_ratio^(1 / 11)) * (ddr_ratio^(1 / 12)) * (ddr_ratio^(1 / 13)) * (ddr_ratio^(1 / 14)),
    incidence_sixteen = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)) * (ddr_ratio^(1 / 5)) * (ddr_ratio^(1 / 6)) * (ddr_ratio^(1 / 7)) * (ddr_ratio^(1 / 8)) * (ddr_ratio^(1 / 9)) * (ddr_ratio^(1 / 10)) * (ddr_ratio^(1 / 11)) * (ddr_ratio^(1 / 12)) * (ddr_ratio^(1 / 13)) * (ddr_ratio^(1 / 14)) * (ddr_ratio^(1 / 15))
  )
}


# calculate_mutation_rate_ratio <- function(int,
#                                           date,
#                                           mut_rate,
#                                           ddr,
#                                           wt,
#                                           anno,
#                                           cancer,
#                                           gene,
#                                           mutation,
#                                           adj_flag = TRUE) {
#   # sample ddr and non ddr samples
#   ddr_sample <- sample(ddr, length(ddr), replace = TRUE)
#   # wt_sample <- sample(wt, length(ddr), replace = TRUE)
#   fn <- paste(date, cancer, gene, mutation,
#     "mutation_rate_ratio.xlsx",
#     sep = "_"
#   )
#   if (adj_flag) {
#     adj_dir <- "adjusted"
#   } else {
#     adj_dir <- "unadjusted"
#   }
#   fn <- here("output", "data", "TCGA", adj_dir, fn)
#   if (cancer == "BRCA") {
#     res <- standardize_clinical_characteristics_breast(
#       anno = anno,
#       wt = wt,
#       ddr = ddr
#       # out_xlsx = fn
#     )
#     wt_sample <- res
#     # if (adj_flag) {
#     #   print("Using adjusted weights")
#     #   wt_sample <- res$adj
#     # } else {
#     #   print("Using unadjusted weights")
#     #   wt_sample <- res$unadj
#     # }
#   } else if (cancer == "OV") {
#     res <- standardize_clinical_characteristics_ovarian(
#       anno = anno,
#       wt = wt,
#       ddr = ddr
#       # out_xlsx = fn
#     )
#     wt_sample <- res
#     # if (adj_flag) {
#     #   print("Using adjusted weights")
#     #   wt_sample <- res$adj
#     # } else {
#     #   print("Using unadjusted weights")
#     #   wt_sample <- res$unadj
#     # }
#   } else {
#     stop("Please specify a valid cancer type. Options are BRCA or OV ...")
#   }
#
#   # calculate median mutation rate of ddr and non ddr samples
#   ddr_median <- median(mut_rate[ddr_sample, "rate"])
#   wt_median <- median(mut_rate[wt_sample, "rate"], na.rm = TRUE)
#
#   # calculate ratio
#   ddr_ratio <- ddr_median / wt_median
#
#   # calculate ratios
#   df <- data.frame(
#     int = int,
#     num = length(ddr),
#     wt_median = wt_median,
#     ddr_median = ddr_median,
#     ratio = ddr_ratio,
#     incidence_two = ddr_ratio * ddr_ratio,
#     incidence_three = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)),
#     incidence_four = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)),
#     incidence_five = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)),
#     incidence_six = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)) * (ddr_ratio^(1 / 5)),
#     incidence_seven = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)) * (ddr_ratio^(1 / 5)) * (ddr_ratio^(1 / 6)),
#     incidence_eight = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)) * (ddr_ratio^(1 / 5)) * (ddr_ratio^(1 / 6)) * (ddr_ratio^(1 / 7)),
#     incidence_nine = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)) * (ddr_ratio^(1 / 5)) * (ddr_ratio^(1 / 6)) * (ddr_ratio^(1 / 7)) * (ddr_ratio^(1 / 8)),
#     incidence_ten = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)) * (ddr_ratio^(1 / 5)) * (ddr_ratio^(1 / 6)) * (ddr_ratio^(1 / 7)) * (ddr_ratio^(1 / 8)) * (ddr_ratio^(1 / 9)),
#     incidence_eleven = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)) * (ddr_ratio^(1 / 5)) * (ddr_ratio^(1 / 6)) * (ddr_ratio^(1 / 7)) * (ddr_ratio^(1 / 8)) * (ddr_ratio^(1 / 9)) * (ddr_ratio^(1 / 10)),
#     incidence_twelve = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)) * (ddr_ratio^(1 / 5)) * (ddr_ratio^(1 / 6)) * (ddr_ratio^(1 / 7)) * (ddr_ratio^(1 / 8)) * (ddr_ratio^(1 / 9)) * (ddr_ratio^(1 / 10)) * (ddr_ratio^(1 / 11)),
#     incidence_thirteen = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)) * (ddr_ratio^(1 / 5)) * (ddr_ratio^(1 / 6)) * (ddr_ratio^(1 / 7)) * (ddr_ratio^(1 / 8)) * (ddr_ratio^(1 / 9)) * (ddr_ratio^(1 / 10)) * (ddr_ratio^(1 / 11)) * (ddr_ratio^(1 / 12)),
#     incidence_fourteen = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)) * (ddr_ratio^(1 / 5)) * (ddr_ratio^(1 / 6)) * (ddr_ratio^(1 / 7)) * (ddr_ratio^(1 / 8)) * (ddr_ratio^(1 / 9)) * (ddr_ratio^(1 / 10)) * (ddr_ratio^(1 / 11)) * (ddr_ratio^(1 / 12)) * (ddr_ratio^(1 / 13)),
#     incidence_fifteen = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)) * (ddr_ratio^(1 / 5)) * (ddr_ratio^(1 / 6)) * (ddr_ratio^(1 / 7)) * (ddr_ratio^(1 / 8)) * (ddr_ratio^(1 / 9)) * (ddr_ratio^(1 / 10)) * (ddr_ratio^(1 / 11)) * (ddr_ratio^(1 / 12)) * (ddr_ratio^(1 / 13)) * (ddr_ratio^(1 / 14)),
#     incidence_sixteen = ddr_ratio * ddr_ratio * (ddr_ratio^(1 / 2)) * (ddr_ratio^(1 / 3)) * (ddr_ratio^(1 / 4)) * (ddr_ratio^(1 / 5)) * (ddr_ratio^(1 / 6)) * (ddr_ratio^(1 / 7)) * (ddr_ratio^(1 / 8)) * (ddr_ratio^(1 / 9)) * (ddr_ratio^(1 / 10)) * (ddr_ratio^(1 / 11)) * (ddr_ratio^(1 / 12)) * (ddr_ratio^(1 / 13)) * (ddr_ratio^(1 / 14)) * (ddr_ratio^(1 / 15))
#   )
#   return(df)
# }


### CALCULATE MEDIAN ESTIMATED INCIDENCES #########################################################
calculate_median_est_incidence <- function() {
  # create grid of cancer, gene and mutation type
  par <- expand.grid(
    cancer = c("BRCA", "OV"), gene = c("BRCA1", "BRCA2"),
    mutation = c("snv", "cna", "deletion", "amplification", "indel")
  ) # ,'insertion', 'indelseg','amplificationseg','cnaseg','deletionseg'))
  # loop over grid and calculate medians
  incidence <- apply(
    par,
    1,
    function(x) {
      filename <- paste(date, x["cancer"], x["gene"], x["mutation"], "incidence_estimates.tsv", sep = "_")
      print(filename)
      # filename <- list.files(pattern = fileflag)
      tmp <- read.delim(filename, as.is = TRUE)
      # read in observed
      # # read in observed values
      observed <- read.delim("observed_incidence_rates.tsv", as.is = TRUE)
      # extract median and CIs
      ob_median <- observed[observed$Cancer == x["cancer"] & observed$Gene == x["gene"], "Median"]
      ob_L95 <- observed[observed$Cancer == x["cancer"] & observed$Gene == x["gene"], "L95"]
      ob_U95 <- observed[observed$Cancer == x["cancer"] & observed$Gene == x["gene"], "U95"]
      ob_n <- observed[observed$Cancer == x["cancer"] & observed$Gene == x["gene"], "Num"]
      # calculate standard deviation
      ob_sd <- sqrt(ob_n) * (ob_U95 - ob_L95) / 3.92
      ob_dist <- rnorm(10000, mean = ob_median, sd = ob_sd)
      # plots <- create_incidence_segplot(tmp, ob_median = ob_median, ob_L95 = ob_L95, ob_U95 = ob_U95, filename = paste(date, x['cancer'], x['gene'], x['mutation'], 'segplot.tiff', sep = '_'),
      # 		main = paste(x['cancer'], x['gene']), driver_max = 15, ylimit = 100, yat = seq(0,100,25))
      # calculate median
      # apply(tmp[,grep('incidence', colnames(tmp))], 2, median)
      # apply KS test
      ks <- run_ks_test(tmp, ob_dist)
      ks$group <- paste(x["cancer"], x["gene"], x["mutation"], sep = "|")
      return(ks)
    }
  )
  colnames(incidence) <- apply(par, 1, paste, collapse = "_")
  rownames(incidence) <- c("two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve", "thirteen")
  return(incidence)
}

calculate_CIs <- function(x) {
  x <- x[order(x)]
  return(x[c(250, 9750)])
}

run_ks_test <- function(df, ob_dist) {
  # run Kolmogorov-Smirnov test for each number of drivers
  stats <- do.call(rbind, apply(
    df[, grep("ratio|incidence", colnames(df))],
    2,
    function(x) {
      stats_tmp <- ks.test(
        x,
        ob_dist
      )
      data.frame(
        D = stats_tmp$statistic,
        P = stats_tmp$p.value
      )
    }
  ))
  stats$driver <- gsub("incidence_", "", grep("ratio|incidence", colnames(df), value = TRUE))
  return(stats)
}

create_incidence_segplot <- function(tmp, mut, ob_median, ob_L95, ob_U95, filename, main,
                                     driver_max = NULL, yat = NULL, ylimits = NULL) {
  # create CI for bootstrap
  mean_inc <- apply(tmp[, grep("ratio|incidence", colnames(tmp))], 2, median)
  # mean_inc <- apply(tmp[, grep("ratio|incidence", colnames(tmp))], 2, median)
  sd_inc <- apply(tmp[, grep("ratio|incidence", colnames(tmp))], 2, sd)
  CI_inc <- apply(tmp[, grep("ratio|incidence", colnames(tmp))], 2, calculate_CIs)


  plot_data <- data.frame(
    num = 1:length(mean_inc), mean = mean_inc,
    L95 = CI_inc[1, ], U95 = CI_inc[2, ]
  )
  # set maximum xlimits
  if (!is.null(driver_max)) {
    plot_data <- plot_data[plot_data$num <= driver_max, ]
  }
  plot_data$num <- factor(plot_data$num)
  max_y <- max(plot_data$U95)
  # set ylimits and yat
  ylimits <- ifelse(
    is.null(ylimits),
    max(max_y + 0.2 * max_y, ob_U95 + ob_U95 * 0.1),
    ylimits
  )
  if (is.null(yat)) {
    yat <- seq(0, ylimits, 10)
    # 	yat <- round(seq(0, ylimit, length.out = 4), digits = -1)
  } else {
    yat <- yat
  }

  # set xlab
  # mut <- gsub("seg", "", strsplit(filename, "_")[[1]][4])
  labels <- c(
    "Number SNV drivers",
    "Number CNA drivers",
    "Number CNA drivers",
    "Number CNA deletion drivers",
    "Number CNA gain drivers",
    "Number small deletion drivers"
  )
  names(labels) <- c(
    "snv",
    "cnaseg",
    "cna",
    "deletion",
    "amplification",
    "indel"
  )
  # create segplot
  create.scatterplot(
    mean ~ num,
    data = plot_data,
    filename = filename,
    main = main,
    ylab.label = "Incidence",
    xlab.label = labels[mut],
    xaxis.lab = 1:length(mean_inc),
    yat = yat,
    ylimits = c(0, ylimits),
    y.error.up = plot_data$U95 - plot_data$mean,
    y.error.down = plot_data$mean - plot_data$L95,
    resolution = 300,
    add.rectangle = TRUE,
    alpha.rectangle = 0.5,
    col.rectangle = "grey",
    xleft.rectangle = 0,
    error.bar.lwd = 0.5,
    xright.rectangle = length(mean_inc) + 1.5,
    ybottom.rectangle = ob_L95,
    ytop.rectangle = ob_U95,
    width = 9,
    # plot.horizontal = FALSE,
    abline.h = ob_median
  )
}

# create_incidence_segplot <- function(tmp, ob_median, ob_L95, ob_U95, filename, main, driver_max = NULL) {
#   # create CI for bootstrap
#   mean_inc <- apply(tmp[,grep('incidence', colnames(tmp))], 2, mean)
#   sd_inc <-  apply(tmp[,grep('incidence', colnames(tmp))], 2, sd)
#   CI_inc <- 1.96*(sd_inc/sqrt(nrow(tmp)))
#
#   plot_data <- data.frame(num = 2:(length(mean_inc)+1), mean = mean_inc, CI = CI_inc)
#   plot_data$L95 <- plot_data$mean-plot_data$CI
#   plot_data$U95 <- plot_data$mean+plot_data$CI
#   # set maximum xlimits
#   if (!is.null(driver_max)) {
#     plot_data <- plot_data[plot_data$num <= driver_max,]
#   }
#   plot_data$num <- factor(plot_data$num)
#   max_y <- max(plot_data$mean+plot_data$CI)
#   ylimits <- max(max_y+0.2*max_y, ob_U95+ob_U95*0.1)
#   # create segplot
#   create.scatterplot(
#     mean ~ num,
#     data = plot_data,
#     filename = filename,
#     main = main,
#     ylab.label = 'Incidence',
#     xlab.label = 'Number Drivers',
#     xaxis.lab = 2:(length(mean_inc)+1),
#     ylimits = c(0, ylimits),
#     yat = seq(0, ylimits, 10),
#     #centers = plot_data$mean,
#     y.error.up = plot_data$CI,
#     resolution = 300,
#     add.rectangle = TRUE,
#     alpha.rectangle = 0.5,
#     col.rectangle = 'grey',
#     xleft.rectangle = 0,
#     error.bar.lwd = 0.5,
#     xright.rectangle = length(mean_inc)+1.5,
#     ybottom.rectangle = ob_L95,
#     ytop.rectangle = ob_U95,
#     width = 9,
#     #plot.horizontal = FALSE,
#     abline.h = ob_median
#   )
# }




#### CREATE BARPLOT OF MEDIANS ####################################################################
create_barplot_of_medians <- function(incidence) {
  # create barplot of medians
  incidence <- as.data.frame(incidence)
  incidence$driver_number <- rownames(incidence)
  plot_data <- gather(incidence, "group", "incidence", -driver_number)

  observed <- read.delim("observed_incidence_rates.tsv", as.is = TRUE)
  observed$group <- paste(observed$Cancer, observed$Gene, sep = "_")
  # reformat observed data
  ob_df <- data.frame(
    driver_number = "observed", group = c(paste(observed$group, "snv", sep = "_"), paste(observed$group, "cna", sep = "_")),
    incidence = rep(observed$Median, 2)
  )
  plot_data <- rbind(plot_data, ob_df)
  # make sure observed is last
  plot_data$index <- NA
  plot_data[plot_data$driver_number == "two", "index"] <- "a"
  plot_data[plot_data$driver_number == "three", "index"] <- "b"
  plot_data[plot_data$driver_number == "four", "index"] <- "c"
  plot_data[plot_data$driver_number == "five", "index"] <- "d"
  plot_data[plot_data$driver_number == "observed", "index"] <- "e"

  create.barplot(
    incidence ~ group,
    data = plot_data,
    filename = paste0(date, "_incidence_means_barplot.tiff"),
    groups = plot_data$index,
    resolution = 300,
    col = c(default.colours(4), "black"),
    # xaxis.rot = 90,
    ylab.label = "Incidence",
    xaxis.lab = c(
      "BRCA\nBRCA1\nCNA", "BRCA\nBRCA1\nSNV", "BRCA\nBRCA2\nCNA", "BRCA\nBRCA2\nSNV",
      "OV\nBRCA1\nCNA", "OV\nBRCA1\nSNV", "OV\nBRCA2\nCNA", "OV\nBRCA1\nSNV"
    ),
    xaxis.cex = 1,
    xlab.label = "",
    width = 8,
    legend = list(
      inside = list(
        fun = draw.key,
        args = list(
          key = list(
            points = list(
              col = c(default.colours(4), "black"),
              pch = 21,
              cex = 1.5,
              fill = c(default.colours(4), "black")
            ),
            text = list(
              lab = c("Two", "Three", "Four", "Five", "Observed")
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
    CI <- med + 1.96 * (stand_dev / sqrt(length(x)))
  } else {
    CI <- med - 1.96 * (stand_dev / sqrt(length(x)))
  }
  return(CI)
}

### TEST MUTATION RATE ############################################################################
test_mutation_rate <- function(df) {
  results <- list()
  for (i in unique(df$group)) {
    tmp <- df[df$group == i, ]
    brca1 <- get.utest.p.and.medianfoldchange(
      tmp[tmp$variant %in% c("WT", "BRCA1"), "rate"],
      tmp[tmp$variant %in% c("WT", "BRCA1"), "variant"] == "WT",
      tmp[tmp$variant %in% c("WT", "BRCA1"), "variant"] == "BRCA1",
      logged = FALSE
    )
    brca2 <- get.utest.p.and.medianfoldchange(
      tmp[tmp$variant %in% c("WT", "BRCA2"), "rate"],
      tmp[tmp$variant %in% c("WT", "BRCA2"), "variant"] == "WT",
      tmp[tmp$variant %in% c("WT", "BRCA2"), "variant"] == "BRCA2",
      logged = FALSE
    )
    results[[i]] <- data.frame(
      group = i,
      variant = c("BRCA1", "BRCA2"),
      FC = c(brca1[2], brca2[2]),
      P = c(brca1[1], brca2[1])
    )
  }
  results <- do.call(rbind, results)
  results$FDR <- p.adjust(results$P, method = "fdr")
  return(results)
}

### GENERATE INCIDENCE PLOTS ##########################################################
get_plot_limits <- function(cancer, mutation, gene) {
  # return <- list(driver_max=15,
  #                ylimits =100,
  #                yat = c(0, 25, 50, 75, 100))
  # # 
  adj_flag <- ifelse(mutation %in% c("cnaseg", "deletionseg"), TRUE, FALSE)

  if (adj_flag) {
    limits_table <- tibble::tibble(
      cancer = c(
        "BRCA", "BRCA", "OV", "OV",
        "BRCA", "BRCA", "OV", "OV"
      ),
      mutation = c(
        "cnaseg", "cnaseg", "cnaseg", "cnaseg",
        "deletionseg", "deletionseg", "deletionseg", "deletionseg"
      ),
      gene = c(
        "BRCA1", "BRCA2",
        "BRCA1", "BRCA2",
        "BRCA1", "BRCA2",
        "BRCA1", "BRCA2"
      ),
      driver_max = c(
        5, 15,
        15, 15,
        15, 15,
        15, 15
      ),
      ylimits = c(
        100, 100, 75, 100,
        100, 100, 100, 100
      ),
      yat = list(
        c(0, 25, 50, 75, 100),
        c(0, 25, 50, 75, 100),
        c(0, 25, 50, 75, 100),
        c(0, 25, 50, 75, 100),
        c(0, 25, 50, 75, 100),
        c(0, 25, 50, 75, 100),
        c(0, 25, 50, 75, 100),
        c(0, 25, 50, 75, 100)
      )
    )
  } else {
    limits_table <- tibble::tibble(
      cancer = c(
        "BRCA", "BRCA", "OV", "OV",
        "BRCA", "BRCA", "OV", "OV"
      ),
      mutation = c(
        "cna", "cna", "cna", "cna",
        "deletion", "deletion", "deletion", "deletion"
      ),
      gene = c(
        "BRCA1", "BRCA2",
        "BRCA1", "BRCA2",
        "BRCA1", "BRCA2",
        "BRCA1", "BRCA2"
      ),
      driver_max = c(
        15, 15, # cna
        15, 15,
        15, 5, # deletion
        5, 15
      ),
      ylimits = c(
        80, 20, 75, 25,
        200, 250, 600, 200
      ),
      yat = list(
        c(0, 20, 40, 60, 80), # cna start
        c(0, 5, 10, 15, 20),
        c(0, 25, 50, 75),
        c(0, 5, 10, 15, 20, 25), # cna end
        c(0, 50, 100, 150, 200), # deletion start
        c(0, 50, 100, 150, 200, 250),
        c(0, 200, 400, 600),
        c(0, 50, 100, 150, 200) # deletion end
      )
    )
  }


  # Match only the first valid row
  idx <- which(limits_table$cancer == cancer &
    limits_table$mutation == mutation &
    limits_table$gene == gene)

  if (length(idx) == 1) {
    driver_max <- limits_table$driver_max[idx]
    ylimits <- limits_table$ylimits[idx]
    yat <- limits_table$yat[[idx]]

    cat("Using", driver_max, ylimits, paste(yat, collapse = " "), "\n")

    return(list(
      driver_max = driver_max,
      ylimits = ylimits,
      yat = yat
    ))
  } else {
    return(list(driver_max = NULL, ylimits = NULL, yat = NULL))
  }
}



calculate_median_est_incidence_detail <- function(date,
                                                  cancer,
                                                  gene,
                                                  mutation) {
  # adj_dir <- if (adj_flag) "adjusted" else "unadjusted"
  # message("Using ", adj_dir, " weights")

  # Define parameter combinations
  params <- expand.grid(cancer = cancer, gene = gene, mutation = mutation)
  print(params)

  # Read observed data once
  observed_data <- readxl::read_xlsx(here("data", "TCGA", "SupplementaryTable01_BRCA_OVCA_SIR.xlsx"), "Formated")

  process_combination <- function(x) {
    cancer <- x["cancer"]
    gene <- x["gene"]
    mutation <- x["mutation"]

    # Construct input/output paths
    infile <- here(
      "output/data", "TCGA/European",
      paste(date, cancer, gene, mutation, "incidence_estimates.tsv",
        sep = "_"
      )
    )
    message("Reading file: ", infile)


    # Read incidence estimate
    est_data <- read.delim(infile, as.is = TRUE)

    # Extract observed statistics
    obs_row <- subset(observed_data, Cancer == cancer & Gene == gene)
    ob_median <- obs_row$Median[[1]]
    ob_L95 <- obs_row$L95[[1]]
    ob_U95 <- obs_row$U95[[1]]
    ob_n <- obs_row$Num[[1]]
    message("n: ", ob_n)

    ob_sd <- sqrt(ob_n) * (ob_U95 - ob_L95) / 3.92
    #set.seed(42)
    ob_dist <- rnorm(10000, mean = ob_median, sd = ob_sd)


    # Plot configuration
    limits <- get_plot_limits(cancer, mutation, gene)
    driver_max <- limits$driver_max
    ylimits <- limits$ylimits
    yat <- limits$yat
    plotfile <- here(
      "output", "figures", "TCGA/European",
      paste(date, cancer, gene, mutation,
        "segplot.tiff",
        sep = "_"
      )
    )
    print(paste("Saving fig.:", plotfile))

    create_incidence_segplot(est_data, mutation,
      ob_median, ob_L95, ob_U95,
      filename = plotfile,
      main = paste(cancer, gene),
      driver_max = driver_max,
      ylimits = ylimits, yat = yat
    )

    ks <- run_ks_test(est_data, ob_dist)
    ks$group <- paste(cancer, gene, mutation, sep = "|")
    return(ks)
  }

  result <- apply(params, 1, process_combination)
  return(result)
}

# calculate_median_est_incidence_detail <- function(date,
#                                                   cancer,
#                                                   gene,
#                                                   mutation,
#                                                   adj_flag = TRUE) {
#   if (adj_flag) {
#     adj_dir <- "adjusted"
#   } else {
#     adj_dir <- "unadjusted"
#   }
#   print(paste("Using", adj_dir, "weights"))
#   # create grid of cancer, gene and mutation type
#   par <- expand.grid(
#     cancer = c(cancer), gene = c(gene),
#     mutation = c(mutation)
#     # mutation = c("snv", "cna", "deletion", "amplification", "indel")
#   ) # ,'insertion', 'indelseg','amplificationseg','cnaseg','deletionseg'))
#   # loop over grid and calculate medians
#   incidence <- apply(
#     par,
#     1,
#     function(x) {
#       filename <- paste(date, x["cancer"], x["gene"], x["mutation"], "incidence_estimates.tsv", sep = "_")
#       filename <- here("data", "TCGA", adj_dir, filename)
#       # filename <- here("output", "data", "TCGA", adj_dir, filename)
#
#       # filename <- paste(date, x["cancer"], x["gene"], x["mutation"], "incidence_estimates.tsv", sep = "_")
#       print(filename)
#       # filename <- list.files(pattern = fileflag)
#       tmp <- read.delim(filename, as.is = TRUE)
#       # read in observed
#       # # read in observed value
#       fn <- here("data", "TCGA", "SupplementaryTable01_BRCA_OVCA_SIR.xlsx")
#       observed <- readxl::read_xlsx(fn, "Formated")
#       # read.delim("SupplementaryTable01_BRCA_OVCA_SIR.xlsx", as.is = TRUE)
#       #      observed <- read.delim("observed_incidence_rates.tsv", as.is = TRUE)
#       # extract median and CIs
#       ob_median <- observed[observed$Cancer == x["cancer"] & observed$Gene == x["gene"], "Median"][[1]]
#       ob_L95 <- observed[observed$Cancer == x["cancer"] & observed$Gene == x["gene"], "L95"][[1]]
#       ob_U95 <- observed[observed$Cancer == x["cancer"] & observed$Gene == x["gene"], "U95"][[1]]
#       ob_n <- observed[observed$Cancer == x["cancer"] & observed$Gene == x["gene"], "Num"][[1]]
#       print(paste0("n:", ob_n))
#       # calculate standard deviation
#       ob_sd <- sqrt(ob_n) * (ob_U95 - ob_L95) / 3.92
#       ob_dist <- rnorm(10000, mean = ob_median, sd = ob_sd)
#       fn <- paste(date, x["cancer"], x["gene"], x["mutation"], "segplot.tiff", sep = "_")
#       fn <- here("output", "figures", "TCGA", adj_dir, fn)
#       print(paste("Saving fig.:", fn))
#
#       if (mutation == "cna" && !adj_flag & gene == "BRCA1") {
#         plots <- create_incidence_segplot(tmp,
#           ob_median = ob_median, ob_L95 = ob_L95, ob_U95 = ob_U95,
#           filename = fn,
#           main = paste(x["cancer"], x["gene"]),
#           driver_max = 15,
#           ylimits = 80, yat = seq(0, 80, 20)
#         )
#       } else if (mutation == "cna" && !adj_flag & gene == "BRCA2") {
#         plots <- create_incidence_segplot(tmp,
#           ob_median = ob_median, ob_L95 = ob_L95, ob_U95 = ob_U95,
#           filename = fn,
#           main = paste(x["cancer"], x["gene"]),
#           driver_max = 15,
#           ylimits = 20, yat = seq(0, 20, 5)
#         )
#       } else {
#         plots <- create_incidence_segplot(tmp,
#           ob_median = ob_median, ob_L95 = ob_L95, ob_U95 = ob_U95,
#           filename = fn,
#           main = paste(x["cancer"], x["gene"])
#         )
#       }
#       # calculate median
#       # apply(tmp[,grep('incidence', colnames(tmp))], 2, median)
#       # apply KS test
#       ks <- run_ks_test(tmp, ob_dist)
#       ks$group <- paste(x["cancer"], x["gene"], x["mutation"], sep = "|")
#       return(ks)
#     }
#   )
#   # colnames(incidence) <- apply(par, 1, paste, collapse = "_")
#   # rownames(incidence) <- c("two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve", "thirteen")
#   return(incidence)
# }

### STANDARDIZE WT ON CLINICAL VARIABLES ##########################################################
standardize_clinical_characteristics_breast <- function(anno, wt, ddr) {
  # read in subtype data
  # subtype <- read.csv(
  # 	'~/GermlineSomaticAssociations/genome-wide/output/somatic_gwas/pancan/driver_rate_comparison/TCGA_gyn_supplementary_table4.csv',
  # 	as.is = TRUE,
  # 	skip = 1
  # 	);
  # # Read subtype & build carrier summary ----
  subtype <- read_xlsx("data/TCGA/mmc4.xlsx", skip = 1) |>
    select(Sample.ID, BRCA_Subtype_PAM50)

  subtype <- subtype[, c("Sample.ID", "BRCA_Subtype_PAM50")]
  # find carrier subtypes and calculate prop
  car_sub <- subtype[subtype$Sample.ID %in% ddr, ]
  car_sub_prop <- table(car_sub$BRCA_Subtype_PAM50) / nrow(car_sub)
  # only keep wt
  wt_sub <- subtype[!subtype$Sample.ID %in% ddr, ]
  # only keep wt that have mutation data
  wt_sub <- wt_sub[wt_sub$Sample.ID %in% wt, ]

  # find carrier pathologic stage
  car_stage <- anno[anno$bcr_patient_barcode %in% ddr, ]
  car_stage_prop <- table(car_stage$pathologic_stage) / nrow(car_stage)
  # only keep wt
  wt_stage <- anno[!anno$bcr_patient_barcode %in% ddr, c("bcr_patient_barcode", "pathologic_stage")]
  # only keep wt that have mutation data
  wt_stage <- wt_stage[wt_stage$bcr_patient_barcode %in% wt, ]
  # merge stage and subtypes
  colnames(wt_stage) <- gsub("bcr_patient_barcode", "Sample.ID", colnames(wt_stage))
  wt_df <- merge(wt_sub, wt_stage, by = "Sample.ID")
  wt_df$group <- paste(wt_df$BRCA_Subtype_PAM50, wt_df$pathologic_stage, sep = "|")
  # get weights of each group
  group_weights <- as.data.frame(table(wt_df$group))
  colnames(group_weights) <- c("group", "prop")
  group_weights$prop <- group_weights$prop / nrow(wt_df)
  # add weights
  wt_df <- merge(wt_df, group_weights, by = "group")
  # randomly sample wt
  wt_sample <- sample_n(wt_df, size = length(ddr), weight = wt_df$prop)
  return(wt_sample$Sample.ID)
}

standardize_clinical_characteristics_ovarian <- function(anno, wt, ddr) {
  # find carrier pathologic stage
  car_stage <- anno[anno$bcr_patient_barcode %in% ddr, ]
  car_stage_prop <- table(car_stage$pathologic_stage) / nrow(car_stage)
  # only keep wt
  wt_stage <- anno[!anno$bcr_patient_barcode %in% ddr, c("bcr_patient_barcode", "clinical_stage")]
  # only keep wt that have mutation data
  wt_stage <- wt_stage[wt_stage$bcr_patient_barcode %in% wt, ]
  # calculate stage weight
  stage_weights <- as.data.frame(table(wt_stage$clinical_stage))
  colnames(stage_weights) <- c("clinical_stage", "prop")
  stage_weights$prop <- stage_weights$prop / nrow(wt_stage)
  # add weights
  wt_df <- merge(wt_stage, stage_weights, by = "clinical_stage")
  # randomly sample wt
  wt_sample <- sample_n(wt_df, size = length(ddr), weight = wt_df$prop)
  return(wt_sample$bcr_patient_barcode)
}

### ADD SUBTYPE AND STAGE #########################################################################
add_subtype_and_stage <- function(df, anno, subtype = TRUE, stage = TRUE, grade = TRUE, stage.type = "pathologic") {
  # read in subtype data
  add_anno_df <- read.csv(
    "~/GermlineSomaticAssociations/genome-wide/output/somatic_gwas/pancan/driver_rate_comparison/TCGA_gyn_supplementary_table4.csv",
    as.is = TRUE,
    skip = 1
  )
  if (subtype) {
    subtype_df <- add_anno_df[, c("Sample.ID", "BRCA_Subtype_PAM50")]
    colnames(subtype_df) <- c("bcr_patient_barcode", "subtype")
    # merge subtype with df
    df <- merge(df, subtype_df, by = "bcr_patient_barcode")
  } else {
    # adding NA columns to keep column names consistent
    df$subtype <- NA
  }
  if (stage) {
    # extract stage
    stage_df <- anno[, c("bcr_patient_barcode", paste0(stage.type, "_stage"))]
    colnames(stage_df) <- c("bcr_patient_barcode", "stage")
    # merge with df
    df <- merge(df, stage_df, by = "bcr_patient_barcode")
  } else {
    # adding NA columns to keep column names consistent
    df$stage <- NA
  }
  if (grade) {
    grade_df <- add_anno_df[, c("Sample.ID", "Tumor_Grade")]
    colnames(grade_df) <- c("bcr_patient_barcode", "grade")
    # merge subtype with df
    df <- merge(df, grade_df, by = "bcr_patient_barcode")
  } else {
    df$grade <- NA
  }
  return(df)
}

### REFROMAT STAGE ################################################################################
reformat_stage <- function(df) {
  # convert roman numerals to numbers
  roman <- 1:5
  names(roman) <- c("I", "II", "III", "IV", "X")
  # reformat stage
  df[df$stage == "[Not Available]", "stage"] <- NA
  df$stage <- gsub("Stage |A|B|C", "", df$stage)
  df$stage <- roman[df$stage]
  return(df)
}

### MUTATION RATE CNA AND SNV #####################################################################
calculate_mut_rate_cna_snv_meth_clones <- function(anno) {
  mut_rate_snv <- get_mutation_rate(type = "snv", anno = anno)
  mut_rate_snv$mutation <- "SNV"
  mut_rate_cna <- get_mutation_rate(type = "cna", anno = anno)
  mut_rate_cna$mutation <- "CNA"
  mut_rate_del <- get_mutation_rate(type = "deletion", anno = anno)
  mut_rate_del$mutation <- "CNA DEL"
  mut_rate_amp <- get_mutation_rate(type = "amplification", anno = anno)
  mut_rate_amp$mutation <- "CNA GAIN"
  mut_rate_indel <- get_mutation_rate(type = "indel", anno = anno)
  mut_rate_indel$mutation <- "INDEL DEL"
  mut_rate_insert <- get_mutation_rate(type = "insertion", anno = anno)
  mut_rate_insert$mutation <- "INSERT"
  # mut_rate_meth 			<- get_mutation_rate(type = 'meth', anno = anno)
  # mut_rate_meth$mutation 	<- 'METH'
  # mut_rate_clones 		<- get_mutation_rate(type = 'clones', anno = anno)
  # mut_rate_clones$mutation <- 'CLONES'
  # mut_rate 				<- rbind(mut_rate_snv, mut_rate_cna, mut_rate_meth, mut_rate_clones)
  mut_rate <- rbind(
    mut_rate_snv, mut_rate_cna, mut_rate_del,
    mut_rate_amp, mut_rate_indel, mut_rate_insert
  )
  return(mut_rate)
}

### FIND BRCA1 and BRCA2 MUTATED SAMPLES ##########################################################
find_brca1_brca2_samples <- function(cancer) {
  ddr <- list()
  ddr[["brca1"]] <- var_anno[var_anno$HUGO_Symbol == "BRCA1" & var_anno$cancer == cancer, "bcr_patient_barcode"]
  ddr[["brca2"]] <- var_anno[var_anno$HUGO_Symbol == "BRCA2" & var_anno$cancer == cancer, "bcr_patient_barcode"]
  ddr[["both"]] <- intersect(ddr$brca1, ddr$brca2)
  return(ddr)
}

### GET KRUSKAL P-VALUE ###########################################################################
get_kruskal_pvalue <- function(plot_data, group, value) {
  pvalue <- kruskal.test(
    plot_data[plot_data$group == group, value],
    factor(plot_data[plot_data$group == group, "variant"])
  )$p.value
  # convert to scientific notation
  sci <- scientific.notation(pvalue, digits = 2)
  return(sci)
}

### FIT LINEAR REG ################################################################################
fit_linear_reg <- function(plot_data, group, value, variant) {
  # subset plot data down to group and variant
  df <- plot_data[plot_data$group == group & plot_data$variant %in% c("WT", variant), ]
  # take log10 of counts and rate
  df$age_at_initial_pathologic_diagnosis <- as.numeric(df$age_at_initial_pathologic_diagnosis)
  df$log10count <- log10(df$count + 1)
  df$log10rate <- log10((df$count + 1) / df$age_at_initial_pathologic_diagnosis)
  # relevel factors
  df$variant <- relevel(factor(df$variant), ref = "WT")
  # fit model depending on cancer type
  if (grepl("BRCA", group)) {
    fit <- lm(
      as.formula(paste0("log10", value, " ~ variant + stage + subtype")),
      data = df
    )
  } else if (grepl("OV", group)) {
    fit <- lm(
      as.formula(paste0("log10", value, " ~ variant + stage + grade")),
      data = df
    )
  }
  # extract beta and pvalue
  beta <- coef(fit)[paste0("variant", variant)]
  pvalue <- summary(fit)$coefficients[paste0("variant", variant), 4]
  ci <- confint(fit)
  stats <- c(beta, pvalue, ci[paste0("variant", variant), 1], ci[paste0("variant", variant), 2])
  names(stats) <- c("beta", "pvalue", "l95", "u95")
  return(stats)
}

### ADD LINEAR REG KEY ############################################################################
add_linear_reg_key <- function(statsdf, group, value, x = 0.02, y = 0.98) {
  # # apply linear model
  # brca1_stats <- fit_linear_reg(plot_data, group, value, 'BRCA1')
  # brca2_stats <- fit_linear_reg(plot_data, group, value, 'BRCA2')
  statsdf <- statsdf[statsdf$group == group & statsdf$value == value, ]
  # format statsdf
  beta_brca1 <- round(statsdf[statsdf$variant == "BRCA1", "beta"], 2)
  l95_brca1 <- round(statsdf[statsdf$variant == "BRCA1", "l95"], 2)
  u95_brca1 <- round(statsdf[statsdf$variant == "BRCA1", "u95"], 2)
  fdr_brca1 <- scientific.notation(statsdf[statsdf$variant == "BRCA1", "fdr"], digits = 2, type = "list")
  beta_brca2 <- round(statsdf[statsdf$variant == "BRCA2", "beta"], 2)
  l95_brca2 <- round(statsdf[statsdf$variant == "BRCA2", "l95"], 2)
  u95_brca2 <- round(statsdf[statsdf$variant == "BRCA2", "u95"], 2)
  fdr_brca2 <- scientific.notation(statsdf[statsdf$variant == "BRCA2", "fdr"], digits = 2, type = "list")
  # create key
  key <- list(
    text = list(
      lab = c(
        as.expression(substitute(
          beta[vrt] * "= " * beta1,
          list(vrt = "BRCA1", beta1 = paste0(beta_brca1, " (", l95_brca1, ",", u95_brca1, ")"))
        )),
        as.expression(substitute(
          "FDR"[vrt] * "= " * base1 %*% 10^exponent1,
          list(vrt = "BRCA1", base1 = fdr_brca1[[1]], exponent1 = fdr_brca1[[2]])
        )),
        as.expression(substitute(
          beta[vrt] * "= " * beta2,
          list(vrt = "BRCA2", beta2 = paste0(beta_brca2, " (", l95_brca2, ",", u95_brca2, ")"))
        )),
        as.expression(substitute(
          "FDR"[vrt] * "= " * base2 %*% 10^exponent2,
          list(vrt = "BRCA2", base2 = fdr_brca2[[1]], exponent2 = fdr_brca2[[2]])
        ))
      )
    ),
    x = x,
    y = y,
    cex = 1.5,
    padding.text = 3,
    corner = c(0, 1)
  )
  return(key)
}
