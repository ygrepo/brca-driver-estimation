##### plot_drivergram.R #####
# Kuan-lin Huang 2019

# # set work dir for testing in dev environ
source("../global_aes_out.R")

somatic_f = "../../../../Huang_lab_data/TCGA_PanCanAtlas_2018/MC3_Ellrott_CellSys2018/mc3.v0.2.8.PUBLIC.maf.gene_vclass_HGVSp_sample.gz"
somatic = read.table(header=T, quote = "", sep="\t", file = gzfile(somatic_f), stringsAsFactors=FALSE)
somatic$bcr_sample_barcode = substr(somatic$Tumor_Sample_Barcode,1,12)

somatic_likelyfunctional_driver_f = "../../../../Huang_lab_data/TCGA_PanCanAtlas_2018/MC3_Ellrott_CellSys2018/mc3.v0.2.8.PUBLIC.maf.gene_vclass_HGVSp_sample_likelyDriver.tsv"
somatic_likelyfunctional_driver = read.table(header=T, quote = "", sep="\t", file = somatic_likelyfunctional_driver_f, stringsAsFactors=FALSE)
# only include one sample per patient
somatic_likelyfunctional_driver_samples = somatic_likelyfunctional_driver[!duplicated(somatic_likelyfunctional_driver$Tumor_Sample_Barcode),c("Tumor_Sample_Barcode","bcr_patient_barcode")]  
somatic_likelyfunctional_driver_samples = somatic_likelyfunctional_driver_samples[order(somatic_likelyfunctional_driver_samples$Tumor_Sample_Barcode),]
non_duplicated_sample = somatic_likelyfunctional_driver_samples$Tumor_Sample_Barcode[!duplicated(somatic_likelyfunctional_driver_samples$bcr_patient_barcode)]
somatic_likelyfunctional_driver_uniq = somatic_likelyfunctional_driver[somatic_likelyfunctional_driver$Tumor_Sample_Barcode %in% non_duplicated_sample,]

fn = "../../../../Huang_lab_data/TCGA_PanCanAtlas_2018/Germline_Huang_Cell2018/PanCan_ClinicalData_V4_wAIM_filtered10389.txt"
clinical = read.table(sep="\t",header=T, quote="",stringsAsFactors = F, file=fn)

clinical_sample = clinical[clinical$bcr_patient_barcode %in% somatic$bcr_sample_barcode,]
count_per_type = data.frame(table(clinical_sample$type))
colnames(count_per_type) = c("Cancer","sample_size")

somatic_likelyfunctional_driver_uniq_clin = merge(somatic_likelyfunctional_driver_uniq,clinical[,c("bcr_patient_barcode","type")],by="bcr_patient_barcode")
driver_per_type = data.frame(table(somatic_likelyfunctional_driver_uniq_clin$Hugo_Symbol,somatic_likelyfunctional_driver_uniq_clin$type))
colnames(driver_per_type) = c("Gene","Cancer","Counts")

overall_counts = merge(driver_per_type,count_per_type,by="Cancer")
overall_counts$avg_count_per_patient = overall_counts$Counts/overall_counts$sample_size
overall_counts = overall_counts[overall_counts$avg_count_per_patient!=0,]

# plot pancan
top_genes = unique(overall_counts$Gene[order(overall_counts$avg_count_per_patient, decreasing = T)])[1:30]
overall_counts$Gene_summary = as.character(overall_counts$Gene)
overall_counts$Gene_summary[!(overall_counts$Gene %in% top_genes)] = "AALL_OTHER_GENES"
overall_counts = overall_counts[order(overall_counts$Gene_summary),]
p = ggplot(data=overall_counts,aes(y=avg_count_per_patient,x=Cancer,fill = Gene_summary))
p = p + facet_grid(.~Cancer,space="free",scale="free")
p = p + geom_col() + theme_bw() + scale_fill_manual(values = col_vector)
p = p + labs(x="Cancer",y= "Avg Count of Driver/Patient")
p 
fn = 'out/driverGram_PanCan.pdf'
ggsave(fn,w=18,h=6,useDingbat=F)

for (cancer in unique(overall_counts$Cancer)){
  overall_counts_c = overall_counts[overall_counts$Cancer == cancer,]
  # plot pancan
  top_genes = unique(overall_counts_c$Gene[order(overall_counts_c$avg_count_per_patient, decreasing = T)])[1:15]
  overall_counts_c$Gene_summary = as.character(overall_counts_c$Gene)
  overall_counts_c$Gene_summary[!(overall_counts_c$Gene %in% top_genes)] = "AALL_OTHER_GENES"
  p = ggplot(data=overall_counts_c,aes(y=avg_count_per_patient,x=Cancer,fill = Gene_summary))
  p = p + facet_grid(.~Cancer,space="free",scale="free")
  p = p + geom_col() + theme_bw() + scale_fill_manual(values = col_vector)
  p = p + labs(x="Cancer",y= "Avg Count of Driver/Patient")
  p 
  fn = paste('out/driverGram_',cancer,'.pdf',sep="")
  ggsave(fn,w=3,h=6,useDingbat=F)
}


fn = "/Users/huangk06/Box\ Sync/Huang_lab_kuan/control_access_data/TCGA_PanCanAtlas_2018/Germline_Huang_Cell2018/PCA_pathVar_integrated_filtered_adjusted_ancestry_PLPonly.tsv" # updated path for protected data
pathVar = read.table(sep="\t",header=T, quote="",stringsAsFactors = F, file=fn)