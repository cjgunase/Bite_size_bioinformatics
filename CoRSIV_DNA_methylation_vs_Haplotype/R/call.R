library(vcfR,quietly = T)
library(tidyr,quietly = T)
library(stringr,quietly = T)

args = commandArgs(trailingOnly=TRUE)

CHR = args[1]

CoRSIV_BED <- read.table("../../GLOBAL_BED_FILES/combined_corsivs.bed")
CoRSIV_Haplotype_BED <- read.table("../LOCAL_BED_FILES/haplotypes_over_CoRSIVs.bed")
methylation_data <- read.csv("../../../CoRSIV_Capture/USC_pediatric_glioma_CoRSIV_Methylation_Data_10x_Depth_merged_reads.csv")
#match the column name in DNA methylation table with genotype files
colnames(methylation_data) <- str_replace(
    str_replace(colnames(methylation_data),'_control',""),'_case',"")
sample_id_table1 <- read.table("../misc_files/sample_id_table1.txt",header=T)
sample_id_table2 <- read.table("../misc_files/subjectIDs.txt",header=T)
SAMPLE_IDs<- merge(sample_id_table1,sample_id_table2,by.x = "ID2",by.y = "blindid")

vcf_chr <- read.vcfR(file = paste0("../VCF/Genotype_post_imputation/",CHR,".subset.dose.vcf.gz"))
gt_chr <- extract.gt(vcf_chr,element = 'GT')
SAMPLE_IDs <- SAMPLE_IDs[match(colnames(gt_chr), SAMPLE_IDs$ID),]#This will endup with sample IDs that will be used in this analysis
colnames(gt_chr) <- SAMPLE_IDs$ID1
gt_chr <- data.frame(gt_chr)
gt_chr$snp_id <- rownames(gt_chr)
gt_chr <- separate(data=gt_chr,col=snp_id,into = c("chr","loc","REF","ALT"),sep = ":" )

haplotype_CoRSIVs <- CoRSIV_Haplotype_BED[CoRSIV_Haplotype_BED$V1==CHR,]
haplotype_CoRSIVs <- haplotype_CoRSIVs[!duplicated(haplotype_CoRSIVs$V10),]
r_values <- c()
i=1
for (CoRSIV_ID in unique(as.character(haplotype_CoRSIVs$V10))){
    print(paste0(CHR,"-",i))
    snps_unformatted <- strsplit(as.character(haplotype_CoRSIVs[haplotype_CoRSIVs$V10==CoRSIV_ID,]$V6),split = "\\|")[[1]]
    snps_formated <- str_replace_all(str_replace(snps_unformatted, "_b38", ""),"_",":")
    snps_formated <- intersect(rownames(gt_chr),snps_formated)
    snps_gt_data <- gt_chr[rownames(gt_chr) %in% snps_formated,]
    snps_gt_data <- snps_gt_data[,1:11]
    
    #sum(as.numeric(unlist(strsplit(as.character(snps_gt_data[,1]),"\\|"))))
    
    allele_sum <- as.numeric(sapply(snps_gt_data, function(x) sum(as.numeric(unlist(strsplit(as.character(x),"\\|"))))))
    
    methy_vector <- as.numeric(methylation_data[methylation_data$CoRSIV_ID==CoRSIV_ID,][colnames(snps_gt_data)[1:11]])
    
    r_values <- c(r_values,cor(allele_sum,methy_vector,method = "pearson"))
    i=i+1
}


write.table(file = paste0("spearman_R_",CHR,".txt"),data.frame(unique(as.character(haplotype_CoRSIVs$V10)),r_values),row.names = FALSE,quote = F)





