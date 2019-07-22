###########################
## script to run meta analysis, clean up output, and save files.
###########################

#meta analysis performed using METAL

cd ./meta_analysis/

module add apps/metal-2011-3-25 

metal

#ALSPAC & NFBC

SCHEME STDERR
STDERR SE
SEPARATOR TAB
MARKER SNP
ALLELE EA NEA
EFFECT BETA
FREQLABEL EAF
PVALUE P
AVERAGEFREQ ON
MINMAXFREQ ON

PROCESS ./GWAS/ALSPAC/step5_GWAS/output/GWAS_output_final.txt

PROCESS ./GWAS/NFBC/NFBC_GWAS_final_parental.txt

OUTFILE ./meta_analysis/meta_analysis_ALSPAC_NFBC .txt

ANALYZE HETEROGENEITY

CLEAR

#meta analysis clean up

cd ./meta_analysis

module add languages/R-3.5-ATLAS-gcc-7.1.0
R

# load GWAS_outut
library(data.table)
data <- fread("meta_analysis_ALSPAC_NFBC1.txt", header = T)
head(data)

colnames(data) <- c("SNP", "EA", "NEA", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "Effect", "SE", "P", "Direction", "HetISq", "HetChiSq", "HetDf", "HetPVal")
head(data)

write.table(data, "meta_analysis_ALSPAC_NFBC1.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

#############
# remove all ? Direction SNPs 

#test on small random sub sample
#data_random <- data[sample(nrow(data), 100), ]

#do whole dataframe
data$qmark <- lapply(data$Direction, function(x)
  if(x== "?+") { 
    1
  }
  else if(x== "+?") {
    1
  }
  else if(x== "?-") {
    1
  } 
  else if(x== "-?") { 
    1
  } 
  else if(x== "??") {
    1
  } 
  else 0
)

dim(data) #11349282
data_new <- subset(data, qmark == 0)
dim(data_new) #7802683
data_new <- data_new[,c(1:15)]

write.table(data_new, "meta_analysis_final.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data_sig <- subset(data_new, P <= 5e-8) 
dim(data_sig) #54
data_sig <- data_sig[order(P),] 
head(data_sig)

write.table(data_sig, "meta_analysis_final_sig.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

#############
# get chr and bp for all rsID for meta-analysis from the snp-stats files

#create combined snp-stats file from step3
head -9 ./GWAS/ALSPAC/step3_create_SNP_stats/output/data_chr01_sample_filter_snpstats.txt | tail -n +9 > ./meta_analysis/snp-stats_combined_rsid_chr_bp.txt
tail -n +10 -q ./GWAS/ALSPAC/step3_create_SNP_stats/output/*snpstats.txt >> ./meta_analysis/snp-stats_combined_rsid_chr_bp.txt
awk '{ print $1, $2, $3 }' snp-stats_combined_rsid_chr_bp.txt > snp-stats_combined_rsid_chr_bp_new.txt
rm snp-stats_combined_rsid_chr_bp.txt
mv snp-stats_combined_rsid_chr_bp_new.txt snp-stats_combined_rsid_chr_bp.txt

#cut columns 1,2,3 and rearrange to SNP, chr, position
cd ./meta_analysis

R

rm(list=ls())

library(data.table)
data <- fread("./meta_analysis/meta_analysis_final.txt", header = T)
head(data)
dim(data) #7802683

snp_stats <- fread("./meta_analysis/snp-stats_combined_rsid_chr_bp.txt", header = T)
head(snp_stats)
dim(snp_stats) #28699532
snp_stats <- snp_stats[,c(2,1,3)]
head(snp_stats)
colnames(snp_stats) <- c("SNP", "chr", "position")
head(snp_stats)

library(dplyr)
data_new <- left_join(data, snp_stats, by = "SNP")
dim(data_new) #7802718
head(data_new)

#remove duplicate SNPs
data_final <- data_new[!duplicated(data_new$SNP), ]
dim(data_final) #7802683
head(data_final)

#calculate OR and CI
data_final$OR <- exp(data_final$Effect)
data_final$ci_lower <- data_final$OR - 1.96 * data_final$SE
data_final$ci_upper <- data_final$OR + 1.96 * data_final$SE
dim(data_final)
head(data_final)

#rearrange columns for ease
names(data_final)
data_final <- data_final[,c(1:17,20:22)]
head(data_final)

write.table(data_final, "meta_analysis_final.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data_sig <- subset(data_final, P <= 5e-8) 
dim(data_sig) #54
head(data_sig)
names(data_sig)
data_sig <- data_sig[order(data_sig$P),] 
head(data_sig)
write.table(data_sig, "meta_analysis_final_sig.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

