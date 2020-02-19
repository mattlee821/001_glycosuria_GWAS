###########################
## script to run meta analysis, clean up output, and save files.
###########################

#meta analysis performed using METAL

cd meta_analysis/

module add languages/R-3.5-ATLAS-gcc-7.1.0
R
library(data.table)

# create sample size weighting column
alspac <- fread("../GWAS/ALSPAC/step5_GWAS/output/GWAS_output_final.txt", header = T)
names(alspac)
alspac$N <- 4/(1/alspac$cases_total + 1/alspac$controls_total)
summary(alspac$N)

nfbc <- fread("../GWAS/NFBC/NFBC_GWAS_final_parental.txt", header = T)
names(nfbc)
nfbc$N <- 4/(1/747 + 1/2991)
summary(nfbc$N)

write.table(alspac, "ALSPAC_GWAS_meta-analysis_prep.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(nfbc, "NFBC_GWAS_meta-analysis_prep.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

q()
n

module add apps/metal-2011-3-25 

metal

#ALSPAC & NFBC

SCHEME SAMPLESIZE
WEIGHT N
STDERR SE
SEPARATOR TAB
MARKER SNP
ALLELE EA NEA
EFFECT BETA
FREQLABEL EAF
PVALUE P
AVERAGEFREQ ON
MINMAXFREQ ON

PROCESS ALSPAC_GWAS_meta-analysis_prep.txt

PROCESS NFBC_GWAS_meta-analysis_prep.txt

OUTFILE meta_analysis_ALSPAC_NFBC_sample-size .txt

ANALYZE HETEROGENEITY

QUIT

#meta analysis clean up
R

# load GWAS_outut
library(data.table)
data <- fread("meta_analysis_ALSPAC_NFBC_sample-size1.txt", header = T)
head(data)

colnames(data) <- c("SNP", "EA", "NEA", "Freq1", "FreqSE", "MinFreq", "MaxFreq", "Weight", "Zscore", "P", "Direction", "HetISq", "HetChiSq", "HetDf", "HetPVal")
head(data)

write.table(data, "meta_analysis_ALSPAC_NFBC_sample-size1.txt", 
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
dim(data_new) #7802690
data_new <- data_new[,c(1:15)]

write.table(data_new, "meta_analysis_final_sample-size.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data_sig <- subset(data_new, P <= 5e-8) 
dim(data_sig) #45
data_sig <- data_sig[order(P),] 
head(data_sig)

write.table(data_sig, "meta_analysis_final_sample-size_sig.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
q()
n
#############
# get chr and bp for all rsID for meta-analysis from the snp-stats files

#create combined snp-stats file from step3
head -9 ../GWAS/ALSPAC/step3_create_SNP_stats/output/data_chr01_sample_filter_snpstats.txt | tail -n +9 > snp-stats_combined_rsid_chr_bp.txt
tail -n +10 -q ../GWAS/ALSPAC/step3_create_SNP_stats/output/*snpstats.txt >> snp-stats_combined_rsid_chr_bp.txt
awk '{ print $1, $2, $3 }' snp-stats_combined_rsid_chr_bp.txt > snp-stats_combined_rsid_chr_bp_new.txt
rm snp-stats_combined_rsid_chr_bp.txt
mv snp-stats_combined_rsid_chr_bp_new.txt snp-stats_combined_rsid_chr_bp.txt

#cut columns 1,2,3 and rearrange to SNP, chr, position
R

rm(list=ls())

library(data.table)
data <- fread("meta_analysis_final_sample-size.txt", header = T)
head(data)
dim(data) #7802690

snp_stats <- fread("snp-stats_combined_rsid_chr_bp.txt", header = T)
head(snp_stats)
dim(snp_stats) #28699532
snp_stats <- snp_stats[,c(2,1,3)]
head(snp_stats)
colnames(snp_stats) <- c("SNP", "chr", "position")
head(snp_stats)

library(dplyr)
data_new <- left_join(data, snp_stats, by = "SNP")
dim(data_new) #7802725
head(data_new)

#remove duplicate SNPs
data_final <- data_new[!duplicated(data_new$SNP), ]
dim(data_final) #7802690
head(data_final)

write.table(data_final, "meta_analysis_final_sample-size.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data_sig <- subset(data_final, P <= 5e-8) 
dim(data_sig) #54
head(data_sig)
names(data_sig)
data_sig <- data_sig[order(data_sig$P),] 
head(data_sig)
write.table(data_sig, "meta_analysis_final_sig_sample-size.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

