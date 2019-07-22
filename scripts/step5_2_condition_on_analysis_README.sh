############################
#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
cd $PBS_O_WORKDIR

cd ./GWAS/ALSPAC/step5_GWAS/conditional

module add apps/snptest.2.5.2

snptest_v2.5.2 -data ./GWAS/ALSPAC/step4_filter_SNPs/output/step2_SNP2KEEP_filter/output/data_chr16_sample_filter_subsetted.bgen ./GWAS/ALSPAC/step1_phenofile/pheno_file_include-list-only.txt -o ./GWAS/ALSPAC/step5_GWAS/conditional/data_chr16_sample_filter_subsetted_conditional.out -log data_chr16_sample_filter_subsetted_conditional.log -frequentist 1 -method expected -pheno e113 -condition_on rs13337037 -cov_names PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 -use_raw_phenotypes 
############################

cd ./GWAS/ALSPAC/step5_GWAS/conditional

module add languages/R-3.5-ATLAS-gcc-7.1.0
R

# load GWAS_outut
library(data.table)
data <- fread("data_chr16_sample_filter_subsetted_conditional.out", header = T, skip = 14, fill=T)
head(data)
dim(data)
summary(data$frequentist_add_pvalue)

# create new data frame with only columns of interest
data_new <- data[,c(1,2,4,5,6,9,18,23,28,29,30,31,39,40,41,42,44,45)]
head(data_new)

str(data_new)

# remove all rows where NA appears in the p-value column
data_new <- data_new[complete.cases(data[ , frequentist_add_pvalue]),]
dim(data_new)
# remove all rows where '.' appears in rsID column
data_new <- subset(data_new, rsid != ".")
dim(data_new)

#change column name for chromosome
head(data_new)
colnames(data_new)[1] <- "chr"


# change column order
colnames(data_new) <- c("chr", "SNP", "position", "NEA", "EA", "info", "all_total", "cases_total", "controls_total", "EAF", "cases_maf", "controls_maf", "all_OR", "all_OR_lower", "all_OR_upper", "P", "BETA", "SE")
data_new2 <- data_new[,c("SNP", "chr", "position", "EA", "NEA", "info", "all_total", "cases_total", "controls_total", "EAF", "cases_maf", "controls_maf", "all_OR", "all_OR_lower", "all_OR_upper", "P", "BETA", "SE")]

dim(data_new2) #278837
summary(data_new2$info)
data_new2 <- subset(data_new2, info >= 0.3)
summary(data_new2$info)
dim(data_new2) #278821

summary(data_new2$EAF)
data_new2 <- subset(data_new2, EAF >= 0.01)
nrow(data_new2) #278419

summary(data_new2$EAF)
summary(data_new2$info)

### final check of data - ensure values make sense
summary(data_new2$all_total)
summary(data_new2$cases_total) 
summary(data_new2$controls_total)


dim(data_new2) #278926
summary(data_new2$P)
summary(data_new2$info)
summary(data_new2$EAF)


write.table(data_new2, "conditional_analysis_rs13337037_output_final.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

###### SNPs reaching GWAS threshol 5e-8
data_sig <- subset(data_new2, P <= 5e-8)
dim(data_sig) #0
head(data_sig)

#write.table(data_sig, "rs13337037_sig.txt", 
#            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")




############################
#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
cd $PBS_O_WORKDIR

cd ./GWAS/ALSPAC/step5_GWAS/conditional

module add apps/snptest.2.5.2

snptest_v2.5.2 -data ./GWAS/ALSPAC/step4_filter_SNPs/output/step2_SNP2KEEP_filter/output/data_chr09_sample_filter_subsetted.bgen ./GWAS/ALSPAC/step1_phenofile/pheno_file_include-list-only.txt -o ./GWAS/ALSPAC/step5_GWAS/conditional/data_chr09_sample_filter_subsetted_conditional.out -log data_chr09_sample_filter_subsetted_conditional.log -frequentist 1 -method expected -pheno e113 -condition_on rs10991823 -cov_names PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 -use_raw_phenotypes 
############################

cd ./GWAS/ALSPAC/step5_GWAS/conditional

module add languages/R-3.5-ATLAS-gcc-7.1.0
R

# load GWAS_outut
library(data.table)
data <- fread("data_chr09_sample_filter_subsetted_conditional.out", header = T, skip = 14, fill=T)
head(data)
dim(data)
summary(data$frequentist_add_pvalue)

# create new data frame with only columns of interest
data_new <- data[,c(1,2,4,5,6,9,18,23,28,29,30,31,39,40,41,42,44,45)]
head(data_new)

str(data_new)

# remove all rows where NA appears in the p-value column
data_new <- data_new[complete.cases(data[ , frequentist_add_pvalue]),]
dim(data_new)
# remove all rows where '.' appears in rsID column
data_new <- subset(data_new, rsid != ".")
dim(data_new)

#change column name for chromosome
head(data_new)
colnames(data_new)[1] <- "chr"


# change column order
colnames(data_new) <- c("chr", "SNP", "position", "NEA", "EA", "info", "all_total", "cases_total", "controls_total", "EAF", "cases_maf", "controls_maf", "all_OR", "all_OR_lower", "all_OR_upper", "P", "BETA", "SE")
data_new2 <- data_new[,c("SNP", "chr", "position", "EA", "NEA", "info", "all_total", "cases_total", "controls_total", "EAF", "cases_maf", "controls_maf", "all_OR", "all_OR_lower", "all_OR_upper", "P", "BETA", "SE")]

dim(data_new2) #386945
summary(data_new2$info)
data_new2 <- subset(data_new2, info >= 0.3)
summary(data_new2$info)
dim(data_new2) #386912

summary(data_new2$EAF)
data_new2 <- subset(data_new2, EAF >= 0.01)
dim(data_new2) #386284

summary(data_new2$EAF)
summary(data_new2$info)

### final check of data - ensure values make sense
summary(data_new2$all_total)
summary(data_new2$cases_total) 
summary(data_new2$controls_total)


dim(data_new2) #386284
summary(data_new2$P)
summary(data_new2$info)
summary(data_new2$EAF)


write.table(data_new2, "conditional_analysis_rs10991823_output_final.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

###### SNPs reaching GWAS threshol 5e-8
data_sig <- subset(data_new2, P <= 5e-8)
dim(data_sig) #0
head(data_sig)

#write.table(data_sig, "rs13337037_sig.txt", 
#            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")








