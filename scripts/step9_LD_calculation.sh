
#chromosome 16 lead SNP

#convert .bgen chr 16 to bed bim fam for PLINK LD calculation
cd ./GWAS/ALSPAC/step9_LD

############################
#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=06:00:00
cd $PBS_O_WORKDIR

cd ./GWAS/ALSPAC/step9_LD/step1

module add apps/qctool-2.0

qctool -g ./GWAS/ALSPAC/step4_filter_SNPs/output/step2_SNP2KEEP_filter/output/data_chr16_sample_filter_subsetted.bgen \
-og ./GWAS/ALSPAC/step9_LD/step1/output/data_chr16_sample_filter_subsetted.bed  \
-ofiletype binary_ped \
-incl-rsids ./GWAS/ALSPAC/step3_create_SNP_stats/output/data_chr16_sample_filter_snpstats.txt \
-s ./GWAS/ALSPAC/step1_phenofile/pheno_file_include-list-only.txt
###########################

#add chr column to table
R
data <- read.table("data_chr16_sample_filter_subsetted.bim")
head(data)
data$chr <- 1
head(data)
data <- data[,c(6,1,2,3,4,5)]
head(data)
write.table(data, "data_chr16_sample_filter_subsetted.bim", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


###########################
#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:30:00
cd $PBS_O_WORKDIR

cd ./GWAS/ALSPAC/step9_LD/step2

module add apps/plink2

plink --bfile ./GWAS/ALSPAC/step9_LD/step1/output/data_chr16_sample_filter_subsetted \
--r2 \
--ld-snp rs13337037 \
--ld-window-kb 1000 \
--ld-window 99999 \
--ld-window-r2 0 \
--allow-extra-chr 
###########################

mv plink.ld data_chr16_sample_filter_subsetted.ld
mv plink.log data_chr16_sample_filter_subsetted.ld.log
mv plink.nosex data_chr16_sample_filter_subsetted.nosex



# chromosome 9 lead SNP

#convert .bgen chr 16 to bed bim fam for PLINK LD calculation
cd ./GWAS/ALSPAC/step9_LD

############################
#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=06:00:00
cd $PBS_O_WORKDIR

cd ./GWAS/ALSPAC/step9_LD/step1

module add apps/qctool-2.0

qctool -g ./GWAS/ALSPAC/step4_filter_SNPs/output/step2_SNP2KEEP_filter/output/data_chr09_sample_filter_subsetted.bgen \
-og ./GWAS/ALSPAC/step9_LD/step1/output/data_chr09_sample_filter_subsetted.bed  \
-ofiletype binary_ped \
-incl-rsids ./GWAS/ALSPAC/step3_create_SNP_stats/output/data_chr09_sample_filter_snpstats.txt \
-s ./GWAS/ALSPAC/step1_phenofile/pheno_file_include-list-only.txt
###########################

#add chr column to table
R
data <- read.table("data_chr09_sample_filter_subsetted.bim")
head(data)
data$chr <- 1
head(data)
data <- data[,c(6,1,2,3,4,5)]
head(data)
write.table(data, "data_chr09_sample_filter_subsetted.bim", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


###########################
#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:30:00
cd $PBS_O_WORKDIR

cd ./GWAS/ALSPAC/step9_LD/step2

module add apps/plink2

plink --bfile ./GWAS/ALSPAC/step9_LD/step1/output/data_chr09_sample_filter_subsetted \
--r2 \
--ld-snp rs10991823 \
--ld-window-kb 1000 \
--ld-window 99999 \
--ld-window-r2 0 \
--allow-extra-chr 
###########################

mv plink.ld data_chr09_sample_filter_subsetted.ld
mv plink.log data_chr09_sample_filter_subsetted.ld.log
mv plink.nosex data_chr09_sample_filter_subsetted.nosex



