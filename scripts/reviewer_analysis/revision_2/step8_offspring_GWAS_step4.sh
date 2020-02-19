##########################
## step9
## Filter Variants for INFO and MAF
##########################

## USE snp-stats data file to filter variants on INFO (>=0.3) and MAF (>=0.01). DO THIS IN R.

##########################
#1. run filter interactively on BC3
##########################
qsub -I -l walltime=6:00:00 -l nodes=1:ppn=1

cd ./001_glycosuria/reviewer_analysis/GWAS_revision2/step8

mkdir output
mkdir output/filter
module add languages/R-3.5-ATLAS-gcc-7.1.0

R

rm(list=ls())

datadir = "./001_glycosuria/reviewer_analysis/GWAS_revision2/step7/output/"
files = list.files(datadir)
w = grep("_snpstats.txt", files)

#### 
dataOUTdir = "./001_glycosuria/reviewer_analysis/GWAS_revision2/step8/output/filter/"

for(file in files){
cat(paste0(file, "\n"))
n = paste0(datadir, file)
snpstats <- read.table(n, header = F, skip = 9, as.is = TRUE)
dim(snpstats)

colnames(snpstats) <- c("chromosome", "rsid", "position", "alleleA", "alleleB", "comment", "HW_exact_p_value", "HW_lrt_p_value", "alleleA_count", "alleleB_count", "alleleA_frequency", "alleleB_frequency", "minor_allele_frequency", "minor_allele", "major_allele", "info", "impute_info", "missing_proportion", "A", "B", "AA", "AB", "BB", "NULL", "total")

snpstats_info <- subset(snpstats, info >= 0.3)

summary(snpstats_info$info)
dim(snpstats_info)

snpstats_info_MAF <- subset(snpstats_info, minor_allele_frequency >= 0.01)

summary(snpstats_info_MAF$minor_allele_frequency)
dim(snpstats_info_MAF)

##. file nameinis: data_chr01_sample_filter_snpstats.txt

filename_out = paste0( gsub(".txt","", file), "_MAF_info.txt"  )
nout = paste0(dataOUTdir, filename_out )

write.table(snpstats_info_MAF, file = nout,
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

}

rm(list=ls())
##########################

#move output to outpout/step1
##########################
# then for each chromosome (include chr# in name of file) make a list of SNPIDs to keep. One SNP on each line.
##########################

## do this in command line
# cut first column from file and create txt file per chromosome of SNP IDs to keep that we then use in the next step
cd ./001_glycosuria/reviewer_analysis/GWAS_revision2/step8/output/filter/
mkdir SNP2Keep
for i in *_snpstats*.txt; do echo ${i}; tail -n -2 ${i} | cut -f 2 ${i} > SNP2Keep/${i%.txt}_SNP2Keep.txt; done

cd ./001_glycosuria/reviewer_analysis/GWAS_revision2/step8/output

##########################
#2. create new bgen data set with newly filtered SNP file
##########################
mkdir step2_SNP2KEEP_filter
cd step2_SNP2KEEP_filter

#create master_sub qsub script as template
###########################
#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -e log/
#PBS -o log/
cd $PBS_O_WORKDIR

module add apps/qctool-2.0

cd ./001_glycosuria/reviewer_analysis/GWAS_revision2/step8/output/step2_SNP2KEEP_filter

dir1=./001_glycosuria/reviewer_analysis/GWAS_revision2/step6/output/
dir2=./001_glycosuria/reviewer_analysis/GWAS_revision2/step8/output/filter/SNP2Keep/

VAR1=data_chr01_sample_filter.bgen
VAR2=data_chr01_sample_filter_snpstats_MAF_info_SNP2Keep.txt

qctool -g ${dir1}${VAR1} -og ${VAR1%.bgen}_subsetted.bgen -incl-rsids ${dir2}${VAR2} -log ${VAR1%.bgen}_subsetted_log.txt
############################

#create .txt file of all names for submission scripts based on files in a directory
ls -1 ./001_glycosuria/reviewer_analysis/GWAS_revision2/step6/output > subset_filenames.txt

#create multiple qsub scripts based on master_sub
cat subset_filenames.txt | while read i; do echo ${i}; awk '{ if (NR == 16) print "VAR1='${i}'"; else print $0}'  master_sub > ${i}.sh; done

#change VAR2 manually in each of the submission files

#
mkdir log

#create run file to run all qsub scripts in one
ls |grep .bgen.sh |sed 's/\ / /g' |awk '{print "qsub "$i" -e log -o log"}' > run

#rename run to run.sh
mv run run.sh

#run run.sh file
sh run.sh

#make run directory and move all .sh files there
mkdir run

#clean up directory
mv *sh run/
mv subset_filenames.txt run/
mkdir output/
mv *bgen output/
mv *log.txt log/
mv master_sub run/

