##########################
## step5
## GWAS and clean up
##########################

#######################################################################################################################################

##########################
#run test GWAS on single chromosome
##########################

###########################
#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
cd $PBS_O_WORKDIR

module add apps/snptest.2.5.2

cd ./001_glycosuria/reviewer_analysis/GWAS_revision2/step5

#run test gwas on single CHR to test
snptest_v2.5.2 -data ./001_glycosuria/reviewer_analysis/GWAS_revision2/step4/output/step2_SNP2KEEP_filter/output/data_chr16_sample_filter_subsetted.bgen ./001_glycosuria/reviewer_analysis/GWAS_revision2/pheno_file/pheno_file_include-list-only.txt -o data_chr16.out -log data_chr16.log -frequentist 1 -method expected -pheno t6001 -cov_names PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 -use_raw_phenotypes
##########################	

# -data <a> <b>...: specify data files for analysis in .gen and .sample pairs. Automatic detection of .gz files.
# -log <a>: name of log file.
# -o <a>: name of output file.
# -frequentist <a> <b>...: specify which Frequentist tests to fit.
# -method <a>: method used to fit model, this can be one of "threshold", "expected", "score", "ml", "newml", or "em".
# -pheno <a>: specify name of phenotype to use.
# -cov_names <a> <b>...: list names of covariates to use.
# -use_raw_phenotypes: Do not normalise continuous phenotypes to have mean 0 and variance 1.

#######################################################################################################################################


#######################################################################################################################################

##########################
# run GWAS
##########################


#create master_sub qsub script as a template
###########################
#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -e log/
#PBS -o log/
cd $PBS_O_WORKDIR

VAR1=data_chr01_sample_filter_subsetted.bgen

module add apps/snptest.2.5.2

cd ./001_glycosuria/reviewer_analysis/GWAS_revision2/step5

snptest_v2.5.2 -data ./001_glycosuria/reviewer_analysis/GWAS_revision2/step4/output/step2_SNP2KEEP_filter/output/${VAR1} ./001_glycosuria/reviewer_analysis/GWAS_revision2/pheno_file/pheno_file_include-list-only.txt -o ./001_glycosuria/reviewer_analysis/GWAS_revision2/step5/output/${VAR1%.bgen}.out -log ${VAR1%.bgen}.log -frequentist 1 -method expected -pheno t6001 -cov_names PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 -use_raw_phenotypes
##########################	

#create .txt file of all names for submission scripts based on files in a directory
ls -1 ./001_glycosuria/reviewer_analysis/GWAS_revision2/step4/output/step2_SNP2KEEP_filter/output > data_file_names.txt

#edit the ouput from this to have a single chr per line

#create multiple qsub scripts based on master_sub
cat data_file_names.txt | while read i; do echo ${i}; awk '{ if (NR == 9) print "VAR1='${i}'"; else print $0}' master_sub > ${i}.sh; done

#
mkdir log
mkdir output
#create run file to run all qsub scripts in one
ls |grep .sh |sed 's/\ / /g' |awk '{print "qsub "$i" -e log -o log"}' > run

#rename run to run.sh
mv run run.sh

#run run.sh file
sh run.sh

mkdir run
mv *sh run/
mv data_file_names.txt run/
mv *.log log/

###########################
# cat files into final GWAS output
# each chromosome file has a header of 13 lines that need skipping and the actual header is on line 14
# make a header file for GWAS output
head -14 data_chr01_sample_filter_subsetted.out | tail -n +14 > GWAS_header.txt

#make GWAS file of all chr with header
head -14 data_chr01_sample_filter_subsetted.out | tail -n +14 > GWAS_output_raw.txt

tail -n +15 -q data* >> GWAS_output_raw.txt

#remove comment column (46) from file GWAS_output_raw.txt
sed -i -r 's/\S+//46' GWAS_output_raw.txt

#check wc -l - there should be 321 lines less in GWAS_output_raw.txt
wc -l *out 
wc -l GWAS_output_raw.txt 
###################################################################################################

###########################
# tidy final GWAS file up for making plots
cd ./001_glycosuria/reviewer_analysis/GWAS_revision2/step5/output

module add languages/R-3.5-ATLAS-gcc-7.1.0
R

# load GWAS_outut
library(data.table)
data <- fread("GWAS_output_raw.txt", header = T, fill = T)
head(data)
dim(data)
summary(data$frequentist_add_pvalue)

# create new data frame with only columns of interest
data_new <- data[,c(1,2,4,5,6,9,18,23,28,29,30,31,39,40,41,42,44,45)]
head(data_new)
str(data_new)
dim(data_new)
names(data_new)

summary(data_new$frequentist_add_pvalue)
# remove all rows where NA appears in the p-value column
dim(data_new)
data_new <- data_new[complete.cases(data_new[ , frequentist_add_pvalue]),]
dim(data_new)
summary(data_new$frequentist_add_pvalue)

# remove all rows where '.' appears in rsID column
dim(data_new)
data_new <- subset(data_new, rsid != ".")
dim(data_new)

# change column order
names(data_new)
colnames(data_new) <- c("chr", "SNP", "position", "NEA", "EA", "info", "all_total", "cases_total", "controls_total", "EAF", "cases_maf", "controls_maf", "all_OR", "all_OR_lower", "all_OR_upper", "P", "BETA", "SE")
names(data_new)
data_new2 <- data_new[,c("SNP", "chr", "position", "EA", "NEA", "info", "all_total", "cases_total", "controls_total", "EAF", "cases_maf", "controls_maf", "all_OR", "all_OR_lower", "all_OR_upper", "P", "BETA", "SE")]

dim(data_new2) 
summary(data_new2$info)
data_new2 <- subset(data_new2, info >= 0.3)
summary(data_new2$info)
dim(data_new2) 

summary(data_new2$EAF)
data_new2 <- subset(data_new2, EAF >= 0.01)
dim(data_new2) 

summary(data_new2$EAF)
summary(data_new2$info)

### final check of data - ensure values make sense
summary(data_new2$all_total)
summary(data_new2$cases_total) 
summary(data_new2$controls_total)

dim(data_new2)
summary(data_new2$P)
summary(data_new2$info)
summary(data_new2$EAF)

write.table(data_new2, "GWAS_output_final.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

###### SNPs reaching GWAS threshol 5e-8
data_sig <- subset(data_new2, P <= 5e-8)
dim(data_sig) 
head(data_sig)
data_sig <- data_sig[order(P),] 
head(data_sig)

write.table(data_sig, "GWAS_output_final_sig.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")








