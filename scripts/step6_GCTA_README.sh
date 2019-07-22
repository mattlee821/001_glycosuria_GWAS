#GCTA analysis
mkdir step6_GCTA
cd ./GWAS/ALSPAC/step6_GCTA/
mkdir PLINK

# GCTA-GERML: estimating of the phenotypic variance explained by the SNPs

# first convert the bgen subsetted files to PLINK

############################
#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=06:00:00
cd $PBS_O_WORKDIR
module add apps/qctool-2.0

cd ./GWAS/ALSPAC/step6_GCTA/PLINK

dir1=./GWAS/ALSPAC/step4_filter_SNPs/output/step2_SNP2KEEP_filter/output/
dir2=./GWAS/ALSPAC/step6_GCTA/PLINK/

VAR1=data_chr01_sample_filter_subsetted.bgen

qctool -g ${dir1}${VAR1} -og ${dir2}${VAR1%.bgen}.bed -ofiletype binary_ped -incl-samples ./GWAS/ALSPAC/step1_phenofile/include.txt -s ./GWAS/ALSPAC/step1_phenofile/pheno_file_include-list-only.txt
############################

#create .txt file of all names for submission scripts based on files in a directory
ls -1 ./GWAS/ALSPAC/step4_filter_SNPs/output/step2_SNP2KEEP_filter/output/ > VAR1.txt

#create multiple qsub scripts based on master_sub
cat VAR1.txt | while read i; do echo ${i}; awk '{ if (NR == 13) print "VAR1='${i}'"; else print $0}'  master_sub > ${i}.sh; done

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

mv *sh run
mv VAR1.txt run/


#####################################################
## Do the bim files have chromosome data in them?? Add chromosome data to bim file 
#####################################################
R

f = list.files(); k = grep("bim",f); f = f[k]
for(i in 2:length(f)){
	file = f[i]
	mydata = read.table(file, h = FALSE, as.is = TRUE, sep = " ")
	chr = rep(i,nrow(mydata))
	mydata[,1] = chr
	write.table(mydata, file = file, col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE)
}

#check that each .bim file has the chromosome added
#chromosome 1 may need to be added manually (add column at start if it does - there may be a column here with NA in, replace NA with 1)
head *bim
R
mydata <- read.table("data_chr01_sample_filter_subsetted.bim")
head(mydata)
mydata$chr <- 1
head(mydata)
mydata <- mydata[,c(6,1,2,3,4,5)]
head(mydata)
write.table(mydata, "data_chr01_sample_filter_subsetted.bim", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
#####################################################################################


#####################################################################################
## perform GCTA analysis
mkdir ./GWAS/ALSPAC/step6_GCTA/step1_GRM

# first create GRM for each chromosome - calculate the GRM for each autosome 
############################
#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=06:00:00
cd $PBS_O_WORKDIR
module add apps/gcta-1.26.0

cd ./GWAS/ALSPAC/step6_GCTA/step1_GRM

dir1=./GWAS/ALSPAC/step6_GCTA/PLINK/
dir2=./GWAS/ALSPAC/step6_GCTA/step1_GRM/output/

VAR1=data_chr03_sample_filter_subsetted

gcta64 --bfile ${dir1}${VAR1%.bgen} --maf 0.01 --make-grm --out ${dir2}${VAR1%.bgen} --thread-num 10
############################
#create multiple qsub scripts based on master_sub
cat ./GWAS/ALSPAC/step6_GCTA/PLINK/run/VAR1.txt | while read i; do echo ${i}; awk '{ if (NR == 13) print "VAR1='${i}'"; else print $0}'  master_sub > ${i}.sh; done

#
mkdir log
mkdir output

#create run file to run all qsub scripts in one
ls |grep .bgen.sh |sed 's/\ / /g' |awk '{print "qsub "$i" -e log -o log"}' > run

#rename run to run.sh
mv run run.sh

#run run.sh file
sh run.sh

#make run directory and move all .sh files there
mkdir run

mv *sh run/
#clean up directories used here
mv *log ./GWAS/ALSPAC/step6_GCTA/step1_GRM/log


#####################################################################################################
# creat a text file for merging the grms together

nano grm_names.txt
data_chr01_sample_filter_subsetted
data_chr02_sample_filter_subsetted
data_chr03_sample_filter_subsetted
data_chr04_sample_filter_subsetted
data_chr05_sample_filter_subsetted
data_chr06_sample_filter_subsetted
data_chr07_sample_filter_subsetted
data_chr08_sample_filter_subsetted
data_chr09_sample_filter_subsetted
data_chr10_sample_filter_subsetted
data_chr11_sample_filter_subsetted
data_chr12_sample_filter_subsetted
data_chr13_sample_filter_subsetted
data_chr14_sample_filter_subsetted
data_chr15_sample_filter_subsetted
data_chr16_sample_filter_subsetted
data_chr17_sample_filter_subsetted
data_chr18_sample_filter_subsetted
data_chr19_sample_filter_subsetted
data_chr20_sample_filter_subsetted
data_chr21_sample_filter_subsetted
data_chr22_sample_filter_subsetted

############################
#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
cd $PBS_O_WORKDIR
cd ./GWAS/ALSPAC/step6_GCTA/step1_GRM/output/
module add apps/gcta-1.26.0

gcta64 --mgrm ./GWAS/ALSPAC/step6_GCTA/step1_GRM/output/grm_names.txt --make-grm --out ./GWAS/ALSPAC/step6_GCTA/step1_GRM/output/test

#####################################################################################################
#Input	phenotype	data	from	a	plain	text	file,	e.g.	test.phen.	If	the	phenotypic	value	is	coded	
# as	0	or	1,	then	it	will	be	recognized	as	a	case-control	study	(0	for	controls	and	1	for	cases).	
# Missing	value	should	be	represented	by "-9"	or	"NA".

#phenotype file has no header and columns are; family ID, individual ID, phenotype (0/1)

module add languages/R-3.5-ATLAS-gcc-7.1.0
R

data <- read.table("./GWAS/ALSPAC/step1_phenofile/pheno_file_include-list-only.txt", header = F, skip =2)
head(data)
dim(data)
data <- data[,c(1,2,28)]

#remove M from col 1
head(data)
data$V1 <- gsub("[a-zA-Z]", "", data$V1)
head(data)

write.table(data, "./GWAS/ALSPAC/step6_GCTA/step1_GRM/output/GCTA_pheno.file", 
            row.names = FALSE, col.names = F, quote = FALSE, sep = "\t")

data <- read.table("./GWAS/ALSPAC/step6_GCTA/step1_GRM/output/test.grm.id", header = F)
head(data)
data$V1 <- gsub("[a-zA-Z]", "", data$V2)
head(data)

write.table(data, "./GWAS/ALSPAC/step6_GCTA/step1_GRM/output/test.grm.id", 
            row.names = FALSE, col.names = F, quote = FALSE, sep = "\t")


#####################################################################################################
## run GCTA
############################
#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
cd $PBS_O_WORKDIR
cd ./GWAS/ALSPAC/step6_GCTA/step1_GRM/output/
module add apps/gcta-1.26.0

gcta64 --reml --grm ./GWAS/ALSPAC/step6_GCTA/step1_GRM/output/test --pheno ./GWAS/ALSPAC/step6_GCTA/step1_GRM/output/GCTA_pheno.file --prevalence 0.20 --out ./GWAS/ALSPAC/step6_GCTA/step1_GRM/output/GCTA_output/prev2022_GCTA_grm.out





