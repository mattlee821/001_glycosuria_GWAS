##########################
## step2
## Subsample Geno data
##########################
## Filter you genetic data set for the individuals with your phenotype and no consent withdrawn.

#create master_sub qsub script as a template
###########################
#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
cd $PBS_O_WORKDIR

module add apps/qctool-2.0

cd ./GWAS/ALSPAC/step2_subsample_genetic_data

dir=./alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data/dosage_bgen
dir2=./alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data

VAR1=data_chr01.bgen

qctool -g ${dir}/${VAR1} -s ${dir2}/data.sample -og ${VAR1%.bgen}_sample_filter.bgen -os ${VAR1%.bgen}.sample -incl-samples include.txt -log ${VAR1%.bgen}_log.txt
############################

cp ./GWAS/ALSPAC/step1_phenofile/include.txt ./GWAS/ALSPAC/step2_subsample_genetic_data

#create .txt file of all names for submission scripts based on files in a directory
ls -1 ./alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data/dosage_bgen > bgen_filenames.txt

#create multiple qsub scripts based on master_sub
cat bgen_filenames.txt | while read i; do echo ${i}; awk '{ if (NR == 14) print "VAR1='${i}'"; else print $0}' master_sub > ${i}.sh; done

#
mkdir log
mkdir output 

#create run file to run all qsub scripts in one
ls |grep .bgen.sh |sed 's/\ / /g' |awk '{print "qsub "$i" -e log -o log"}' > run

#rename run to run.sh
mv run run.sh

#run run.sh file
sh run.sh

#
mkdir run
mv *.sh run
mv bgen_filenames.txt run

############################
# wait till the jobs have finished then clean up the directory

#clean up
cd ./GWAS/ALSPAC/step2_subsample_genetic_data/
mv data_chr01.sample data.sample.new
rm *sample
mv data.sample.new output/data.sample
mv *txt ./GWAS/ALSPAC/step2_subsample_genetic_data/log
mv data* output/


