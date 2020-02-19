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
#PBS -e log/
#PBS -o log/
cd $PBS_O_WORKDIR

module add apps/qctool-2.0

cd ./001_glycosuria/reviewer_analysis/GWAS_revision2/step6

dir=./alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data/dosage_bgen
dir1=./001_glycosuria/reviewer_analysis/GWAS_revision2
dir2=./alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data

VAR1=data_chr01.bgen

qctool -g ${dir}/${VAR1} -s ${dir2}/data.sample -og ${VAR1%.bgen}_sample_filter.bgen -os ${VAR1%.bgen}.sample -incl-samples ${dir1}/pheno_file/offspring_GWAS/include.txt -log ${VAR1%.bgen}_log.txt
############################

#create .txt file of all names for submission scripts based on files in a directory
ls -1 ./alspac/studies/latest/alspac/genetic/variants/arrays/gwas/imputed/1000genomes/released/2015-10-30/data/dosage_bgen > bgen_filenames.txt

#create multiple qsub scripts based on master_sub
cat bgen_filenames.txt | while read i; do echo ${i}; awk '{ if (NR == 17) print "VAR1='${i}'"; else print $0}' master_sub > ${i}.sh; done

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
mv data_chr01.sample data.sample.new
rm *sample
mv data.sample.new output/data.sample
mv *txt log
mv data* output/
mv master_sub run/

