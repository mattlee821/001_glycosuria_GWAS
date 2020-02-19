##########################
## step3
## Estimate SNP Summary Stats
##########################

## Create a new SNPSTATS file on your newly subsetted data set from step2

#create master_sub qsub script as template
###########################
#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -e log/
#PBS -o log/
cd $PBS_O_WORKDIR

module add apps/qctool-2.0

cd ./001_glycosuria/reviewer_analysis/GWAS_revision2/step3

dir=./001_glycosuria/reviewer_analysis/GWAS_revision2/step2/output

VAR1=data_chr01_sample_filter.bgen

qctool -g ${dir}/${VAR1} -snp-stats -osnp ${VAR1%.bgen}_snpstats.txt
############################

#create .txt file of all names for submission scripts based on files in a directory
ls -1 ./001_glycosuria/reviewer_analysis/GWAS_revision2/step2/output > sample_filter_filenames.txt
# remove last line of this new .txt file - data.sample

#create multiple qsub scripts based on master_sub
cat sample_filter_filenames.txt | while read i; do echo ${i}; awk '{ if (NR == 15) print "VAR1='${i}'"; else print $0}' master_sub > ${i}.sh; done

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
mv *.sh run
mv sample_filter_filenames.txt run
mv master_sub run
############################
# wait till the jobs have finished then clean up the directory

#clean up
mv data* output/
#####################################################################################
