##########################
## step3
## Estimate SNP Summary Stats
##########################

## Create a new SNPSTATS file on your newly subsetted data set from step2

#create master_sub qsub script as template
############################
#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
cd $PBS_O_WORKDIR

module add apps/qctool-2.0

cd ./GWAS/ALSPAC/step3_create_SNP_stats

dir=./GWAS/ALSPAC/step2_subsample_genetic_data/output                                

VAR1=data_chr01_sample_filter.bgen

qctool -g ${dir}/${VAR1} -snp-stats -osnp ${VAR1%.bgen}_snpstats.txt
############################

#create .txt file of all names for submission scripts based on files in a directory
ls -1 ./GWAS/ALSPAC/step2_subsample_genetic_data/output > sample_filter_filenames.txt
# remove last line of this new .txt file

#create multiple qsub scripts based on master_sub
cat sample_filter_filenames.txt | while read i; do echo ${i}; awk '{ if (NR == 13) print "VAR1='${i}'"; else print $0}' master_sub > ${i}.sh; done

#
mkdir log
mkdir output

#create run file to run all qsub scripts in one
ls |grep .bgen.sh |sed 's/\ / /g' |awk '{print "qsub "$i" -e log -o log"}' > run

#rename run to run.sh
mv run run.sh

#rucdn run.sh file
sh run.sh

#make run directory and move all .sh files there
mkdir run
mv *.sh run
mv sample_filter_filenames.txt run

############################
# wait till the jobs have finished then clean up the directory

#clean up
mv data* output/
#####################################################################################
