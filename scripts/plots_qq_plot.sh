###########################
#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
cd $PBS_O_WORKDIR

cd ./GWAS/ALSPAC/plots

module add languages/R-3.5.1-ATLAS-gcc-6.1

Rscript ./GWAS/ALSPAC/plots/qq_plot.R
