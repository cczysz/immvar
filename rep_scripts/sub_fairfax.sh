#!/bin/bash

#PBS -N czysz_fairfax
#PBS -S /bin/bash
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=32gb

#PBS -o $HOME/fairfax.out
#PBS -e $HOME/fairfax.err

module load R/3.1.0

#DATA_DIR=/home/t.cczysz/mRNA_expression/CEL_files/CD4
R CMD BATCH --no-save --no-restore /group/stranger-lab/immvar_rep/emtab2232/fairfax.R
