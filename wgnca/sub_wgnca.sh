#!/bin/bash

#PBS -N czysz_wgnca
#PBS -S /bin/bash
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=20gb

#PBS -o $HOME/
#PBS -e $HOME/

module load R/3.1.0

R CMD BATCH --no-save --no-restore /group/stranger-lab/immvar/wgnca/wgnca.R
