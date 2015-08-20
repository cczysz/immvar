#!/bin/bash

#PBS -N czysz_go
#PBS -S /bin/bash
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb

#PBS -o $HOME/
#PBS -e $HOME/go.err

module load R/3.1.0

R CMD BATCH --no-save --no-restore /group/stranger-lab/immvar/DE_analysis/GO.R
