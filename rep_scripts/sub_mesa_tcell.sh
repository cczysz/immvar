#!/bin/bash

#PBS -N czysz_mesa_tcell
#PBS -S /bin/bash
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=32gb

#PBS -o $HOME/mesa_tcell.out
#PBS -e $HOME/mesa_tcell.err

module load R/3.1.0

R CMD BATCH --no-save --no-restore /group/stranger-lab/immvar/rep_scripts/mesa_tcells.R
