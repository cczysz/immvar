#!/bin/bash

#PBS -N czysz_mesa_mono
#PBS -S /bin/bash
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=32gb

#PBS -o $HOME/mesa_mono.out
#PBS -e $HOME/mesa_mono.err

module load R/3.1.0

R CMD BATCH --no-save --no-restore /group/stranger-lab/immvar/rep_scripts/mesa_monocytes.R
