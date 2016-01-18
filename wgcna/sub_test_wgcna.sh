#!/bin/bash

#PBS -N czysz_test_wgnca
#PBS -S /bin/bash
#PBS -l walltime=25:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=36gb

#PBS -o $HOME/
#PBS -e $HOME/

module load R/3.1.0
ulimit -s unlimited

Rscript --vanilla /group/stranger-lab/immvar/wgcna/test_wgcna.R
