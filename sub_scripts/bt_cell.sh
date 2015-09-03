#!/bin/bash

#PBS -N czysz_de_bt
#PBS -S /bin/bash
#PBS -l walltime=5:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb

#PBS -o $HOME/
#PBS -e $HOME/de_bt.err

module load R/3.1.0

R CMD BATCH --no-save --no-restore /group/stranger-lab/immvar/DE_analysis/bt_celltype.R
#. /group/stranger-lab/immvar/DE_analysis/generateFiltBTCell.R -p 'Caucasian'
