#!/bin/bash

#PBS -N czysz_de
#PBS -S /bin/bash
#PBS -l walltime=30:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=16gb

#PBS -o $HOME/
#PBS -e $HOME/de.err

module load R/3.1.0

#DATA_DIR=/home/t.cczysz/mRNA_expression/CEL_files/CD4
R CMD BATCH --no-save --no-restore /group/stranger-lab/immvar/DE_analysis/de_analysis.R
