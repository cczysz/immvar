#!/bin/bash

#PBS -N czysz_gsea
#PBS -S /bin/bash
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb

#PBS -o $HOME/
#PBS -e $HOME/de.err

module load java

bash /group/stranger-lab/immvar/gsea/runGSEA.preRank.bash