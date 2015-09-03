#PBS -l mem=8gb
#PBS -l nodes=1:ppn=1
#PBS -N bt_cell_filter

#PBS -o /home/t.cri.cczysz/bt_cell_filt.out
#PBS -e /home/t.cri.cczysz/bt_cell_filt.err

module load python/3.4.3
module load R/3.1.0
#/apps/compilers/R/3.1.0/bin/Rscript
# /home/t.cri.cczysz/bin/Rscript --no-save --no-restore /group/stranger-lab/immvar/DE_analysis/generateFiltBTCell.R -p 'Caucasian'
/apps/compilers/R/3.1.0/bin/Rscript --no-save --no-restore /group/stranger-lab/immvar/DE_analysis/generateFiltBTCell.R -p 'Caucasian'
