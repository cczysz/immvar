#!/bin/bash

bash generate.matricies.bash
rm ./matricies/*.txt
cp /scratch/t.cczysz/gsea_out/*matrix.txt ./matricies/

cd ./matricies/
perl ../generate_final_matrix.pl
R CMD BATCH --no-restore --no-save ../generateHeatmap.R
