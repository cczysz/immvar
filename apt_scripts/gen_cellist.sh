#!/bin/bash

DIR=/group/stranger-lab/nicolel/mRNA_expression/CEL_files/CD4

echo "cel_files" > celfiles.txt;
while read line
do
	echo $DIR/$line >> celfiles.txt;
done < cd4_cau_celfiles.txt
