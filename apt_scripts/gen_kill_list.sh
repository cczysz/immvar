#!/bin/bash

# Note: The output file needs to be filtered for lines with only 3 columns, and the transcript_cluster_ids need to be incremented by 1 to be converted to probeset_ids

MAP=/group/stranger-lab/moliva/ImmVar/probes_mapping/Hugene_info/HuGene-1_0-st-v1.r4.mps
CLF=/group/stranger-lab/moliva/ImmVar/probes_mapping/Hugene_info/HuGene-1_0-st-v1.r4.clf

F=/group/stranger-lab/immvar/apt_scripts/probes

echo -e "probe_id\tprobeset_id\tx\ty" > /group/stranger-lab/immvar/apt_scripts/kill_list.txt.tmp;
while read line
do
	#pset=$(grep -e '\s$line\s' $MAP | cut -f 1);
	#loc=$(grep -e "^$line\s" $CLF | cut -f 2,3);
	entry=$(grep -w -m 1 $line $F)
	echo $entry >> /group/stranger-lab/immvar/apt_scripts/kill_list.txt.tmp;
done < /group/stranger-lab/immvar/apt_scripts/probes_to_kill.txt
