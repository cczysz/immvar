## gather top 20 ###

FDIR=/scratch/t.cczysz/gsea_out/
WDIR=/scratch/t.cczysz/gsea_test/

cd $FDIR;
for e in $(echo "c1.all c2.cgp c2.cp.biocarta c2.cp.kegg c2.cp.reactome c3.mir c3.tft c4.cgn c4.cm c5.bp c5.cc c5.mf c6.all c7.all"); do
	echo $e; 
	head -n 21 *.$e.*/gsea_report_for_na_pos_*.xls | grep 'tag' | awk '{print $1}'  | sort | uniq > $e.top20_gene_sets.na_pos.txt ;
	for u in $(cat $e.top20_gene_sets.na_pos.txt); do
		echo $u;
		#echo $(grep -w $u *$e*/gsea_report_for_na_pos_*xls | awk '{print $9}' ) >> $e.top20_gene_sets.na_pos.values.txt;
		echo $(grep -cw $u *$e*/gsea_report_for_na_pos_*xls | grep -oP ":[01]" | sed 's/://g' | sed 's/\s/\t/g' ) >> $e.top20_gene_sets.na_pos.values.txt;
	done;
	paste $e.top20_gene_sets.na_pos.txt $e.top20_gene_sets.na_pos.values.txt > $e.top20_gene_sets.na_pos.matrix.txt;
done;

mkdir top20; mv *.top20* top20/;

#for i in $(echo "cd4 cd14 intersect joint"); do
#for i in $(echo "immvar.CD14 immvar.CD4 rep.emtab2232 rep.GenCord rep.GSE56580 rep.GSE56045 meta.cd14 meta.cd4"); do
for i in $(echo "immvar.CD14 immvar.CD4 rep.emtab2232 rep.GSE56045 meta.cd14 rep.GenCord rep.GSE56580 meta.cd4"); do
	for e in $(echo "c1.all c2.cgp c2.cp.biocarta c2.cp.kegg c2.cp.reactome c3.mir c3.tft c4.cgn c4.cm c5.bp c5.cc c5.mf c6.all c7.all"); do
		cat $i.rnk.$e*/gsea_report_for_na_pos_*xls | awk -F "\t" '{if ($7 < 0.05) {print $1}}' > $i.$e.na_pos.filter.nompval.0.05.txt ;
	done;
done;

#for i in $(echo "cd4 cd14 intersect joint"); do
#for i in $(echo "immvar.CD14 immvar.CD4 rep.emtab2232 rep.GenCord rep.GSE56580 rep.GSE56045 meta.cd14 meta.cd4"); do
for i in $(echo "immvar.CD14 immvar.CD4 rep.emtab2232 rep.GSE56045 meta.cd14 rep.GenCord rep.GSE56580 meta.cd4"); do
	for e in $(echo "c1.all c2.cgp c2.cp.biocarta c2.cp.kegg c2.cp.reactome c3.mir c3.tft c4.cgn c4.cm c5.bp c5.cc c5.mf c6.all c7.all"); do
		cat $i.rnk.$e*/gsea_report_for_na_pos_*xls | awk -F "\t" '{if ($8 < 0.25) {print $1}}' > $i.$e.na_pos.filter.fdr.0.25.txt;
	done;
done;

mkdir significant; mv *txt significant;

for e in $(echo "c1.all c2.cgp c2.cp.biocarta c2.cp.kegg c2.cp.reactome c3.mir c3.tft c4.cgn c4.cm c5.bp c5.cc c5.mf c6.all c7.all"); do
	echo $e; 
	for u in $(cat top20/$e.top20_gene_sets.na_pos.txt); do
		echo $(grep -cw $u significant/*$e*.na_pos.filter.fdr.0.25.txt | grep -oP ":[01]" | sed 's/://g' | sed 's/\s/\t/g' ) >> $e.top20_gene_sets.na_pos.filter.fdr.0.25.txt;
	done;
	paste top20/$e.top20_gene_sets.na_pos.txt $e.top20_gene_sets.na_pos.filter.fdr.0.25.txt > $e.top20_gene_sets.na_pos.filter.fdr.0.25.matrix.txt; 
done;

for e in $(echo "c1.all c2.cgp c2.cp.biocarta c2.cp.kegg c2.cp.reactome c3.mir c3.tft c4.cgn c4.cm c5.bp c5.cc c5.mf c6.all c7.all"); do
	echo $e;
	for u in $(cat top20/$e.top20_gene_sets.na_pos.txt); do
		#echo $(grep -cw $u significant/*$e*.na_pos.filter.nompval.0.05.txt | grep -oP ":[01]" | sed 's/://g' | sed 's/\s/\t/g' ) >> $e.top20_gene_sets.na_pos.filter.nompval.0.05.txt;
		echo $(grep -cw $u significant/*$e*.na_pos.filter.nompval.0.05.txt | grep -oP ":[01]" | sed 's/://g' | sed 's/\s/\t/g' ) >> $e.top20_gene_sets.na_pos.filter.nompval.0.05.txt;
	done;
	paste top20/$e.top20_gene_sets.na_pos.txt $e.top20_gene_sets.na_pos.filter.nompval.0.05.txt > $e.top20_gene_sets.na_pos.filter.nompval.0.05.matrix.txt;
done;
