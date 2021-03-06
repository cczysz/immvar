## gather top 20 ###

# for i in $(echo "rep.emtab2232 rep.GSE56045 meta.cd14 meta.cd14.all1.txt immvar.CD14 immvar.CD4 meta.cd4_all1.txt meta.cd4 rep.GenCord rep.GSE56580"); do
#for e in $(echo "c1.all c2.cgp c2.cp.biocarta c2.cp.kegg c2.cp.reactome c3.mir c3.tft c4.cgn c4.cm c5.bp c5.cc c5.mf c6.all c7.all"); do
FDIR=/scratch/t.cczysz/gsea_out/
WDIR=/scratch/t.cczysz/gsea_test/

sets="c1.all c2.cgp c2.cp.biocarta c2.cp.kegg c2.cp.reactome c3.mir c3.tft c4.cgn c4.cm c5.bp c5.cc c5.mf c6.all c7.all hr.all"
#study="rep.emtab2232 rep.GSE56045 meta.cd14 meta.cd14.all1.txt immvar.CD14 immvar.CD4 meta.cd4_all1.txt meta.cd4 rep.GenCord rep.GSE56580"
study="meta.all.cd14_all1.txt meta.all.cd14_immvar1.txt immvar.CD14 immvar.CD4 meta.all.cd4_immvar1.txt meta.all.cd4_all1.txt"

cd $FDIR;
for e in $(echo $sets); do
	echo $e;
	head -n 21 *.$e.*/gsea_report_for_na_pos_*.xls | grep 'tag' | awk '{print $1}'  | sort | uniq > $e.top20_gene_sets.na_pos.txt ;
	for u in $(cat $e.top20_gene_sets.na_pos.txt); do
		echo $u;
		for i in $(echo $study); do
#echo $(grep -cw $u *$e*/gsea_report_for_na_pos_*xls | grep -oP ":[01]" | sed 's/://g' | sed 's/\s/\t/g' ) >> $e.top20_gene_sets.na_pos.values.txt;
			echo $(grep -cw "$u" <(head -21 $i*$e*/gsea_report_for_na_pos_*xls)) >> tmp.txt
		done;
		echo $(cat tmp.txt | tr '\n' '\t') >> $e.top20_gene_sets.na_pos.values.txt;
		rm tmp.txt;
	done;
	paste $e.top20_gene_sets.na_pos.txt $e.top20_gene_sets.na_pos.values.txt > $e.top20_gene_sets.na_pos.matrix.txt;
done;

mkdir top20; mv *.top20* top20/;

# For each collection/study pair, print name of top 20 gene set with p-value < 0.05
for i in $(echo $study); do
	for e in $(echo $sets); do
		awk 'NR>=2&&NR<=21' $i.rnk.$e*/gsea_report_for_na_pos_*xls | awk -F "\t" '{if ($7 < 0.05) {print $1}}' > $i.$e.na_pos.filter.nompval.0.05.txt ;
	done;
done;

# For each collection/study pair, print name of top 20 gene set with q-value < 0.25
for i in $(echo $study); do
	for e in $(echo $sets); do
		awk 'NR>=2&&NR<=21' $i.rnk.$e*/gsea_report_for_na_pos_*xls | awk -F "\t" '{if ($8 < 0.25) {print $1}}' > $i.$e.na_pos.filter.fdr.0.25.txt;
	done;
done;

mkdir significant; mv *txt significant;

echo $PWD 

for e in $(echo $sets); do
	echo $e;
	for u in $(cat top20/$e.top20_gene_sets.na_pos.txt); do
		for i in $(echo $study); do
			echo $(grep -cw "$u" significant/$i.$e.na_pos.filter.fdr.0.25.txt) >> tmp.txt
		done;
		echo $(cat tmp.txt | tr '\n' '\t') >> $e.top20_gene_sets.na_pos.filter.fdr.0.25.txt;
		rm tmp.txt;
	done;
	paste top20/*$e.top20_gene_sets.na_pos.txt $e.top20_gene_sets.na_pos.filter.fdr.0.25.txt > $e.top20_gene_sets.na_pos.filter.fdr.0.25.matrix.txt;
done;

for e in $(echo $sets); do
	echo $e;
	for u in $(cat top20/$e.top20_gene_sets.na_pos.txt); do
		for i in $(echo $study); do
			echo $(grep -cw "$u" significant/$i.$e.na_pos.filter.nompval.0.05.txt) >> tmp.txt
		done;
		echo $(cat tmp.txt | tr '\n' '\t') >> $e.top20_gene_sets.na_pos.filter.nompval.0.05.txt;
		rm tmp.txt;
	done;
	paste top20/*$e.top20_gene_sets.na_pos.txt $e.top20_gene_sets.na_pos.filter.nompval.0.05.txt > $e.top20_gene_sets.na_pos.filter.nompval.0.05.matrix.txt;
done;
cp top20/*matrix.txt .;
