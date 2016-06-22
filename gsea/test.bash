for f in /home/t.cri.cczysz/redo/*.rnk; do
for set in $(ls /scratch/t.cczysz/gene_symbols/hr* | grep -v c[2-5].all); do
geneset_short=$(basename $set | grep -oP "[ch].*?\.v5\.0");

fname=$(basename $f);
sname=$(basename $set);
echo "$fname";
echo "$sname";
echo "$geneset_short";

java -Xms10g -Xmx15g -cp /home/t.cri.cczysz/gsea2-2.2.0.jar xtools.gsea.GseaPreranked -gmx $set -rnk $f -rpt_label $fname.$geneset_short -out /scratch/t.cczysz/gsea_out/$fname_out/ -collapse false -mode Max_probe -norm meandiv -nperm 1000 -scoring_scheme weighted -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 800 -set_min 15 -zip_report false -gui false &> /scratch/t.cczysz/gsea_out/$fname.$sname.preRank.log;
done
done
