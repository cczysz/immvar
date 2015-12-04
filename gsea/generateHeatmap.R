library(gplots)

my_palette <- colorRampPalette(c("white","red"))(n = 4);
gene_set_classes=c("c1.all","c2.cgp","c2.cp.biocarta","c2.cp.kegg","c2.cp.reactome","c3.mir","c3.tft","c4.cgn","c4.cm","c5.bp","c5.cc","c5.mf","c6.all","c7.all");
sexes=c("na_pos");
pdf("heatmaps.pdf");
for (gene_set_class in gene_set_classes) {
	for (sex in sexes) {
		file=paste(c(gene_set_class,'top20_gene_sets',sex,'total.matrix.txt'),collapse='.');
		mat=as.matrix(read.table(file,sep=" ",row.names=1,header=T));
		heatmap.2(mat, Rowv = FALSE, Colv = FALSE, dendrogram = "none", notecol = "black", notecex = 0.5,trace = "none", key = F, margins = c(7, 11),col=my_palette,cexCol = 1,cexRow=0.5,lhei=c(1,10));
		title(file, cex.main = 1);
	}
}
dev.off()
