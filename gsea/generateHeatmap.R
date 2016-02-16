library(gplots)

my_palette <- colorRampPalette(c("white","red"))(n = 10);
gene_set_classes=c("c1.all","c2.cgp","c2.cp.biocarta","c2.cp.kegg","c2.cp.reactome","c3.mir","c3.tft","c4.cgn","c4.cm","c5.bp","c5.cc","c5.mf","c6.all","c7.all", "hr.all");
sexes=c("na_pos");
pdf("heatmaps.pdf", width=18, height=13)
#par(mai=c(2,0.5,1,10))
for (gene_set_class in gene_set_classes) {
	for (sex in sexes) {
		file=paste(c(gene_set_class,'top20_gene_sets',sex,'total.matrix.txt'),collapse='.');
		mat=as.matrix(read.table(file,sep=" ",row.names=1,header=T));
		if (gene_set_class%in%c("c2.cp.biocarta","c2.cp.kegg","c2.cp.reactome")) {
			rn <- matrix(rownames(mat), ncol=1)
			new.rn <- apply(rn, 1, function(x) {paste(unlist(strsplit(x, split="_"))[-1], collapse='_')})
			rownames(mat) <- new.rn
		}
		heatmap.2(mat, Rowv = FALSE, Colv = FALSE,cellnote=mat, dendrogram = "none", notecol = "black", notecex = 0.7,trace = "none", key = F, margins = c(10,29),col=my_palette,cexCol = 0.75,cexRow=0.75,lhei=c(1,6));
		par(lend = 1)           
		legend("topleft",      
    			legend = c("1: Present in top 20", "2: Passes Nominal Pvalue < 0.05", "4: Passes FDR Threshold < 25%"),
    			col = c("white", "white", "white"),  
    			lty= 0,             
			cex=0.5,
			title="Heatmap score breakdown"
		)
		title(gene_set_class, cex.main = 1);
	}
}
dev.off()
