library(limma)
library(stats)
library(VennDiagram)
library(gtools)
library(ggplot2)

setwd('/group/stranger-lab/immvar/meta')
load('/group/stranger-lab/immvar_rep/mappings/annot.Robj')

cd4.meta <- read.table(file='cd4_meta1.txt', header=T)
cd14.meta <- read.table(file='cd14_meta1.txt', header=T)

cd4_gene_counts <- table(annots[as.character(cd4.meta$MarkerName),"chr"])
cd14_gene_counts <- table(annots[as.character(cd14.meta$MarkerName),"chr"])

cd4.meta$Q.value <- p.adjust(cd4.meta$P.value, method="BH")
cd4.meta <- cd4.meta[order(cd4.meta$Q.value),]
cd14.meta$Q.value <- p.adjust(cd14.meta$P.value, method="BH")
cd14.meta <- cd14.meta[order(cd14.meta$Q.value),]


meta_sig_cd4_bh <- cd4.meta[cd4.meta$Q.value<=0.05, ]
meta_sig_cd4_bh <- meta_sig_cd4_bh[order(meta_sig_cd4_bh$P.value),]
rownames(meta_sig_cd4_bh) <- as.character(meta_sig_cd4_bh$MarkerName)

cd4_meta_percent <- table(annots[rownames(meta_sig_cd4_bh), "chr"]) / cd4_gene_counts
cd4_gene_pers <- data.frame(chr_name=names(cd4_meta_percent), per = as.numeric(cd4_meta_percent))

g <- ggplot(cd4_gene_pers, aes(x=chr_name, y=per))
pdf(file='/group/stranger-lab/czysz/cd4.meta.locs.pdf')
g + geom_bar(stat="identity")
dev.off()

meta_sig_cd14_bh <- cd14.meta[cd14.meta$Q.value<=0.05, ]
meta_sig_cd14_bh <- meta_sig_cd14_bh[order(meta_sig_cd14_bh$P.value),]
rownames(meta_sig_cd14_bh) <- as.character(meta_sig_cd14_bh$MarkerName)
cd14_meta_percent <- table(annots[rownames(meta_sig_cd14_bh), "chr"]) / cd14_gene_counts
cd14_gene_pers <- data.frame(chr=mixedorder(names(cd14_meta_percent)), chr_name = names(cd14_meta_percent), per = as.numeric(cd14_meta_percent))

setwd('/group/stranger-lab/immvar_data')

cd14.fits<-list()
cd4.fits<-list()
for ( pop in c('Caucasian', 'African-American', 'Asian')) {
for ( cell in c('CD14', 'CD4') ) {
	if (cell=='CD14') { load(paste('fit',pop, cell, 'Robj', sep='.'))
		eb.fit$Q.value <- p.adjust(eb.fit$p.value, method='fdr')
		cd14.fits[[pop]]<-data.frame(eb.fit)
	} else { load(paste('fit',pop, cell, 'Robj', sep='.'))
		eb.fit$Q.value <- p.adjust(eb.fit$p.value, method='fdr')
		cd4.fits[[pop]]<-data.frame(eb.fit) }
}
}

# CD14 Venn Diagram
cd14.venn <- list(cau=rownames(subset(cd14.fits[[1]], Q.value<0.05)),
	afr=rownames(subset(cd14.fits[[2]], Q.value<0.05)),
	asn=rownames(subset(cd14.fits[[3]], Q.value<0.05)))

cd14.shared <- intersect(cd14.venn[[1]], intersect(cd14.venn[[2]], cd14.venn[[3]]))
print(table(annots[cd14.shared, "chr"]))
venn.diagram(cd14.venn,
	 filename='/group/stranger-lab/czysz/cd14_separate_venn.tiff',
	 fontfamily="Helvetica",
	 main.fontfamily="Helvetica",
	 cat.fontfamily="Helvetica",
	 sub.fontfamily="Helvetica",
	 fill=topo.colors(3),
	 main="CD14 - Separate VennDiagram", sub.cex=1.5,
	 width=10, height=10, units="in")

# CD4
cd4.venn <- list(cau=rownames(subset(cd4.fits[[1]], Q.value<0.05)),
	afr=rownames(subset(cd4.fits[[2]], Q.value<0.05)),
	asn=rownames(subset(cd4.fits[[3]], Q.value<0.05)))
cd4.shared <- intersect(cd4.venn[[1]], intersect(cd4.venn[[2]], cd4.venn[[3]]))
print(table(annots[cd4.shared, "chr"]))

venn.diagram(cd4.venn,
	 filename='/group/stranger-lab/czysz/cd4_separate_venn.tiff',
	 fontfamily="Helvetica",
	 main.fontfamily="Helvetica",
	 cat.fontfamily="Helvetica",
	 sub.fontfamily="Helvetica",
	 fill=topo.colors(3),
	 main="CD4 - Separate VennDiagram", sub.cex=1.5,
	 width=10, height=10, units="in")
