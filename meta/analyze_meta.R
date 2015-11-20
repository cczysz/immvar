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

bf.cd4 <- 0.05/nrow(cd4.meta)
bf.cd14 <- 0.05/nrow(cd14.meta)

qval.cd4 <- p.adjust(cd4.meta$P.value, method="BH")
cd4.meta <- cbind(cd4.meta, Q.value = qval.cd4)
cd4.meta <- cd4.meta[order(cd4.meta$Q.value),]
qval.cd14 <- p.adjust(cd14.meta$P.value, method="BH")
cd14.meta <- cbind(cd14.meta, Q.value = qval.cd14)
cd14.meta <- cd14.meta[order(cd14.meta$Q.value),]

meta_sig_cd4_bon <- cd4.meta[cd4.meta$P.value<=bf.cd4, ]
meta_sig_cd14_bon <- cd14.meta[cd14.meta$P.value<=bf.cd14, ]

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

	if (cell=='CD14') { load(paste('fit',pop, cell, 'Robj', sep='.'));cd14.fits[[pop]]<-eb.fit }
	else { load(paste('fit',pop, cell, 'Robj', sep='.'));cd4.fits[[pop]]<-eb.fit }
}
}

sig_cd14_bh<-list()
sig_cd14_bon<-list()

sig_cd4_bh<-list()
sig_cd4_bon<-list()

for (i in seq(3)) {
	sig_cd14_bh[[i]]<-topTable(cd14.fits[[i]], number=Inf, p.value=0.05)
	sig_cd14_bon[[i]]<-topTable(cd14.fits[[i]], number=Inf, adjust.method="bonferroni", p.value=0.05)
	
	sig_cd4_bh[[i]]<-topTable(cd4.fits[[i]], number=Inf, p.value=0.05)
	sig_cd4_bon[[i]]<-topTable(cd4.fits[[i]], number=Inf, adjust.method="bonferroni", p.value=0.05)
}

# CD14 Venn Diagram

venn.diagram(list(cau=rownames(sig_cd14_bh[[1]]), afr=rownames(sig_cd14_bh[[2]]), asn=rownames(sig_cd14_bh[[3]])),
	 filename='/group/stranger-lab/czysz/cd14_separate_venn.tiff',
	 fontfamily="Helvetica",
	 main.fontfamily="Helvetica",
	 cat.fontfamily="Helvetica",
	 sub.fontfamily="Helvetica",
	 fill=topo.colors(3),
	 main="CD14 - Separate VennDiagram", sub.cex=1.5,
	 width=10, height=10, units="in")

# CD4
venn.diagram(list(cau=rownames(sig_cd4_bh[[1]]), afr=rownames(sig_cd4_bh[[2]]), asn=rownames(sig_cd4_bh[[3]])),
	 filename='/group/stranger-lab/czysz/cd4_separate_venn.tiff',
	 fontfamily="Helvetica",
	 main.fontfamily="Helvetica",
	 cat.fontfamily="Helvetica",
	 sub.fontfamily="Helvetica",
	 fill=topo.colors(3),
	 main="CD4 - Separate VennDiagram", sub.cex=1.5,
	 width=10, height=10, units="in")
load('/group/stranger-lab/immvar_data/fit.joint.CD14.Robj')
joint.cd14.fit <- eb.fit
joint.cd14.sig <- topTable(joint.cd14.fit, number=Inf, p.value=0.05)
joint.cd14.fit$q.value <- p.adjust(joint.cd14.fit$p.value, method="fdr")

load('/group/stranger-lab/immvar_data/fit.joint.CD4.Robj')
joint.cd4.fit <- eb.fit
joint.cd4.fit$q.value <- p.adjust(joint.cd4.fit$p.value, method="fdr")
joint.cd4.sig <- topTable(joint.cd4.fit, number=Inf, p.value=0.05)

venn.diagram(list(meta=meta_sig_cd4_bh$MarkerName, joint=rownames(joint.cd4.sig)),
	 filename='/group/stranger-lab/czysz/cd4_venn.tiff',
	 fontfamily="Helvetica",
	 main.fontfamily="Helvetica",
	 cat.fontfamily="Helvetica",
	 sub.fontfamily="Helvetica",
	 fill=topo.colors(2),
	 main="CD4 VennDiagram", sub.cex=1.5,
	 width=10, height=10, units="in")

venn.diagram(list(meta=meta_sig_cd14_bh$MarkerName, joint=rownames(joint.cd14.sig)),
	 filename='/group/stranger-lab/czysz/cd14_venn.tiff', main.fontfamily="Helvetica", sub.fontfamily="Helvetica",
	 fontfamily="Helvetica",
	 cat.fontfamily="Helvetica",
	 fill=topo.colors(2),
	 main="CD14 VennDiagram", sub.cex=1.5,
	 width=10, height=10, units="in")
