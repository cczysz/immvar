setwd('/group/stranger-lab/czysz/ImmVar')
load('/group/stranger-lab/immvar_rep/mappings/annot.Robj')
library(gplots)

cells <- c('CD14', 'CD4')
sexes <- c('female1_test', 'female2_test')

go.suffix <- "GO.Robj" # CD14_male.GO.Robj; GOenr
col.suffix <- "merge.Robj" # CD14_malemerge.Robj

go.results <- list()
mod.results <- list()
col.results <- list()

for (cell in cells) {
for (sex in sexes) {

	# Load expression data for gene names
	expr.in <- paste('/group/stranger-lab/moliva/ImmVar/Robjects/', paste(cell, '.joint.norm.exp_genes.Robj', sep=''), sep='')
	if (cell == 'CD14') { load(file=expr.in);expr<-exp_genes.cd14.joint.norm } else {load(file=expr.in);expr<-exp_genes.cd4.joint.norm}

	load(file=paste(paste(cell, sex, sep='_'), go.suffix, sep='.'))
	go.results[[cell]][[sex]] <- GOenr

	load(file=paste(paste(cell, sex, sep='_'), col.suffix, sep=''))
	names(moduleLabels) <- rownames(expr)
	names(moduleColors) <- rownames(expr)
	mod.results[[cell]][[sex]] <- moduleLabels
	col.results[[cell]][[sex]] <- moduleColors
}
}

gene.colors <- data.frame(CD14.female1=mod.results$CD14$female1_test, CD14.female2=mod.results$CD14$female2_test, 
	CD4.female1=mod.results$CD4$female1_test, CD4.female2=mod.results$CD4$female2_test)

cd14.intersect.matrix <- matrix(0, nrow=length(unique(gene.colors[,2])),ncol=length(unique(gene.colors[,1])))
cd14.union.matrix <- matrix(0, nrow=length(unique(gene.colors[,2])),ncol=length(unique(gene.colors[,1])))
rownames(cd14.intersect.matrix) <- names(table(gene.colors$CD14.female2))
rownames(cd14.union.matrix) <- names(table(gene.colors$CD14.female2))
colnames(cd14.intersect.matrix) <- names(table(gene.colors$CD14.female1))
colnames(cd14.union.matrix) <- names(table(gene.colors$CD14.female1))

for (i in seq(nrow(cd14.intersect.matrix))) {
for (j in seq(ncol(cd14.intersect.matrix))) {

	mal.col <- rownames(cd14.intersect.matrix)[i]
	fem.col <- colnames(cd14.intersect.matrix)[j]

	mal.genes <- rownames(subset(gene.colors, CD14.female2==mal.col))
	fem.genes <- rownames(subset(gene.colors, CD14.female1==fem.col))

	cd14.intersect.matrix[i,j] <- length(intersect(mal.genes, fem.genes))
	cd14.union.matrix[i,j] <- length(union(mal.genes, fem.genes))

}
}

cd14.femnames = paste(names(table(gene.colors$CD14.female1)), table(gene.colors$CD14.female1), sep=" N=")
cd14.fem1.order <- order(table(gene.colors$CD14.female1), decreasing=T)
cd14.malnames = paste(names(table(gene.colors$CD14.female2)), table(gene.colors$CD14.female2), sep=" N=")
cd14.fem2.order <- order(table(gene.colors$CD14.female2), decreasing=T)
cd4.femnames = paste(names(table(gene.colors$CD4.female1)), table(gene.colors$CD4.female1), sep=" N=")
cd4.fem1.order <- order(table(gene.colors$CD4.female1), decreasing=T)
cd4.malnames = paste(names(table(gene.colors$CD4.female2)), table(gene.colors$CD4.female2), sep=" N=")
cd4.fem2.order <- order(table(gene.colors$CD4.female2), decreasing=T)

cd14.sharing.matrix <- cd14.intersect.matrix/cd14.union.matrix
cd14.sharing.matrix <- cd14.sharing.matrix[cd14.fem2.order, cd14.fem1.order]
pdf('/group/stranger-lab/czysz/ImmVar/plots/test_wgcna_sharing.pdf', width=8, height=8)
heatmap.2(cd14.sharing.matrix, margins=c(7,7), key=T, density.info='none', key.title='', key.xlab='Overlap', Rowv=F, Colv=F, cexRow=0.8, cexCol=0.8, notecol='black', cellnote=trunc(100*cd14.sharing.matrix), trace='none', notecex=1,
	ylab='Female2', xlab='Female1', main='CD14', col=heat.colors(100, alpha=0.5)[seq(100,1)], labRow=cd14.malnames[cd14.fem2.order], labCol=cd14.femnames[cd14.fem1.order])

cd4.intersect.matrix <- matrix(0, nrow=length(unique(gene.colors[,4])),ncol=length(unique(gene.colors[,3])))
cd4.union.matrix <- matrix(0, nrow=length(unique(gene.colors[,4])),ncol=length(unique(gene.colors[,3])))
rownames(cd4.intersect.matrix) <- names(table(gene.colors$CD4.female2))
rownames(cd4.union.matrix) <- names(table(gene.colors$CD4.female2))
colnames(cd4.intersect.matrix) <- names(table(gene.colors$CD4.female1))
colnames(cd4.union.matrix) <- names(table(gene.colors$CD4.female1))

for (i in seq(length(rownames(cd4.intersect.matrix)))) {
for (j in seq(length(colnames(cd4.intersect.matrix)))) {

	mal.col <- rownames(cd4.intersect.matrix)[i]
	fem.col <- colnames(cd4.intersect.matrix)[j]

	mal.genes <- rownames(subset(gene.colors, CD4.female2==mal.col))
	fem.genes <- rownames(subset(gene.colors, CD4.female1==fem.col))

	cd4.intersect.matrix[i,j] <- length(intersect(mal.genes, fem.genes))
	cd4.union.matrix[i,j] <- length(union(mal.genes, fem.genes))

}
}

cd4.sharing.matrix <- cd4.intersect.matrix/cd4.union.matrix
cd4.sharing.matrix <- cd14.sharing.matrix[cd4.fem2.order, cd4.fem1.order]
heatmap.2(cd4.sharing.matrix, margins=c(7,7), key=T, density.info='none', key.title='', key.xlab='Overlap', Rowv=F, Colv=F, cexRow=0.8, cexCol=0.8, notecol='black', cellnote=trunc(100*cd4.sharing.matrix), trace='none', notecex=1,
	ylab='Female2', xlab='Female1', main='CD4', col=heat.colors(100, alpha=0.5)[seq(100,1)], labRow=cd4.malnames[cd4.fem2.order], labCol=cd4.femnames[cd4.fem1.order])
#heatmap.2(cd4.sharing.matrix, margins=c(7,7), key=T, density.info='none', key.title='', key.xlab='Overlap', Rowv=F, Colv=F, cexRow=0.8, cexCol=0.8, notecol='black', cellnote=trunc(100*cd4.sharing.matrix), trace='none', notecex=1,
	#ylab='Female2', xlab='Female1', main='CD4', col=heat.colors(100, alpha=0.5)[seq(100,1)], labRow=cd4.fem2.order, labCol=cd4.fem1.order)

dev.off()
