setwd('/group/stranger-lab/czysz/ImmVar')
load('/group/stranger-lab/immvar_rep/mappings/annot.Robj')
library(gplots)

cells <- c('CD14', 'CD4')
sexes <- c('male', 'female')

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

gene.colors <- data.frame(CD14.female=mod.results$CD14$female, CD14.male=mod.results$CD14$male, 
	CD4.female=mod.results$CD4$female, CD4.male=mod.results$CD4$male)

cd14.intersect.matrix <- matrix(0, nrow=length(unique(gene.colors[,2])),ncol=length(unique(gene.colors[,1])))
cd14.union.matrix <- matrix(0, nrow=length(unique(gene.colors[,2])),ncol=length(unique(gene.colors[,1])))
rownames(cd14.intersect.matrix) <- names(table(gene.colors$CD14.male))
rownames(cd14.union.matrix) <- names(table(gene.colors$CD14.male))
colnames(cd14.intersect.matrix) <- names(table(gene.colors$CD14.female))
colnames(cd14.union.matrix) <- names(table(gene.colors$CD14.female))

for (i in seq(length(rownames(cd14.intersect.matrix)))) {
for (j in seq(length(colnames(cd14.intersect.matrix)))) {

	mal.col <- rownames(cd14.intersect.matrix)[i]
	fem.col <- colnames(cd14.intersect.matrix)[j]

	mal.genes <- rownames(subset(gene.colors, CD14.male==mal.col))
	fem.genes <- rownames(subset(gene.colors, CD14.female==fem.col))

	cd14.intersect.matrix[i,j] <- length(intersect(mal.genes, fem.genes))
	cd14.union.matrix[i,j] <- length(union(mal.genes, fem.genes))

}
}

cd14.femnames = paste(names(table(gene.colors$CD14.female)), table(gene.colors$CD14.female), sep=" N=")
cd14.malnames = paste(names(table(gene.colors$CD14.male)), table(gene.colors$CD14.male), sep=" N=")
cd4.femnames = paste(names(table(gene.colors$CD4.female)), table(gene.colors$CD4.female), sep=" N=")
cd4.malnames = paste(names(table(gene.colors$CD4.male)), table(gene.colors$CD4.male), sep=" N=")

cd14.sharing.matrix <- cd14.intersect.matrix/cd14.union.matrix
pdf('/group/stranger-lab/czysz/ImmVar/plots/wgnca_sharing.pdf', width=8, height=8)
heatmap.2(cd14.sharing.matrix, margins=c(7,7), key=T, density.info='none', key.title='', key.xlab='Overlap', Rowv=F, Colv=F, cexRow=0.8, cexCol=0.8, notecol='black', cellnote=trunc(100*cd14.sharing.matrix), trace='none', notecex=0.5,
	ylab='Male', xlab='Female', main='CD14', col=heat.colors(100, alpha=0.5)[seq(100,1)], labRow=cd14.malnames, labCol=cd14.femnames)

cd4.intersect.matrix <- matrix(0, nrow=length(unique(gene.colors[,4])),ncol=length(unique(gene.colors[,3])))
cd4.union.matrix <- matrix(0, nrow=length(unique(gene.colors[,4])),ncol=length(unique(gene.colors[,3])))
rownames(cd4.intersect.matrix) <- names(table(gene.colors$CD4.male))
rownames(cd4.union.matrix) <- names(table(gene.colors$CD4.male))
colnames(cd4.intersect.matrix) <- names(table(gene.colors$CD4.female))
colnames(cd4.union.matrix) <- names(table(gene.colors$CD4.female))

for (i in seq(length(rownames(cd4.intersect.matrix)))) {
for (j in seq(length(colnames(cd4.intersect.matrix)))) {

	mal.col <- rownames(cd4.intersect.matrix)[i]
	fem.col <- colnames(cd4.intersect.matrix)[j]

	mal.genes <- rownames(subset(gene.colors, CD4.male==mal.col))
	fem.genes <- rownames(subset(gene.colors, CD4.female==fem.col))

	cd4.intersect.matrix[i,j] <- length(intersect(mal.genes, fem.genes))
	cd4.union.matrix[i,j] <- length(union(mal.genes, fem.genes))

}
}

cd4.sharing.matrix <- cd4.intersect.matrix/cd4.union.matrix
heatmap.2(cd4.sharing.matrix, margins=c(7,7), key=T, density.info='none', key.title='', key.xlab='Overlap', Rowv=F, Colv=F, cexRow=0.8, cexCol=0.8, notecol='black', cellnote=trunc(100*cd4.sharing.matrix), trace='none', notecex=0.5,
	ylab='Male', xlab='Female', main='CD4', col=heat.colors(100, alpha=0.5)[seq(100,1)], labRow=cd4.malnames, labCol=cd4.femnames)

dev.off()
