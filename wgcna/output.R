setwd('/group/stranger-lab/czysz/ImmVar')
load('/group/stranger-lab/immvar_rep/mappings/annot.Robj')
library(gplots)

cells <- c('CD14', 'CD4')
sexes <- c('male', 'female')

go.suffix <- "GO.Robj" # CD14_male.GO.Robj; GOenr
col.suffix <- "merge.Robj" # CD14_malemerge.Robj

go.results <- list()
col.results <- list()

for (cell in cells) {
for (sex in sexes) {

	# Load expression data for gene names
	expr.in <- paste('/group/stranger-lab/moliva/ImmVar/Robjects/', paste(cell, '.joint.norm.exp_genes.Robj', sep=''), sep='')
	if (cell == 'CD14') { load(file=expr.in);expr<-exp_genes.cd14.joint.norm } else {load(file=expr.in);expr<-exp_genes.cd4.joint.norm}

	load(file=paste(paste(cell, sex, sep='_'), go.suffix, sep='.'))
	go.results[[cell]][[sex]] <- GOenr

	load(file=paste(paste(cell, sex, sep='_'), col.suffix, sep=''))
	names(moduleColors) <- rownames(expr)
	col.results[[cell]][[sex]] <- moduleColors
}
}

gene.colors <- data.frame(CD14.female=col.results$CD14$female, CD14.male=col.results$CD14$male, 
	CD4.female=col.results$CD4$female, CD4.male=col.results$CD4$male)

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

cd14.sharing.matrix <- cd14.intersect.matrix/cd14.union.matrix
heatmap.2(cd14.sharing.matrix, margins=c(7,7), key=F, Rowv=F, Colv=F, cexRow=0.8, cexCol=0.8, notecol='black', cellnote=trunc(100*cd14.sharing.matrix), trace='none', notecex=0.5,
	ylab='Male', xlab='Female', main='CD14')

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
heatmap.2(cd4.sharing.matrix, margins=c(7,7), key=F, Rowv=F, Colv=F, cexRow=0.8, cexCol=0.8, notecol='black', cellnote=trunc(100*cd4.sharing.matrix), trace='none', notecex=0.5,
	ylab='Male', xlab='Female', main='CD4')
dev.off()
