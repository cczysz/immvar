load('/home/t.cri.cczysz/mendeliangenes.Robj')
load('/group/stranger-lab/immvar_rep/mappings/annot.Robj')
mendel.genes <- genes

mendel <- read.csv(file='/home/t.cri.cczysz/mmc3.csv', header=T)
mendel <- mendel[,-2]

mendel.disorder <- list()
for (i in seq(nrow(mendel))) {
	genes <- unlist(strsplit(mendel[i,2], split=';'))
	mendel.disorder[[mendel[i,1]]] = genes
}
cd4.rep.meta <- read.table('/group/stranger-lab/immvar/meta/cd4_meta_immvar1.txt', header=T)
cd4.rep.meta <- cbind(cd4.rep.meta,
        q.value=p.adjust(cd4.rep.meta$P.value, method='fdr'),
        rank=rank(cd4.rep.meta$P.value) / nrow(cd4.rep.meta),
        chr=annots[cd4.rep.meta$MarkerName, 'chr'],
        gene=annots[cd4.rep.meta$MarkerName, "symbol_id"])
cd4.rep.meta$chr <- as.character(cd4.rep.meta$chr)
#cd4.rep.meta$xesc <- cd4.rep.meta$gene%in%xesc.genes
cd4.rep.sig.5 <- subset(cd4.rep.meta, q.value<0.05)

cd4.sig.mendel <- cd4.rep.sig.5[cd4.rep.sig.5$gene%in%mendel.genes,]
cd4.sig.mendel <- cd4.sig.mendel[order(cd4.sig.mendel$q.value),]
cd4.disorders <- c()
for (gene in cd4.sig.mendel$gene) {
	for (disease in names(which(unlist(lapply(mendel.disorder, '%in%', gene))))) {
		cd4.disorders <- rbind(cd4.disorders, c(gene, disease, cd4.rep.sig.5[cd4.rep.sig.5$gene==gene,c('chr', 'Zscore', 'P.value', 'q.value')]))
}
}
cd4.disorders <- matrix(unlist(cd4.disorders), ncol=6)
cd4.disorders[,6] <- signif(as.numeric(cd4.disorders[,6]),2)
cd4.disorders[,5] <- signif(as.numeric(cd4.disorders[,5]),2)
write.table(cd4.disorders, file='/group/stranger-lab/czysz/cd4.mendel.csv', quote=F, row.names=F, sep=',')

cd14.rep.meta <- read.table('/group/stranger-lab/immvar/meta/cd14_meta_immvar1.txt', header=T)
cd14.rep.meta <- cbind(cd14.rep.meta,
        q.value=p.adjust(cd14.rep.meta$P.value, method='fdr'),
        rank=rank(cd14.rep.meta$P.value) / nrow(cd14.rep.meta),
        chr=annots[cd14.rep.meta$MarkerName, 'chr'],
        gene=annots[cd14.rep.meta$MarkerName, "symbol_id"])
cd14.rep.meta$chr <- as.character(cd14.rep.meta$chr)
#cd14.rep.meta$xesc <- cd14.rep.meta$gene%in%xesc.genes
cd14.rep.sig.5 <- subset(cd14.rep.meta, q.value<0.05)

cd14.sig.mendel <- cd14.rep.sig.5[cd14.rep.sig.5$gene%in%mendel.genes,]
cd14.sig.mendel <- cd14.sig.mendel[order(cd14.sig.mendel$q.value),]
cd14.disorders <- c()
for (gene in cd14.sig.mendel$gene) {
	for (disease in names(which(unlist(lapply(mendel.disorder, '%in%', gene))))) {
		cd14.disorders <- rbind(cd14.disorders, c(gene, disease, cd14.rep.sig.5[cd14.rep.sig.5$gene==gene,c('chr', 'Zscore', 'P.value', 'q.value')]))
		#print(gene)
		#print(disease) 
}
}

cd14.disorders <- matrix(unlist(cd14.disorders), ncol=6)
cd14.disorders[,6] <- signif(as.numeric(cd14.disorders[,6]),2)
cd14.disorders[,5] <- signif(as.numeric(cd14.disorders[,5]),2)
write.table(cd14.disorders, file='/group/stranger-lab/czysz/cd14.mendel.csv', quote=F, row.names=F, sep=',')
