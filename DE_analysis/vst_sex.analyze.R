load('/group/stranger-lab/czysz/ImmVar/vst.sex.cd14.Robj')
load('/group/stranger-lab/czysz/ImmVar/vst.sex.cd4.Robj')

annot <- read.table(file='/group/stranger-lab/moliva/ImmVar/probes_mapping/annotations/gencode.v22.TSS.txt', header=T, as.is=T)

top.cd14 <- list(cd14.vst.list[[1]][cd14.vst.list[[1]]>0.2],
	cd14.vst.list[[2]][cd14.vst.list[[2]]>0.2],
	cd14.vst.list[[3]][cd14.vst.list[[3]]>0.2])

top.cd4 <- list(cd4.vst.list[[1]][cd4.vst.list[[1]]>0.2],
	cd4.vst.list[[2]][cd4.vst.list[[2]]>0.2],
	cd4.vst.list[[3]][cd4.vst.list[[3]]>0.2])

all.cd14.genes <- names(unlist(top.cd14))
all.cd4.genes <- names(unlist(top.cd4))

cd14.multi <- table(all.cd14.genes)[table(all.cd14.genes)>1]
cd4.multi <- table(all.cd4.genes)[table(all.cd4.genes)>1]

top.cd14.shared <- c()
for (gene in names(cd14.multi)) {
	symbol <- annot[annot$ensembl_id==gene, "symbol_id"]
	top.cd14.shared <- rbind(top.cd14.shared, matrix(c(cd14.vst.list[[1]][gene], cd14.vst.list[[2]][gene], cd14.vst.list[[3]][gene]), ncol=3, dimnames=list(symbol, c('EU','AA','EA'))))
}

write.table(top.cd14.shared, file='/home/t.cri.cczysz/cd14.vst.sex.txt', quote=F, sep='\t')
top.cd4.shared <- c()
for (gene in names(cd4.multi)) {
	symbol <- annot[annot$ensembl_id==gene, "symbol_id"]
	top.cd4.shared <- rbind(top.cd4.shared, matrix(c(cd4.vst.list[[1]][gene], cd4.vst.list[[2]][gene], cd4.vst.list[[3]][gene]), ncol=3, dimnames=list(symbol, c('EU','AA','EA'))))
}
write.table(top.cd4.shared, file='/home/t.cri.cczysz/cd4.vst.sex.txt', quote=F, sep='\t')
