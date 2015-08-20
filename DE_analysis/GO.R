library(limma)

load('/group/stranger-lab/czysz/ImmVar/annot.Robj')
fit.dir <- "/group/stranger-lab/immvar_data"

populations <- c("Caucasian","Asian","African-American")
cell.types <- c("CD14","CD4")

for (population in populations) {
for (cell.type in cell.types) {
	fit.f <- paste("fit",cell.type,population,"Robj",sep=".")
	load(file=pate(fit.dir,fit.f,sep="/"))

	ttable <- topTable(eb.fit, number = Inf)
	symbols <- c()
	for (id in rownames(ttable)) {
		symbol <- gencode_annot[gencode_annot$gene_id==id, "gene_symbol"]
		symbols <- c(symbols,symbol)
	}
	out.dir <- '/home/t.cri.cczysz/'
	out.f <- paste('sorted_genes',cell.type,population,'txt',sep='.')
	write.table(symbols,file=paste(out.dir,out.f,sep=""),row.names=F,col.names=F,quote=F)
}
}
