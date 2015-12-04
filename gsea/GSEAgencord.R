load('/group/stranger-lab/immvar_rep/mappings/annot.Robj')
gene_annot <- annots
library(biomaRt)
if (F) {
rownames(gene_annot)=gene_annot[,1]
gene_annot=gene_annot[order(as.numeric(gsub("M","25",gsub("Y","24",gsub("X","23",substring(gene_annot$chr,4)))))),]
gene_annot$chr=factor(gene_annot$chr,levels=unique(gene_annot$chr))
}

#load("Robjects/residuals.CD14.Caucasian.Robj")
#load("Robjects/phen.Robj")
if (F) {
cd14.cau.res=expr.residuals;
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
#ensembl=useMart("ensembl")
#ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl");
biomart.results=getBM(ensembl,attributes=c("ensembl_gene_id","hgnc_symbol","entrezgene"),filters="ensembl_gene_id",values=substring(rownames(cd14.cau.res),1,15))
uni=biomart.results[!duplicated(biomart.results[,1]),]
uni=uni[!duplicated(uni[,2]),]
uni=uni[!nchar(uni$hgnc_symbol)<1,]
a=cd14.cau.res[substring(rownames(cd14.cau.res),1,15)%in%uni[,1],]
a=cbind(uni$hgnc_symbol,a)
colnames(a)[1]='Name'
write.table(file="data/mRNAexp/RES_files/expr.residuals.cd14.cau.2.txt",a,quote=F,row.names=F,sep="\t")
## Manually open *.cls file to write the header after saving the group info.
## Change CellType for CD14-CD4
write.table(file="data/pheno/CD14/sex.cd14.cau.cls",t(as.matrix(as.numeric(phen[phen$Race%in%"Caucasian" & phen$CellType%in%"CD14+16-Mono","Sex"])-1)),quote=F,col.names=F,row.names=F)
dim(cd14.cau.res);
}

## preRANKED

immvar.fit.objects <- c('fit.joint.CD14.Robj', 'fit.joint.CD4.Robj')
immvar.dir <- '/group/stranger-lab/immvar_data/'
rep.fit.objects <- c('emtab2232/fairfax_fit.Robj', 'GenCord/gencord_fit.Robj', 'GSE56580/mesa_tcells_fit.Robj', 'GSE56045/gencord_fit.Robj')
rep.dir <- '/group/stranger-lab/immvar_rep/'
if (T) {
for (obj in immvar.fit.objects) {
	load(paste(immvar.dir, obj, sep=''))
	pvals <- eb.fit$p.value
	ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
	biomart.results=getBM(ensembl,attributes=c("ensembl_gene_id","hgnc_symbol","entrezgene"),filters="ensembl_gene_id",values=substring(rownames(pvals),1,15))
	uni=biomart.results[!duplicated(biomart.results[,1]),]
	uni=uni[!duplicated(uni[,2]),]
	uni=uni[!nchar(uni$hgnc_symbol)<1,]
	a=pvals[substring(rownames(pvals),1,15)%in%uni[,1],]
	rownames(uni)=names(a)
	to.write <- cbind(uni[names(a[order(a)]), "hgnc_symbol"], -log10(a[order(a)]))

	f.name <- paste('immvar', unlist(strsplit(obj, split='[.]'))[3], 'rnk', sep='.')
	write.table(file=paste("/home/t.cri.cczysz/", f.name, sep=''), cbind(uni[names(a[order(a)]),"hgnc_symbol"], -log10(a[order(a)])), quote=F, row.names=F, col.names=F, sep="\t")
}
}
if (T) {
for (obj in rep.fit.objects) {
	load(paste(rep.dir, obj, sep=''))

	pvals <- fit$p.value
	ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
	biomart.results=getBM(ensembl,attributes=c("ensembl_gene_id","hgnc_symbol","entrezgene"),filters="ensembl_gene_id",values=substring(rownames(pvals),1,15))
	uni=biomart.results[!duplicated(biomart.results[,1]),]
	uni=uni[!duplicated(uni[,2]),]
	uni=uni[!nchar(uni$hgnc_symbol)<1,]
	a=pvals[substring(rownames(pvals),1,15)%in%uni[,1],]
	rownames(uni)=names(a)

	f.name <- paste('rep', unlist(strsplit(obj, split='[/]'))[1], 'rnk', sep='.')
	write.table(file=paste("/home/t.cri.cczysz/", f.name, sep=''), cbind(uni[names(a[order(a)]),"hgnc_symbol"], -log10(a[order(a)])), quote=F, row.names=F, col.names=F, sep="\t")
}
}

if (T) {
setwd('/group/stranger-lab/immvar/meta/')
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
for (f in c('cd14_meta_rep1.txt', 'cd4_meta_rep1.txt')) {
	meta_res <- read.table(file=f, header=T, row.names=1)
	pvals <- data.frame(meta_res$P.value, row.names=rownames(meta_res))
		
	biomart.results=getBM(ensembl,attributes=c("ensembl_gene_id","hgnc_symbol","entrezgene"),filters="ensembl_gene_id",values=substring(rownames(meta_res),1,15))
	uni=biomart.results[!duplicated(biomart.results[,1]),]
	uni=uni[!duplicated(uni[,2]),]
	uni=uni[!nchar(uni$hgnc_symbol)<1,]
	a=pvals[substring(rownames(pvals),1,15)%in%uni[,1],]
	names(a) <- rownames(pvals)[substring(rownames(pvals), 1, 15)%in%uni[,1]]
	rownames(uni) <- uni[, 1]
	a.order <- a[order(a)]	
	to.write <- cbind(uni[substr(names(a[order(a)]),1,15), "hgnc_symbol"], -log10(a[order(a)]))
	f.name <- paste('meta', unlist(strsplit(f, split='_'))[1], 'rnk', sep='.')
	write.table(to.write, file=paste('/home/t.cri.cczysz/', f.name, sep=''), quote=F, row.names=F, col.names=F, sep="\t")
}
}
