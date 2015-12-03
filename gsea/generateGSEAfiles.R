gene_annot=read.table("probes_mapping/annotations/gencode.v22.TSS.txt",header=T)

library(biomaRt)
rownames(gene_annot)=gene_annot[,1]
gene_annot=gene_annot[order(as.numeric(gsub("M","25",gsub("Y","24",gsub("X","23",substring(gene_annot$chr,4)))))),]
gene_annot$chr=factor(gene_annot$chr,levels=unique(gene_annot$chr))

load("Robjects/residuals.CD14.Caucasian.Robj")
load("Robjects/phen.Robj")
cd14.cau.res=expr.residuals;
ensembl=useMart("ensembl")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl");
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

q()


## preRANKED

cd14.cau.pvals=as.matrix(read.table("DE_gender/summary_pvals.CD14.Caucasian.txt",header=T))
ensembl=useMart("ensembl")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl");
biomart.results=getBM(ensembl,attributes=c("ensembl_gene_id","hgnc_symbol","entrezgene"),filters="ensembl_gene_id",values=substring(rownames(cd14.cau.pvals),1,15))
uni=biomart.results[!duplicated(biomart.results[,1]),]
uni=uni[!duplicated(uni[,2]),]
uni=uni[!nchar(uni$hgnc_symbol)<1,]
a=cd14.cau.pvals[substring(rownames(cd14.cau.pvals),1,15)%in%uni[,1],]
rownames(uni)=rownames(a)
write.table(file="GSEA/CD14/CAU/cd14.cau.rank_candidates_perm.rnk",cbind(uni[rownames(a[order(a[,"P.Value"]),]),"hgnc_symbol"],-log10(a[order(a[,"P.Value"]),"P.Value"])),quote=F,row.names=F,col.names=F,sep="\t")
