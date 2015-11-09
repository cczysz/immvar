library(limma)

wd='/group/stranger-lab/immvar_data/'
setwd(wd)

for (pop in c('Caucasian', 'African-American', 'Asian')) {
for (cell in c('CD14', 'CD4')) {

f_in <- paste('fit', pop, cell, 'Robj', sep='.')
load(file=f_in)

if (cell == 'CD14') { 
	ttable<-topTable(eb.fit, number=Inf)
	sig<-topTable(eb.fit, number=Inf, p.value=0.05)
	n<-nrow(cd14.eb.fit$design)
	n_genes<-nrow(ttable) } 
else {
	ttable<-topTable(eb.fit, number=Inf) 
	sig<-topTable(eb.fit, number=Inf, p.value=0.05)
	n<-nrow(cd4.eb.fit$design)
	n_genes<-nrow(ttable)
}

out_matrix <- matrix(
	c(GENES=rownames(ttable),
	REF=rep('G', times=n_genes),
	ALT=rep('A', times=n_genes), 
	STRAND=rep('+', times=n_genes), 
	EFFECT=ttable$logFC,
	PVALUE=ttable$P.Value, 
	N=eb.fit$N,
	ncol=7)
colnames(out_matrix) <- c('GENES', 'REF', 'ALT', 'STRAND', 'EFFECT', 'PVALUE', 'N')

f_out <- paste('meta_f', pop, cell, 'txt', sep='.')
write.table(out_matrix, file=paste('/home/t.cri.cczysz/',f_out,sep=''), quote=F, row.names=F)

f_out <- paste('sig_de', pop, cell, 'txt', sep='.')
write.table(sig, file=paste('/home/t.cri.cczysz/',f_out,sep=''), quote=F, sep="\t")

}
}
