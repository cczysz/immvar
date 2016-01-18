library(limma)
library(ggplot2)

load('/group/stranger-lab/immvar_rep/mappings/annot.Robj')
wd='/group/stranger-lab/immvar_data/'
setwd(wd)

#for (pop in c('Caucasian', 'African-American', 'Asian')) {
for (cell in c('CD14', 'CD4')) {

f_in <- paste('fit.joint', cell, 'Robj', sep='.')
load(file=f_in)

if (cell == 'CD14') { 
	ttable<-topTable(eb.fit, number=Inf)
	cd14.genes <- rownames(ttable)
	sig<-topTable(eb.fit, number=Inf, p.value=0.05)
	n<-nrow(eb.fit$design)
	n_genes<-nrow(ttable) } 
else {
	ttable<-topTable(eb.fit, number=Inf) 
	cd4.genes <- rownames(ttable)
	sig<-topTable(eb.fit, number=Inf, p.value=0.05)
	n<-nrow(eb.fit$design)
	n_genes<-nrow(ttable)
}

out_matrix <- matrix(
	c(GENES=rownames(ttable),
	REF=rep('G', times=n_genes),
	ALT=rep('A', times=n_genes), 
	STRAND=rep('+', times=n_genes), 
	EFFECT=ttable$logFC,
	PVALUE=ttable$P.Value, 
	N=eb.fit$N),
	ncol=7)
colnames(out_matrix) <- c('GENES', 'REF', 'ALT', 'STRAND', 'EFFECT', 'PVALUE', 'N')

f_out <- paste('meta_f.joint', cell, 'txt', sep='.')
write.table(out_matrix, file=paste('/home/t.cri.cczysz/',f_out,sep=''), quote=F, row.names=F)
}
#}

cd14.files.in <- c('emtab2232/fairfax_fit.Robj', 'GSE56045/gencord_fit.Robj')
cd4.files.in <- c('GenCord/gencord_fit.Robj', 'GSE56580/mesa_tcells_fit.Robj')
cd14.studies <- c('Fairfax', 'MesaM')
cd4.studies <- c('Gencord', 'MesaT')
cd14.sample_size <- as.numeric(c(432, 1264))
cd4.sample_size <- as.numeric(c(85, 227))
# fit,fit
fit_dir <- "/group/stranger-lab/immvar_rep/"

for (i in seq(length(cd14.files.in))) {
	load(paste(fit_dir, cd14.files.in[i], sep=''))
	
	ttable <- topTable(fit, number=Inf)
#	ttable <- ttable[rownames(ttable)%in%cd14.genes,]
	n_genes <- nrow(ttable)
	rm(out_matrix)
	out_matrix <- data.frame(
		GENES=rownames(ttable),
		REF=rep('G', times=n_genes),
		ALT=rep('A', times=n_genes), 
		STRAND=rep('+', times=n_genes), 
		EFFECT=ttable$logFC,
		PVALUE=ttable$P.Value, 
		N=cd14.sample_size[i])
	colnames(out_matrix) <- c('GENES', 'REF', 'ALT', 'STRAND', 'EFFECT', 'PVALUE', 'N')
	print(head(out_matrix))

	f_out <- paste('meta_f', cd14.studies[i], 'txt', sep='.')
	write.table(out_matrix, file=paste('/home/t.cri.cczysz/',f_out,sep=''), quote=F, row.names=F)
}
for (i in seq(length(cd4.files.in))) {
	load(paste(fit_dir, cd4.files.in[i], sep=''))
	
	ttable <- topTable(fit, number=Inf)
	#ttable <- ttable[rownames(ttable)%in%cd4.genes,]
	n_genes <- nrow(ttable)
	rm(out_matrix)
	out_matrix <- data.frame(
		GENES=rownames(ttable),
		REF=rep('G', times=n_genes),
		ALT=rep('A', times=n_genes), 
		STRAND=rep('+', times=n_genes), 
		EFFECT=ttable$logFC,
		PVALUE=ttable$P.Value, 
		N=cd4.sample_size[i])
	colnames(out_matrix) <- c('GENES', 'REF', 'ALT', 'STRAND', 'EFFECT', 'PVALUE', 'N')
	print(head(out_matrix))

	f_out <- paste('meta_f', cd4.studies[i], 'txt', sep='.')
	write.table(out_matrix, file=paste('/home/t.cri.cczysz/',f_out,sep=''), quote=F, row.names=F)
}
	if (F) {
	rep.fit <- data.frame(fit)
	rep.fit$chr <- as.character(annots[rownames(fit),'chr'])
	rep.fit[!(rep.fit$chr=='chrY' | rep.fit$chr=='chrX'), 'chr'] <- 'auto'
	# QQ Plot	
	pdf(file=paste('/group/stranger-lab/czysz/ImmVar/plots/', paste(studies[i], 'qq.pdf', sep='_'), sep=''))
	#print(qplot(data=rep.fit, x=-log10(seq(nrow(rep.fit))/nrow(rep.fit)), y=-log10(sort(p.value)), main=studies[i], xlab='expected', ylab='observed') + geom_abline())
	rep.fit.ordered <- cbind(rep.fit[order(rep.fit$p.value),], null.rank=seq(nrow(rep.fit))/nrow(rep.fit))
	g <- ggplot(data=rep.fit.ordered, aes(x=-log10(null.rank), y=-log10(p.value), fill=chr, color=chr))
	print(g + geom_point() + labs(title=studies[i], x='expected', y='observed') + geom_abline() + guides(fill=FALSE) + scale_color_manual(values=c("black", "red", "blue"), 
                       name="Chromosome",
                       breaks=c("auto", "chrX", "chrY"),
                       labels=c("Autosome", "X", "Y"))
		#geom_point(data=rep.fit.ordered[rownames(rep.fit.ordered)%in%rownames(subset(annots, chr=='chrY')), ], col='blue') +
		#geom_point(data=rep.fit.ordered[rownames(rep.fit.ordered)%in%rownames(subset(annots, chr=='chrX')), ], col='red') 
	)
	#print(qplot(data=rep.fit, x=coefficients, y=-log10(p.value+1e-350), xlab='log2FC', ylab='-log10(pvalue)', main=studies[i]))
	dev.off()
	# Volcano Plot
	q.vals <- p.adjust(rep.fit$p.value, method='fdr')
	sig.5 <- rep.fit$p.value[which(q.vals==min(q.vals[q.vals>=0.05]))][1]
	sig.10 <- rep.fit$p.value[which(q.vals==min(q.vals[q.vals>=0.1]))][1]
	pdf(file=paste('/group/stranger-lab/czysz/ImmVar/plots/', paste(studies[i], 'volcano.pdf', sep='_'), sep=''))
	print(ggplot(rep.fit, aes(x=coefficients, y=-log10(p.value), color=chr)) + geom_point() + 
		labs(title=studies[i], x='log2FC', y= '-log10(pvalue)') + geom_abline(intercept=-log10(sig.10), slope=0) + 
		guides(fill=FALSE) + scale_color_manual(values=c("black", "red", "blue"), 
                       name="Chromosome",
                       breaks=c("auto", "chrX", "chrY"),
                       labels=c("Autosome", "X", "Y"))

	)
	print(ggplot(rep.fit, aes(x=coefficients, y=-log10(p.value))) + geom_point() + 
		geom_point(data=rep.fit[rownames(rep.fit)%in%rownames(subset(annots, chr=='chrY')), ], col='blue') +
		geom_point(data=rep.fit[rownames(rep.fit)%in%rownames(subset(annots, chr=='chrX')), ], col='red') +
		labs(title=studies[i], x='log2FC', y= '-log10(pvalue)') + scale_x_continuous(limits = c(-1, 1)) + scale_y_continuous(limits = c(0, 30)) +
		geom_abline(intercept=-log10(sig.10), slope=0)

	)
	#print(qplot(data=rep.fit, x=coefficients, y=-log10(p.value+1e-350), xlab='log2FC', ylab='-log10(pvalue)', main=studies[i]))
	dev.off() }
#+ geom_point(data=na.omit(fit[rownames(subset(annots, chr=='chrY')), ]), col='red') + 
		#geom_point(data=na.omit(fit[rownames(subset(annots, chr=='chrX')), ]), col='blue')

#dev.off()
