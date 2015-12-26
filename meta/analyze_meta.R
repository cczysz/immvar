library(limma)
library(stats)
library(VennDiagram)
library(gtools)
library(ggplot2)
library(reshape2)

setwd('/group/stranger-lab/immvar/meta')
load('/group/stranger-lab/immvar_rep/mappings/annot.Robj')

makePlots <- function(fit, cell.type, annots) {
	
	qqplot.data <- data.frame(obs=sort(fit$p.value), exp=seq(length(fit$p.value))/length(fit$p.value))
	g <- ggplot(data=qqplot.data, aes(x=exp, y=obs))
	g + geom_point()
}
########################################
# Process jointly normalized ImmVar data
load('/group/stranger-lab/immvar_data/fit.joint.CD14.Robj')
joint.cd14.fit <- data.frame(eb.fit)
joint.cd14.fit$q.value <- p.adjust(joint.cd14.fit$p.value, method="fdr")
joint.cd14.fit$rank <- rank(joint.cd14.fit$p.value) / nrow(joint.cd14.fit)
joint.cd14.fit$chr <- as.character(annots[rownames(joint.cd14.fit),'chr'])
joint.cd14.sig <- subset(joint.cd14.fit, q.value<=0.05)

sig.5 <- max(joint.cd14.sig$p.value)
pdf(file=paste('/group/stranger-lab/czysz/ImmVar/plots/','qq_CD14.pdf',sep=''))
qqplot.data <- data.frame(obs=-log10(sort(joint.cd14.fit$p.value)),
	exp=-log10(seq(length(joint.cd14.fit$p.value))/length(joint.cd14.fit$p.value)),
	chr=as.character(joint.cd14.fit$chr[order(joint.cd14.fit$rank)]),
	row.names=row.names(joint.cd14.fit[order(joint.cd14.fit$p.value),]))
qqplot.data$chr <- as.character(qqplot.data$chr)
qqplot.data[!(qqplot.data$chr=='chrY' | qqplot.data$chr=='chrX'), 'chr'] <- 'auto'
g <- ggplot(data=qqplot.data, aes(x=exp, y=obs, color=chr))
g + geom_point() + labs(title="CD14", x='Expected', y='Observed') + 
	scale_color_manual(values=c('black', 'red','blue'),
		name='Chromosome', breaks=c('auto', 'chrX', 'chrY'), labels=c('Auto', 'X', 'Y')) + 
	geom_abline()
dev.off()

pdf(file=paste('/group/stranger-lab/czysz/ImmVar/plots/','volcano_CD14.pdf',sep=''))
volplot.data <- data.frame(logFC=joint.cd14.fit$coefficients,
	pval=-log10(joint.cd14.fit$p.value),
	chr=as.character(joint.cd14.fit$chr),
	row.names=rownames(joint.cd14.fit))
volplot.data$chr <- as.character(volplot.data$chr)
volplot.data[!(volplot.data$chr=='chrY' | volplot.data$chr=='chrX'),'chr'] <- 'auto'
g <- ggplot(data=volplot.data, aes(x=logFC, y=pval, color=chr))
g + geom_point() + labs(title='ImmVar CD14', x='log2FC', y='-log10(pvalue)') + 
	scale_color_manual(values=c('black', 'red','blue'),
		name='Chromosome', breaks=c('auto','chrX', 'chrY'), labels=c('Auto','X', 'Y')) + 
	geom_abline(slope=0, intercept=-log10(sig.5))

if (F) {
g + geom_point() + labs(title='ImmVar CD14-Zoomed', x='log2FC', y='-log10(pvalue)') + 
	scale_color_manual(values=c('black', 'red','blue'),
		name='Chromosome', breaks=c('auto','chrX', 'chrY'), labels=c('Auto','X', 'Y')) + 
	scale_x_continuous(limits = c(-1, 1)) + scale_y_continuous(limits = c(0, 30)) +
	geom_abline(slope=0, intercept=-log10(sig.5))
}
dev.off()

load('/group/stranger-lab/immvar_data/fit.joint.CD4.Robj')
print(topTable(eb.fit))
joint.cd4.fit <- data.frame(eb.fit)
joint.cd4.fit$q.value <- p.adjust(joint.cd4.fit$p.value, method="fdr")
joint.cd4.fit$rank <- rank(joint.cd4.fit$p.value) / nrow(joint.cd4.fit)
joint.cd4.fit$chr <- as.character(annots[rownames(joint.cd4.fit),'chr'])
joint.cd4.sig <- subset(joint.cd4.fit, q.value<=0.05)

sig.5 <- max(joint.cd14.sig$p.value)

pdf(file=paste('/group/stranger-lab/czysz/ImmVar/plots/','qq_CD4.pdf',sep=''))
qqplot.data <- data.frame(obs=-log10(sort(joint.cd4.fit$p.value)),
	exp=-log10(seq(length(joint.cd4.fit$p.value))/length(joint.cd4.fit$p.value)),
	chr=as.character(joint.cd14.fit$chr[order(joint.cd14.fit$rank)]),
	row.names=row.names(joint.cd14.fit[order(joint.cd14.fit$p.value),]))
qqplot.data$chr <- as.character(qqplot.data$chr)
qqplot.data[!(qqplot.data$chr=='chrY' | qqplot.data$chr=='chrX'),'chr'] <- 'auto'

g <- ggplot(data=qqplot.data, aes(x=exp, y=obs, color=chr))
g + geom_point() + labs(title="CD4", x='Expected', y='Observed') + 
	scale_color_manual(values=c('black','red','blue'),
		name='Chromosome', breaks=c('auto', 'chrX', 'chrY'), labels=c('Auto','X', 'Y')) + 
	geom_abline()
dev.off()

pdf(file=paste('/group/stranger-lab/czysz/ImmVar/plots/','volcano_CD4.pdf',sep=''))
qqplot.data <- data.frame(logFC=joint.cd4.fit$coefficients,
	pval=-log10(joint.cd4.fit$p.value),
	chr=joint.cd4.fit$chr,
	row.names=rownames(joint.cd4.fit))
qqplot.data$chr <- as.character(qqplot.data$chr)
qqplot.data[!(qqplot.data$chr=='chrY' | qqplot.data$chr=='chrX'),'chr'] <- 'auto'
g <- ggplot(data=qqplot.data, aes(x=logFC, y=pval, color=chr))
g + geom_point() + labs(title='ImmVar CD4', x='log2FC', y='-log10(pvalue)') + 
	scale_color_manual(values=c('black','red','blue'),
		name='Chromosome', breaks=c('auto','chrX', 'chrY'), labels=c('Auto','X', 'Y')) + 
	geom_abline(slope=0, intercept=-log10(sig.5[1]))

g + geom_point() + labs(title='ImmVar CD4-Zoomed', xlab='log2FC', ylab='-log10(pvalue)') + 
	scale_color_manual(values=c('black','red','blue'),
		name='Chromosome', breaks=c('auto','chrX', 'chrY'), labels=c('Auto','X', 'Y')) + 
	scale_x_continuous(limits = c(-1, 1)) + scale_y_continuous(limits = c(0, 30)) +
	geom_abline(slope=0, intercept=-log10(sig.5[1]))
dev.off()

chr.counts <- data.frame(cd4=table(annots[rownames(joint.cd4.sig), 'chr']), cd14=table(annots[rownames(joint.cd14.sig), 'chr']))
	colnames(chr.counts) <- c('Chr', 'CD4', 'Chr', 'CD14')
	chr.counts <- chr.counts[mixedorder(chr.counts$Chr), ]
	chr.counts[, c(1,3)] <- factor(as.character(chr.counts[,1]), levels=as.character(chr.counts[, 1]))
	chr.counts <- melt(chr.counts, id.vars='Chr', measure.vars=c('CD4', 'CD14'))
	chr.counts[, 1] <- factor(as.character(chr.counts[,1]), levels=as.character(chr.counts[,1]))
pdf(file='/group/stranger-lab/czysz/ImmVar/plots/joint.chr.counts.pdf', width=11, height=8)
g <- ggplot(chr.counts, aes(x=Chr, y=value, fill=variable, levels=Chr))
g + geom_bar(stat="identity", position='dodge')
dev.off()

chr.percents <- data.frame(cd4=table(annots[rownames(joint.cd4.sig), 'chr'])/table(annots[rownames(joint.cd4.fit), 'chr']),
		cd14=table(annots[rownames(joint.cd14.sig), 'chr'])/table(annots[rownames(joint.cd14.fit), 'chr']))
	colnames(chr.percents) <- c('Chr', 'CD4', 'Chr', 'CD14')
	chr.percents <- chr.percents[mixedorder(chr.percents$Chr), ]
	chr.percents[, c(1,3)] <- factor(as.character(chr.percents[,1]), levels=as.character(chr.percents[, 1]))
	chr.percents <- melt(chr.percents, id.vars='Chr', measure.vars=c('CD4', 'CD14'))
	chr.percents[, 1] <- factor(as.character(chr.percents[,1]), levels=as.character(chr.percents[,1]))

pdf(file='/group/stranger-lab/czysz/ImmVar/plots/joint.chr.percents.pdf', width=11, height=8)
g <- ggplot(chr.percents, aes(x=Chr, y=value, fill=variable, levels=Chr))
g + geom_bar(stat="identity", position='dodge')
dev.off()

joint.shared.sig <- intersect(rownames(joint.cd14.sig), rownames(joint.cd4.sig))

joint.venn <- list(CD4 = rownames(joint.cd4.sig), CD14=rownames(joint.cd14.sig))

###############################################
# Process meta-analysis of replication datasets
setwd('/group/stranger-lab/immvar/meta/')

cd14.rep.meta <- read.table(file='cd14_meta_all1.txt', header=T, row.names=1)
cd14.rep.meta <- subset(cd14.rep.meta, Weight==2097)
cd4.rep.meta <- read.table(file='cd4_meta_all1.txt', header=T, row.names=1)
cd4.rep.meta <- subset(cd4.rep.meta, Weight==719)

cd14.rep.meta <- cbind(cd14.rep.meta, 
	q.value=p.adjust(cd14.rep.meta$P.value, method='fdr'),
	rank=rank(cd14.rep.meta$P.value) / nrow(cd14.rep.meta),
	chr=annots[rownames(cd14.rep.meta), 'chr'],
	gene=annots[rownames(cd14.rep.meta), "symbol_id"])
cd14.rep.meta$chr <- as.character(cd14.rep.meta$chr)
cd14.rep.sig.5 <- subset(cd14.rep.meta, q.value<0.05)

pdf(file='/group/stranger-lab/czysz/ImmVar/plots/cd14.meta.qq.pdf')
g <- ggplot(data=cd14.rep.meta, aes(x=-log10(seq(nrow(cd14.rep.meta))/nrow(cd14.rep.meta)), y=-log10(sort(P.value)), color=replace(cd14.rep.meta$chr, !(cd14.rep.meta$chr=='chrX' | cd14.rep.meta$chr=='chrY'), 'auto')[order(cd14.rep.meta$rank)]))
g + geom_point() + labs(title="CD14 All Meta", x='Expected', y='Observed') + 
	scale_color_manual(values=c('black', 'red','blue'),
		name='Chromosome', breaks=c('auto', 'chrX', 'chrY'), labels=c('Auto', 'X', 'Y')) + 
	geom_abline()
dev.off()

pdf(file='/group/stranger-lab/czysz/ImmVar/plots/cd14.meta.volcano.pdf')
g <- ggplot(data=cd14.rep.meta, aes(x=Zscore, y=-log10(P.value), color=replace(cd14.rep.meta$chr, !(cd14.rep.meta$chr=='chrX' | cd14.rep.meta$chr=='chrY'), 'auto'))) 
g + geom_point() + labs(main='CD14 All Meta', x='Expected', y='Observed') +
	scale_color_manual(values=c('black', 'red','blue'),
		name='Chromosome', breaks=c('auto', 'chrX', 'chrY'), labels=c('Auto', 'X', 'Y')) + 
	geom_abline(slope=0, intercept=max(cd14.rep.sig.5$P.value))
dev.off()

cd4.rep.meta <- cbind(cd4.rep.meta, 
	q.value=p.adjust(cd4.rep.meta$P.value, method='fdr'),
	rank=rank(cd4.rep.meta$P.value) / nrow(cd4.rep.meta),
	chr=annots[rownames(cd4.rep.meta), 'chr'],
	gene=annots[rownames(cd4.rep.meta), "symbol_id"])
cd4.rep.meta$chr <- as.character(cd4.rep.meta$chr)
cd4.rep.sig.5 <- subset(cd4.rep.meta, q.value<0.05)

pdf(file='/group/stranger-lab/czysz/ImmVar/plots/cd4.meta.qq.pdf')
g <- ggplot(data=cd4.rep.meta, aes(x=-log10(seq(nrow(cd4.rep.meta))/nrow(cd4.rep.meta)), y=-log10(sort(P.value)), color=replace(cd4.rep.meta$chr, !(cd4.rep.meta$chr=='chrX' | cd4.rep.meta$chr=='chrY'), 'auto')[order(cd4.rep.meta$rank)]))
g + geom_point() + labs(title="CD4 All Meta", x='Expected', y='Observed') + 
	scale_color_manual(values=c('black', 'red','blue'),
		name='Chromosome', breaks=c('auto', 'chrX', 'chrY'), labels=c('Auto', 'X', 'Y')) + 
	geom_abline()
dev.off()

pdf(file='/group/stranger-lab/czysz/ImmVar/plots/cd4.meta.volcano.pdf')
#g <- ggplot(data=cd4.rep.meta, aes(x=Zscore, y=-log10(P.value), color=replace(cd4.rep.meta$chr, !(cd4.rep.meta$chr=='chrX' | cd4.rep.meta$chr=='chrY'), 'auto'))) 
g <- ggplot(data=cd4.rep.meta, aes(x=Zscore, y=-log10(P.value), color=replace(chr, !(chr=='chrX' | chr=='chrY'), 'auto'))) 
g + geom_point() + labs(main='CD4 All Meta', x='Zscore', y='-log10(P.value)') +
	scale_color_manual(values=c('black', 'red','blue'),
		name='Chromosome', breaks=c('auto', 'chrX', 'chrY'), labels=c('Auto', 'X', 'Y')) + 
	geom_abline(slope=0, intercept=max(cd14.rep.sig.5$P.value))
dev.off()

rep.cd4.sig <- list()
rep.cd14.sig <- list()

cd14.overlap <- c()
cd4.overlap <- c()
for (fdr in c(0.05, 0.1, 0.15, 0.2)) {
	rep.cd14.sig[[as.character(fdr)]] <- rownames(subset(cd14.rep.meta,q.value <= fdr))
	rep.cd4.sig[[as.character(fdr)]]  <- rownames(subset(cd4.rep.meta,q.value <= fdr))

	cd14.overlap <- c(cd14.overlap, length(intersect(rep.cd14.sig[[fdr]], rownames(joint.cd14.sig))))	
	cd4.overlap <- c(cd4.overlap, length(intersect(rep.cd4.sig[[fdr]], rownames(joint.cd4.sig))))	
}

rep.genes.cd4 <- list()
rep.genes.cd14 <- list()
for (i in seq(length(rep.cd14.sig))) {
	rep.genes.cd4[[i]] <- intersect(rep.cd4.sig[[i]], rownames(joint.cd4.sig))
	rep.genes.cd14[[i]] <- intersect(rep.cd14.sig[[i]], rownames(joint.cd14.sig))
}

pdf(file='/group/stranger-lab/czysz/ImmVar/plots/cd4.volcano.replicating.pdf')
ggplot(data=cd4.rep.meta, aes(x=Zscore, y=-log10(P.value))) + geom_point() + 
	geom_point(data=cd4.rep.meta[rownames(cd4.rep.meta)%in%rep.genes.cd4[[2]],], col='red') +
	labs(title='CD4 Meta Overlap', x='log2FC', ylab='-log10(pvalue)')
dev.off()

pdf(file='/group/stranger-lab/czysz/ImmVar/plots/cd14.volcano.replicating.pdf')
ggplot(data=cd14.rep.meta, aes(x=Zscore, y=-log10(P.value))) + geom_point() + 
	geom_point(data=cd14.rep.meta[rownames(cd14.rep.meta)%in%rep.genes.cd14[[2]],], col='red') +
	labs(title='CD14 Meta Overlap', x='log2FC', ylab='-log10(pvalue)')
dev.off()
# Number of genes significant in both replication and discovery
rep.counts <- data.frame(CD14=unlist(lapply(rep.genes.cd14, length)),
	CD4=unlist(lapply(rep.genes.cd4, length)),
	fdr=c("0.05", "0.1", "0.15", "0.2"))

# Percent of significant genes in discovery set also significant in replication sets
rep.rep.percents <- data.frame(CD14=unlist(lapply(rep.genes.cd14, length))/unlist(lapply(rep.cd14.sig, length)),
	CD4=unlist(lapply(rep.genes.cd4, length))/unlist(lapply(rep.cd4.sig, length)),
	fdr=c("0.05", "0.1", "0.15", "0.2"))

# Percent of significant genes in discovery set also significant in replication sets
rep.dis.percents <- data.frame(CD14=unlist(lapply(rep.genes.cd14, length))/nrow(joint.cd14.sig),
	CD4=unlist(lapply(rep.genes.cd4, length))/nrow(joint.cd4.sig),
	fdr=c("0.05", "0.1", "0.15", "0.2"))

pdf(file='/group/stranger-lab/czysz/ImmVar/plots/rep.counts.pdf')
ggplot(data=melt(rep.counts), aes(x=fdr, y=value, fill=variable)) +
	geom_bar(stat="identity", position='dodge') +
	ggtitle('Counts of replicating genes at different FDR')

ggplot(data=melt(rep.rep.percents), aes(x=fdr, y=value, fill=variable)) + 
	geom_bar(stat="identity", position='dodge') +
	ggtitle('Intersect significant discovery genes over total significant in replication at different FDR')

ggplot(data=melt(rep.dis.percents), aes(x=fdr, y=value, fill=variable)) +
	geom_bar(stat="identity", position='dodge') +
	ggtitle('Intersect significant discovery genes over total significant in discovery set at different FDR')
dev.off()

if (T) {

#joint.venn <- list(CD4 = rownames(joint.cd4.sig), CD14=rownames(joint.cd14.sig))
venn.diagram(joint.venn,
	 filename='/group/stranger-lab/czysz/joint_venn.tiff',
	 fontfamily="Helvetica",
	 main.fontfamily="Helvetica",
	 cat.fontfamily="Helvetica",
	 sub.fontfamily="Helvetica",
	 fill=topo.colors(length(joint.venn)),
	 main="Joint CD4 vs CD14 VennDiagram", sub.cex=1.5,
	 width=10, height=10, units="in")

meta.cd4.venn <- list(meta=rep.cd4.sig[[2]], joint=rownames(joint.cd4.sig))

venn.diagram(meta.cd4.venn,
	 filename='/group/stranger-lab/czysz/cd4_venn.tiff',
	 fontfamily="Helvetica",
	 main.fontfamily="Helvetica",
	 cat.fontfamily="Helvetica",
	 sub.fontfamily="Helvetica",
	 fill=topo.colors(length(meta.cd4.venn)),
	 main="CD4 VennDiagram\nSex-biased genes (FDR<5%)", sub.cex=1.5,
	 width=10, height=10, units="in")

meta.cd14.venn <- list(meta=rep.cd14.sig[[2]], joint=rownames(joint.cd14.sig))
venn.diagram(meta.cd14.venn,
	 filename='/group/stranger-lab/czysz/cd14_venn.tiff', main.fontfamily="Helvetica", sub.fontfamily="Helvetica",
	 fontfamily="Helvetica",
	 cat.fontfamily="Helvetica",
	 fill=topo.colors(length(meta.cd14.venn)),
	 main="CD14 VennDiagram\nSex-biased genes (FDR<5%)", sub.cex=1.5,
	 width=10, height=10, units="in")
}
