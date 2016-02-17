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

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
	# Code from R cookbook
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

xesc.genes <- read.table(file='/home/t.cri.cczysz/likely_escape_genes.txt', header=F)
xesc.genes <- as.character(xesc.genes[,1])
########################################
# Process jointly normalized ImmVar data
load('/group/stranger-lab/immvar_data/fit.joint.CD14.Robj')
joint.cd14.fit <- data.frame(eb.fit)
joint.cd14.fit$q.value <- p.adjust(joint.cd14.fit$p.value, method="fdr")
joint.cd14.fit$rank <- rank(joint.cd14.fit$p.value) / nrow(joint.cd14.fit)
joint.cd14.fit$chr <- as.character(annots[rownames(joint.cd14.fit),'chr'])
joint.cd14.sig <- subset(joint.cd14.fit, q.value<=0.05)

sig.5 <- max(joint.cd14.sig$p.value)

pdf(file='/group/stranger-lab/czysz/ImmVar/plots/immvar.cd14.fc.pdf', width=11, height=8)
#plot(density(joint.cd14.fit$coefficients), main='CD14 density of log2 fold change for all genes', xlab='log2FC', ylab='')
#plot(density(joint.cd14.sig$coefficients), main='CD14 density of log2 fold change for significant genes', xlab='log2FC', ylab='')
ggplot(data.frame(log2fc = joint.cd14.fit$coefficients), aes(x=log2fc)) + geom_density() + labs(title='', x='log2FC', y='')
ggplot(data.frame(log2fc = joint.cd14.sig$coefficients), aes(x=log2fc)) + geom_density() + labs(title='', x='log2FC', y='')

dev.off()

pdf(file='/group/stranger-lab/czysz/ImmVar/plots/cd14.sig.chr.pdf', width=11, height=8)
cd14.bar.df <- data.frame(chr=names(table(joint.cd14.sig$chr)), count=as.numeric(table(joint.cd14.sig$chr)))
cd14.bar.df$count = cd14.bar.df$count[mixedorder(cd14.bar.df$chr)]
cd14.bar.df$chr = factor(mixedsort(cd14.bar.df$chr), levels=as.character(mixedsort(cd14.bar.df$chr)))
g <- ggplot(data=cd14.bar.df[-23,], aes(x=chr, y=count, fill=chr))

cd14.male <- c(table(subset(joint.cd14.sig,coefficients>0)$chr), chrM=0) 
cd14.male <- cd14.male[mixedorder(names(cd14.male))]
cd14.female <- c(table(subset(joint.cd14.sig,coefficients<0)$chr), chrM=0, chrY=0)
cd14.female <- cd14.female[mixedorder(names(cd14.female))]
cd14.bar.df <- data.frame(chr=names(cd14.male), male=as.numeric(cd14.male), female=as.numeric(cd14.female))
cd14.bar.df$chr = factor(mixedsort(cd14.bar.df$chr), levels=as.character(mixedsort(cd14.bar.df$chr)))
g <- ggplot(data=melt(cd14.bar.df), aes(x=chr, y=value, fill=variable))
g + geom_bar(stat='identity') + labs(title='', x='chromosome', y='count') +
	guides(fill=guide_legend(title="Sex")) + 
	scale_fill_manual(name="Sex", breaks=c("female", "male"), labels=c("Female", "Male"), values=c('blue', 'red'))

cd14.bar.allcount <- as.numeric(table(joint.cd14.fit$chr))[mixedorder(names(table(joint.cd14.fit$chr)))]
cd14.bar.df$male <- cd14.bar.df$male / cd14.bar.allcount
cd14.bar.df$female <- cd14.bar.df$female / cd14.bar.allcount
g <- ggplot(data=melt(cd14.bar.df), aes(x=chr, y=value, fill=variable))
g + geom_bar(stat='identity') + labs(title='', x='chromosome', y='percent') + guides(fill=guide_legend(title='Sex')) + 
	scale_fill_manual(name="Sex", breaks=c("female", "male"), labels=c("Female", "Male"), values=c('blue', 'red'))
dev.off()

pdf(file=paste('/group/stranger-lab/czysz/ImmVar/plots/','qq_CD14.pdf',sep=''), width=11, height=8)
qqplot.data <- data.frame(obs=-log10(sort(joint.cd14.fit$p.value)),
	exp=-log10(seq(length(joint.cd14.fit$p.value))/length(joint.cd14.fit$p.value)),
	chr=as.character(joint.cd14.fit$chr[order(joint.cd14.fit$rank)]),
	row.names=row.names(joint.cd14.fit[order(joint.cd14.fit$p.value),]))
qqplot.data$chr <- as.character(qqplot.data$chr)
qqplot.data[!(qqplot.data$chr=='chrY' | qqplot.data$chr=='chrX'), 'chr'] <- 'auto'
g <- ggplot(data=qqplot.data, aes(x=exp, y=obs, color=chr))
g + geom_point() + labs(title="", x='Expected -log10(p-value)', y='Observed -log10(p-value)') + 
	scale_color_manual(values=c('black', 'red','blue'),
		name='Chromosome', breaks=c('auto', 'chrX', 'chrY'), labels=c('Auto', 'X', 'Y')) + 
	geom_abline()
dev.off()

pdf(file=paste('/group/stranger-lab/czysz/ImmVar/plots/','volcano_CD14.pdf',sep=''), width=12, height=8)
volplot.data <- data.frame(logFC=joint.cd14.fit$coefficients,
	pval=-log10(joint.cd14.fit$p.value),
	chr=as.character(joint.cd14.fit$chr),
	row.names=rownames(joint.cd14.fit),
	xesc=joint.cd14.fit$symbol%in%xesc.genes)
volplot.data$chr <- as.character(volplot.data$chr)
volplot.data[!(volplot.data$chr=='chrY' | volplot.data$chr=='chrX'),'chr'] <- 'auto'
g <- ggplot(data=volplot.data, aes(x=logFC, y=pval, color=chr))
p1 <- g + geom_point() + labs(title='', x='log2FC', y='-log10(pvalue)') + 
	scale_color_manual(values=c('black', 'red','blue'),
		name='Chromosome', breaks=c('auto','chrX', 'chrY'), labels=c('Auto','X', 'Y')) + 
	geom_abline(slope=0, intercept=-log10(sig.5)) + annotate('rect',xmin=-1, xmax=1, ymin=0, ymax=30, alpha=0.2) + theme(legend.position="none") 
g <- ggplot(data=volplot.data, aes(x=logFC, y=pval, color=chr))
p2 <- g + geom_point() + labs(title='', x='log2FC', y='-log10(pvalue)') + 
	scale_color_manual(values=c('black', 'red','blue'),
		name='Chromosome', breaks=c('auto','chrX', 'chrY'), labels=c('Auto','X', 'Y')) + 
	scale_x_continuous(limits = c(-1, 1)) + scale_y_continuous(limits = c(0, 30)) +
	geom_abline(slope=0, intercept=-log10(sig.5))
multiplot(p1, p2, cols=2)
g <- ggplot(data=volplot.data, aes(x=logFC, y=pval, color=xesc))
g + geom_point() + labs(title='', x='log2FC', y='-log10(pvalue)') + 
	geom_abline(slope=0, intercept=-log10(sig.5)) + geom_point(data=subset(volplot.data, xesc==T),aes(color=xesc), col='green')
dev.off()

library(gtools)

load('/group/stranger-lab/immvar_data/fit.joint.CD4.Robj')
joint.cd4.fit <- data.frame(eb.fit)
joint.cd4.fit$q.value <- p.adjust(joint.cd4.fit$p.value, method="fdr")
joint.cd4.fit$rank <- rank(joint.cd4.fit$p.value) / nrow(joint.cd4.fit)
joint.cd4.fit$chr <- as.character(annots[rownames(joint.cd4.fit),'chr'])
joint.cd4.sig <- subset(joint.cd4.fit, q.value<=0.05)

sig.5 <- max(joint.cd14.sig$p.value)

pdf(file='/group/stranger-lab/czysz/ImmVar/plots/immvar.cd4.fc.pdf', width=11, height=8)
#plot(density(joint.cd14.fit$coefficients), main='CD14 density of log2 fold change for all genes', xlab='log2FC', ylab='')
#plot(density(joint.cd14.sig$coefficients), main='CD14 density of log2 fold change for significant genes', xlab='log2FC', ylab='')
ggplot(data.frame(log2fc = joint.cd4.fit$coefficients), aes(x=log2fc)) + geom_density() + labs(title='', x='log2FC', y='')
ggplot(data.frame(log2fc = joint.cd4.sig$coefficients), aes(x=log2fc)) + geom_density() + labs(title='', x='log2FC', y='')

dev.off()

pdf(file='/group/stranger-lab/czysz/ImmVar/plots/cd4.sig.chr.pdf', width=11, height=8)
cd4.bar.df <- data.frame(chr=names(table(joint.cd4.sig$chr)), count=as.numeric(table(joint.cd4.sig$chr)))
cd4.bar.df$count = cd4.bar.df$count[mixedorder(cd4.bar.df$chr)]
cd4.bar.df$chr = factor(mixedsort(cd4.bar.df$chr), levels=as.character(mixedsort(cd4.bar.df$chr)))
#cd4.bar.df$count = cd4.bar.df$count[mixedorder(cd4.bar.df$chr)]
#g <- ggplot(data=cd4.bar.df, aes(x=chr, y=count, fill=chr))
#g + geom_bar(stat='identity') + labs(title='CD14 - Per-chromosome count of significant genes', x='chromosome', y='count')

cd4.male <- c(table(subset(joint.cd4.sig,coefficients>0)$chr), chrM=0) 
cd4.male <- cd4.male[mixedorder(names(cd4.male))]
cd4.female <- c(table(subset(joint.cd4.sig, coefficients<0)$chr), chrY=0, chrM=0)
cd4.female <- cd4.female[mixedorder(names(cd4.female))]
rm(cd4.bar.df)
cd4.bar.df <- data.frame(chr=names(cd4.male), male=as.numeric(cd4.male), female=as.numeric(cd4.female))
cd4.bar.df.per <- cd4.bar.df
cd4.bar.df.per$chr = factor(mixedsort(cd4.bar.df$chr), levels=as.character(mixedsort(cd4.bar.df$chr)))
g <- ggplot(data=melt(cd4.bar.df), aes(x=chr, y=value, fill=variable))
g + geom_bar(stat='identity') + labs(title='', x='chromosome', y='count') + guides(fill=guide_legend(title="Sex")) + 
	scale_fill_manual(name="Sex", breaks=c("female", "male"), labels=c("Female", "Male"), values=c('blue', 'red'))

cd4.bar.allcount <- table(joint.cd4.fit$chr)
cd4.bar.allcount <- as.numeric(table(joint.cd4.fit$chr))[mixedorder(names(table(joint.cd4.fit$chr)))]
cd4.bar.df.per$male <- cd4.bar.df.per$male / cd4.bar.allcount
cd4.bar.df.per$female <- cd4.bar.df.per$female / cd4.bar.allcount
g <- ggplot(data=melt(cd4.bar.df.per), aes(x=chr, y=value, fill=variable))
g + geom_bar(stat='identity') + labs(title='', x='chromosome', y='percent') + guides(fill=guide_legend(title='Sex')) + 
	scale_fill_manual(name="Sex", breaks=c("female", "male"), labels=c("Female", "Male"), values=c('blue', 'red'))
dev.off()
pdf(file=paste('/group/stranger-lab/czysz/ImmVar/plots/','qq_CD4.pdf',sep=''), width=12, height=8)
qqplot.data <- data.frame(obs=-log10(sort(joint.cd4.fit$p.value)),
	exp=-log10(seq(length(joint.cd4.fit$p.value))/length(joint.cd4.fit$p.value)),
	chr=as.character(joint.cd14.fit$chr[order(joint.cd14.fit$rank)]),
	row.names=row.names(joint.cd14.fit[order(joint.cd14.fit$p.value),]))
qqplot.data$chr <- as.character(qqplot.data$chr)
qqplot.data[!(qqplot.data$chr=='chrY' | qqplot.data$chr=='chrX'),'chr'] <- 'auto'

g <- ggplot(data=qqplot.data, aes(x=exp, y=obs, color=chr))
g + geom_point() + labs(title="", x='Expected -log10(p-value)', y='Observed -log10(p-value)') + 
	scale_color_manual(values=c('black','red','blue'),
		name='Chromosome', breaks=c('auto', 'chrX', 'chrY'), labels=c('Auto','X', 'Y')) + 
	geom_abline()
dev.off()

pdf(file=paste('/group/stranger-lab/czysz/ImmVar/plots/','volcano_CD4.pdf',sep=''), width=12, height=10)
qqplot.data <- data.frame(logFC=joint.cd4.fit$coefficients,
	pval=-log10(joint.cd4.fit$p.value),
	chr=joint.cd4.fit$chr,
	row.names=rownames(joint.cd4.fit),
	xesc=joint.cd4.fit$symbol%in%xesc.genes)
qqplot.data$chr <- as.character(qqplot.data$chr)
qqplot.data[!(qqplot.data$chr=='chrY' | qqplot.data$chr=='chrX'),'chr'] <- 'auto'
g <- ggplot(data=qqplot.data, aes(x=logFC, y=pval, color=chr))
p1 <- g + geom_point() + labs(title='', x='log2FC', y='-log10(pvalue)') + 
	scale_color_manual(values=c('black','red','blue'),
		name='Chromosome', breaks=c('auto','chrX', 'chrY'), labels=c('Auto','X', 'Y')) + 
	geom_abline(slope=0, intercept=-log10(sig.5[1])) + annotate('rect', xmin=-1, xmax=1, ymin=0, ymax=30, alpha=0.2) + theme(legend.position="none") 

p2 <- g + geom_point() + labs(title='', xlab='log2FC', ylab='-log10(pvalue)') + 
	scale_color_manual(values=c('black','red','blue'),
		name='Chromosome', breaks=c('auto','chrX', 'chrY'), labels=c('Auto','X', 'Y')) + 
	scale_x_continuous(limits = c(-1, 1)) + scale_y_continuous(limits = c(0, 30)) +
	geom_abline(slope=0, intercept=-log10(sig.5[1])) 
multiplot(p1, p2, cols=2)
g <- ggplot(data=qqplot.data, aes(x=logFC, y=pval, color=xesc))
g + geom_point() + labs(title='', x='log2FC', y='-log10(pvalue)') + 
	geom_abline(slope=0, intercept=-log10(sig.5)) + geom_point(data=subset(qqplot.data, xesc==T),aes(color=xesc), col='green')
dev.off()


joint.shared.sig <- intersect(rownames(joint.cd14.sig), rownames(joint.cd4.sig))

joint.venn <- list(CD4 = rownames(joint.cd4.sig), CD14=rownames(joint.cd14.sig))

###############################################
# Process meta-analysis of replication datasets
setwd('/group/stranger-lab/immvar/meta/')

cd14.rep.meta <- read.table(file='cd14_meta_immvar1.txt', header=T, row.names=1)
#cd14.rep.meta <- subset(cd14.rep.meta, Weight==2097)
cd4.rep.meta <- read.table(file='cd4_meta_immvar1.txt', header=T, row.names=1)
#cd4.rep.meta <- subset(cd4.rep.meta, Weight==719)

cd14.rep.meta <- cbind(cd14.rep.meta, 
	q.value=p.adjust(cd14.rep.meta$P.value, method='fdr'),
	rank=rank(cd14.rep.meta$P.value) / nrow(cd14.rep.meta),
	chr=annots[rownames(cd14.rep.meta), 'chr'],
	gene=annots[rownames(cd14.rep.meta), "symbol_id"])
cd14.rep.meta$chr <- as.character(cd14.rep.meta$chr)
cd14.xesc <- (cd14.rep.meta$gene%in%xesc.genes & cd14.rep.meta$chr=='chrX')
cd14.x.notxesc <- (!(cd14.rep.meta$gene%in%xesc.genes) & cd14.rep.meta$chr=='chrX')

tmp <- character(nrow(cd14.rep.meta))
tmp <- replace(tmp, cd14.xesc, 'xesc')
tmp <- replace(tmp, cd14.x.notxesc, 'nxesc')
tmp <- replace(tmp, tmp=="", 'notx')
#cd14.rep.meta$xesc <- (cd14.rep.meta$gene%in%xesc.genes & cd14.rep.meta$chr=='chrX')
cd14.rep.meta$xesc <- tmp

cd14.rep.sig.5 <- subset(cd14.rep.meta, q.value<0.05)
sig.5 <- max(cd14.rep.sig.5$P.value)

pdf(file='/group/stranger-lab/czysz/ImmVar/plots/metacd14.sig.chr.pdf', width=11, height=8)
cd14.bar.df <- data.frame(chr=names(table(cd14.rep.sig.5$chr)), count=as.numeric(table(cd14.rep.sig.5$chr)))
cd14.bar.df$count = cd14.bar.df$count[mixedorder(cd14.bar.df$chr)]
cd14.bar.df$chr = factor(mixedsort(cd14.bar.df$chr), levels=as.character(mixedsort(cd14.bar.df$chr)))
#cd14.bar.df$count = cd14.bar.df$count[mixedorder(cd14.bar.df$chr)]
g <- ggplot(data=cd14.bar.df, aes(x=chr, y=count, fill=chr))
#g + geom_bar(stat='identity') + labs(title='CD14 - Per-chromosome count of significant genes', x='chromosome', y='count')

#cd14.male <- c(table(subset(cd14.rep.sig.5,Zscore<0)$chr), chrM=0) 
cd14.male <- table(subset(cd14.rep.sig.5,Zscore<0)$chr) 
cd14.male <- cd14.male[mixedorder(names(cd14.male))]
cd14.female <- c(table(subset(cd14.rep.sig.5,Zscore>0)$chr), chrY=0)
cd14.female <- cd14.female[mixedorder(names(cd14.female))]
cd14.bar.df <- data.frame(chr=mixedsort(names(table(cd14.rep.sig.5$chr))), male=as.numeric(cd14.male), female=as.numeric(cd14.female))
cd14.bar.df$chr = factor(mixedsort(cd14.bar.df$chr), levels=as.character(mixedsort(cd14.bar.df$chr)))
library(reshape2)
g <- ggplot(data=melt(cd14.bar.df), aes(x=chr, y=value, fill=variable))
g + geom_bar(stat='identity') + labs(title='', x='chromosome', y='count') + guides(fill=guide_legend(title="Sex")) + 
	scale_fill_manual(name="Sex", breaks=c("female", "male"), labels=c("Female", "Male"), values=c('blue', 'red'))

cd14.bar.allcount <- as.numeric(table(cd14.rep.meta$chr))[mixedorder(names(table(cd14.rep.meta$chr)))]
cd14.bar.df$male <- cd14.bar.df$male / cd14.bar.allcount
cd14.bar.df$female <- cd14.bar.df$female / cd14.bar.allcount
cd14.bar.df$female[23] <- 0
#cd14.bar.df <- data.frame(chr=names(table(cd14.rep.sig.5$chr)), count=as.numeric(table(cd14.meta.sig.5$chr))/as.numeric(table(cd14.rep.meta$chr))[-23])
#cd14.bar.df$chr = factor(mixedsort(cd14.bar.df$chr), levels=as.character(mixedsort(cd4.bar.df$chr)))
#cd14.bar.df$count = cd14.bar.df$count[mixedorder(cd14.bar.df$chr)]
g <- ggplot(data=melt(cd14.bar.df), aes(x=chr, y=value, fill=variable))
g + geom_bar(stat='identity') + labs(title='', x='chromosome', y='percent') + guides(fill=guide_legend(title='Sex')) + 
	scale_fill_manual(name="Sex", breaks=c("female", "male"), labels=c("Female", "Male"), values=c('blue', 'red'))
dev.off()

pdf(file='/group/stranger-lab/czysz/ImmVar/plots/cd14.meta.qq.pdf', width=12, height=8)
g <- ggplot(data=cd14.rep.meta, aes(x=-log10(seq(nrow(cd14.rep.meta))/nrow(cd14.rep.meta)), y=-log10(sort(P.value)), color=replace(cd14.rep.meta$chr, !(cd14.rep.meta$chr=='chrX' | cd14.rep.meta$chr=='chrY'), 'auto')[order(cd14.rep.meta$rank)]))
g + geom_point() + labs(title="", x='Expected -log10(p-value)', y='Observed -log10(p-value)') + 
	scale_color_manual(values=c('black', 'red','blue'),
		name='Chromosome', breaks=c('auto', 'chrX', 'chrY'), labels=c('Auto', 'X', 'Y')) + 
	geom_abline()
dev.off()

pdf(file='/group/stranger-lab/czysz/ImmVar/plots/cd14.meta.volcano.pdf', width=12, height=10)
g <- ggplot(data=cd14.rep.meta, aes(x=Zscore, y=-log10(P.value), color=replace(cd14.rep.meta$chr, !(cd14.rep.meta$chr=='chrX' | cd14.rep.meta$chr=='chrY'), 'auto'))) 
p1 <- g + geom_point() + labs(title ='', x='Zscore', y='-log10(p-value)') +
	scale_color_manual(values=c('black', 'red','blue'),
		name='Chromosome', breaks=c('auto', 'chrX', 'chrY'), labels=c('Auto', 'X', 'Y')) + 
	geom_abline(slope=0, intercept=-log10(max(cd14.rep.sig.5$P.value))) + annotate("rect", xmin=-10, xmax=10, ymin=0, ymax=30, alpha=0.2) + theme(legend.position="none") 

p2 <- g + geom_point() + labs(title ='', x='Zscore', y='-log10(p-value)') +
	scale_color_manual(values=c('black', 'red','blue'),
		name='Chromosome', breaks=c('auto', 'chrX', 'chrY'), labels=c('Auto', 'X', 'Y')) + 
	geom_abline(slope=0, intercept=-log10(max(cd14.rep.sig.5$P.value))) +
	scale_x_continuous(limits = c(-10, 10)) + scale_y_continuous(limits = c(0, 30)) 
multiplot(p1,p2, cols=2)

g <- ggplot(data=cd14.rep.meta, aes(x=Zscore, y=-log10(P.value), color=xesc))
g + geom_point() + labs(title ='', x='Zscore', y='-log10(p-value)') + geom_abline(slope=0, intercept=max(cd14.rep.sig.5$P.value)) + scale_color_manual(name='X-inactivation escape',breaks=c('notx', 'nxesc', 'xesc'), values=c('grey50', '#E69F00', 'green'), labels=c('Non-X', 'X-linked, not XIE', 'X-linked, XIE'))
dev.off()
cd4.rep.meta <- cbind(cd4.rep.meta, 
	q.value=p.adjust(cd4.rep.meta$P.value, method='fdr'),
	rank=rank(cd4.rep.meta$P.value) / nrow(cd4.rep.meta),
	chr=annots[rownames(cd4.rep.meta), 'chr'],
	gene=annots[rownames(cd4.rep.meta), "symbol_id"])
cd4.rep.meta$chr <- as.character(cd4.rep.meta$chr)
tmp <- character(nrow(cd4.rep.meta))
cd4.xesc <- (cd4.rep.meta$gene%in%xesc.genes & cd4.rep.meta$chr=='chrX')
cd4.x.notxesc <- (!(cd4.rep.meta$gene%in%xesc.genes) & cd4.rep.meta$chr=='chrX')
tmp <- replace(tmp, cd4.xesc, 'xesc')
tmp <- replace(tmp, cd4.x.notxesc, 'nxesc')
tmp <- replace(tmp, tmp=="", 'notx')
cd4.rep.meta$xesc <- tmp
#cd4.rep.meta$xesc <- cd4.rep.meta$gene%in%xesc.genes
cd4.rep.sig.5 <- subset(cd4.rep.meta, q.value<0.05)

pdf(file='/group/stranger-lab/czysz/ImmVar/plots/cd4.meta.qq.pdf', width=12, height=8)
g <- ggplot(data=cd4.rep.meta, aes(x=-log10(seq(nrow(cd4.rep.meta))/nrow(cd4.rep.meta)), y=-log10(sort(P.value)), color=replace(cd4.rep.meta$chr, !(cd4.rep.meta$chr=='chrX' | cd4.rep.meta$chr=='chrY'), 'auto')[order(cd4.rep.meta$rank)]))
g + geom_point() + labs(title="", x='Expected -log10(p-value)', y='Observed -log10(p-value)') + 
	scale_color_manual(values=c('black', 'red','blue'),
		name='Chromosome', breaks=c('auto', 'chrX', 'chrY'), labels=c('Auto', 'X', 'Y')) + 
	geom_abline()

dev.off()

pdf(file='/group/stranger-lab/czysz/ImmVar/plots/metacd4.sig.chr.pdf', width=11, height=8)
cd4.bar.df <- data.frame(chr=names(table(cd4.rep.sig.5$chr)), count=as.numeric(table(cd4.rep.sig.5$chr)))
cd4.bar.df$count = cd4.bar.df$count[mixedorder(cd4.bar.df$chr)]
cd4.bar.df$chr = factor(mixedsort(cd4.bar.df$chr), levels=as.character(mixedsort(cd4.bar.df$chr)))
#cd4.bar.df$count = cd4.bar.df$count[mixedorder(cd4.bar.df$chr)]
g <- ggplot(data=cd4.bar.df, aes(x=chr, y=count, fill=chr))
#g + geom_bar(stat='identity') + labs(title='CD4 - Per-chromosome count of significant genes', x='chromosome', y='count')

cd4.male <- c(table(subset(cd4.rep.sig.5,Zscore<0)$chr), chrM=0) 
cd4.male <- cd4.male[mixedorder(names(cd4.male))]
cd4.female <- c(table(subset(cd4.rep.sig.5,Zscore>0)$chr), chrY=0, chrM=0)
cd4.female <- cd4.female[mixedorder(names(cd4.female))]
rm(cd4.bar.df)
cd4.bar.df <- data.frame(chr=names(cd4.male), male=as.numeric(cd4.male), female=as.numeric(cd4.female))
cd4.bar.df$chr = factor(mixedsort(cd4.bar.df$chr), levels=as.character(mixedsort(cd4.bar.df$chr)))
g <- ggplot(data=melt(cd4.bar.df), aes(x=chr, y=value, fill=variable))
g + geom_bar(stat='identity') + labs(title='', x='chromosome', y='count') + guides(fill=guide_legend(title='Sex')) + 
	scale_fill_manual(name="Sex", breaks=c("female", "male"), labels=c("Female", "Male"), values=c('blue', 'red'))

cd4.bar.allcount <- table(cd4.rep.meta$chr)
cd4.bar.allcount <- cd4.bar.allcount[mixedorder(names(cd4.bar.allcount))]
cd4.bar.df$male <- as.numeric(cd4.bar.df$male / cd4.bar.allcount)
cd4.bar.df$female <- as.numeric(cd4.bar.df$female / cd4.bar.allcount)
#cd4.bar.df$female[23] <- 0
#cd4.bar.df <- data.frame(chr=names(table(cd4.rep.sig.5$chr)), count=as.numeric(table(cd4.meta.sig.5$chr))/as.numeric(table(cd4.rep.meta$chr))[-23])
#cd4.bar.df$chr = factor(mixedsort(cd4.bar.df$chr), levels=as.character(mixedsort(cd4.bar.df$chr)))
#cd4.bar.df$count = cd4.bar.df$count[mixedorder(cd4.bar.df$chr)]
g <- ggplot(data=melt(cd4.bar.df), aes(x=chr, y=value, fill=variable))
g + geom_bar(stat='identity') + labs(title='', x='chromosome', y='percent') + guides(fill=guide_legend(title='Sex')) + 
	scale_fill_manual(name="Sex", breaks=c("female", "male"), labels=c("Female", "Male"), values=c('blue', 'red'))
dev.off()

pdf(file='/group/stranger-lab/czysz/ImmVar/plots/cd4.meta.volcano.pdf', width=12, height=10)
#g <- ggplot(data=cd4.rep.meta, aes(x=Zscore, y=-log10(P.value), color=replace(cd4.rep.meta$chr, !(cd4.rep.meta$chr=='chrX' | cd4.rep.meta$chr=='chrY'), 'auto'))) 
g <- ggplot(data=cd4.rep.meta, aes(x=Zscore, y=-log10(P.value), color=replace(chr, !(chr=='chrX' | chr=='chrY'), 'auto'))) 
p1 <- g + geom_point() + labs(title ='', x='Zscore', y='-log10(P.value)') +
	scale_color_manual(values=c('black', 'red','blue'),
		name='Chromosome', breaks=c('auto', 'chrX', 'chrY'), labels=c('Auto', 'X', 'Y')) +
	geom_abline(slope=0, intercept=max(cd14.rep.sig.5$P.value)) + annotate("rect", xmin=-10, xmax=10, ymin=0, ymax=30, alpha=0.2) + theme(legend.position="none") 

p2 <- g + geom_point() + labs(title ='', x='Zscore', y='-log10(p-value)') +
	scale_color_manual(values=c('black', 'red','blue'),
		name='Chromosome', breaks=c('auto', 'chrX', 'chrY'), labels=c('Auto', 'X', 'Y')) + 
	geom_abline(slope=0, intercept=-log10(max(cd14.rep.sig.5$P.value))) + 
	scale_x_continuous(limits = c(-10, 10)) + scale_y_continuous(limits = c(0, 30)) 
multiplot(p1, p2, cols=2)
g <- ggplot(data=cd4.rep.meta, aes(x=Zscore, y=-log10(P.value), color=xesc))
g + geom_point() + labs(title ='', x='Zscore', y='-log10(p-value)') + geom_abline(slope=0, intercept=max(cd4.rep.sig.5$P.value)) + scale_color_manual(name='X-inactivation escape',breaks=c('notx', 'nxesc', 'xesc'), values=c('grey50', '#E69F00', 'green'), labels=c('Non-X', 'X-linked, not XIE', 'X-linked, XIE'))
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

if (F) {

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

load('/group/stranger-lab/moliva/ImmVar/Robjects/CD4.joint.norm.exp_genes.Robj')
cd4.exp <- exp_genes.cd4.joint.norm
rm(exp_genes.cd4.joint.norm)
load('/group/stranger-lab/moliva/ImmVar/Robjects/CD14.joint.norm.exp_genes.Robj')
cd14.exp <- exp_genes.cd14.joint.norm
rm(exp_genes.cd14.joint.norm)

cd4.means <- rowMeans(cd4.exp) / apply(cd4.exp, 1, sd)
cd14.means <- rowMeans(cd14.exp) / apply(cd14.exp, 1, sd)

cd4.quantiles <- quantile(cd4.means, probs=seq(10)/10)
cd14.quantiles <- quantile(cd14.means, probs=seq(10)/10)

cd4.quants <- apply(as.matrix(cd4.means), 1, function(x, cd4.quantiles) {min(which(x<=cd4.quantiles))}, cd4.quantiles=cd4.quantiles)
cd14.quants <- apply(as.matrix(cd14.means), 1, function(x, cd14.quantiles) {min(which(x<=cd14.quantiles))}, cd14.quantiles=cd14.quantiles)

df <- data.frame(gene=rownames(joint.cd4.fit), quantile=cd4.quants, CD4=joint.cd4.fit$p.value, CD14=joint.cd14.fit$p.value)
cd4.mean <- mean(df$CD4)
cd4.sd <- sd(df$CD4)
cd14.mean <- mean(df$CD14)
cd14.sd <- sd(df$CD14)
df <- df[order(df$quantile, df$gene),]

cd4.p <- c()
cd14.p <- c()
for (i in seq(10)) {
	tmp <- subset(df, quantile==i)
	cd4.p <- c(cd4.p, (tmp$CD4 - cd4.mean) / cd4.sd)
	cd14.p <- c(cd14.p, (tmp$CD4 - cd14.mean) / cd14.sd)
}

#df$CD4 <- cd4.p
#df$CD14 <- cd14.p
df.melt <- melt(df, id=c('gene','quantile'), measure.vars=c('CD4','CD14'))
#df.melt <- df.melt[order(df.melt$quantile, df.melt$gene),] # Order by quantile, then by gene
# Gene quantile variable(cell type) value(p.value)

cd4.melt <- subset(df.melt, variable=='CD4')
cd14.melt <- subset(df.melt, variable=='CD14')

if (F) {
pdf(width=12, height=10)
g <- ggplot(data=df.melt, aes(x=as.factor(quantile), y=value, fill=variable)) + geom_violin(draw_quantiles=c(0.25,0.5,0.75), position=position_dodge(), scale='width') + geom_boxplot(position=position_dodge(width=1), width=0.2)
print(g)
g <- ggplot(data=df.melt, aes(x=as.factor(quantile), y=value, fill=variable)) + geom_boxplot(position=position_dodge())
print(g)
dev.off()
}

meta.shared <- intersect(rownames(cd4.rep.sig.5), rownames(cd14.rep.sig.5))
meta.union <- union(rownames(cd4.rep.sig.5), rownames(cd14.rep.sig.5))
meta.union.annot <- annots[meta.union,]
meta.shared.annot <- annots[meta.shared, ]

horm.prod <- read.table(file='/home/t.cri.cczysz/hormone.producing.txt', header=F)
horm.prod <- t(horm.prod)
horm.prod <- horm.prod[c(-1:-2),]
horm.prod <- unique(horm.prod)
horm.recpt <- read.table(file='/home/t.cri.cczysz/hormone.receptors.txt', header=F)
horm.recpt <- t(horm.recpt)
horm.recpt <- horm.recpt[c(-1:-2),]
horm.recpt <- unique(horm.recpt)
