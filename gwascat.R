if (!(require(gwascat))) {
	library(BiocInstaller)
	biocLite('gwascat', suppressUpdates=T)
}

library(gwascat)
library(ggplot2)
library(reshape2)
load('/group/stranger-lab/immvar_rep/mappings/annot.Robj')

if (file.exists(file='/group/stranger-lab/czysz/disease_genes.Robj')) {

data(gwrngs38) # load gwas data

#print(sort(unique(getTraits(gwrngs38)))) # print list of all traits in data

inflam.toi <- c("Ankylosing spondylitis", "Crohn's disease", "Ulcerative colitis", "Multiple sclerosis", "Type 1 diabetes", "Rheumatoid arthritis", "Primary biliary cirrhosis", "Systemic lupus erythematosus", "Psoriasis") # save trait of interest
metab.toi <- c('Type 2 diabetes', 'Obesity', 'Weight', "HDL cholesterol", "LDL cholesterol", "Cholesterol, total", "Body mass index")
neuro.toi <- c("Alzheimer's disease", 'Autism', 'Bipolar disorder', "Bipolar disorder and schizophrenia", "Cognitive performance", "Cognitive test performance", "Parkinson's disease", "Schizophrenia")
cancer.toi <- c("Bladder cancer", "Breast cancer", "Colorectal cancer", "Endometrial cancer", "Esophageal cancer", "Esophageal cancer and gastric cancer", "Glioma", "Lung adenocarcinoma", "Lung cancer", "Non-small cell lung cancer", "Ovarian cancer", "Pancreatic cancer", "Prostate cancer", "Thyroid cancer", "Testicular cancer", "Testicular germ cell cancer", "Testicular germ cell tumor", "Upper aerodigestive tract cancers", "Urinary bladder cancer")

getGenes <- function(toi.list) {
	gene.list <- list()
	for (i in toi.list) {	
	genes <- subsetByTraits(gwrngs38, tr=i)$Mapped_gene # pull out entries with trait of interest
	genes <- unique(genes)
	genes <- genes[(unlist(lapply(gregexpr('-', genes),'[', 1))==-1 & unlist(lapply(gregexpr(';', genes),'[',1))==-1)]
	gene.list[[i]] <- genes }
	return(gene.list)

}
if (F) {
inflam.genes <- list()
for (i in inflam.toi) {
	genes <- subsetByTraits(gwrngs38, tr=i)$Mapped_gene # pull out entries with trait of interest
	genes <- unique(genes)
	#x <- (sapply(genes, FUN=function(j) {gregexpr('-', j)[[1]]==1}))
	#genes <- genes[!(sapply(genes, FUN=function(x) {grep('[-;]', x)}))]
	genes <- genes[(unlist(lapply(gregexpr('-', genes),'[', 1))==-1 & unlist(lapply(gregexpr(';', genes),'[',1))==-1)]
	inflam.genes[[i]] <- genes
}

metab.genes <- list()
for (i in metab.toi) {
	metab.genes[[i]] <- subsetByTraits(gwrngs38, tr=i)$Mapped_gene # pull out entries with trait of interest
}
neuro.genes <- list()
for (i in neuro.toi) {
	neuro.genes[[i]] <- subsetByTraits(gwrngs38, tr=i)$Mapped_gene # pull out entries with trait of interest
}
cancer.genes <- list()
for (i in cancer.toi) {
	cancer.genes[[i]] <- subsetByTraits(gwrngs38, tr=i)$Mapped_gene # pull out entries with trait of interest
}
}
inflam.genes <- getGenes(inflam.toi)
metab.genes <- getGenes(metab.toi)
neuro.genes <- getGenes(neuro.toi)
cancer.genes <- getGenes(cancer.toi)

disease.genes = list(Inflammatory=inflam.genes,
	Metabolism=metab.genes, 
	Neurological=neuro.genes, 
	Cancer=cancer.genes)

save(disease.genes, file='/group/stranger-lab/czysz/disease_genes.Robj') } else {load(file='/group/stranger-lab/czysz/disease_genes.Robj') }

disease.gene.pairs <- melt(unlist(disease.genes))
library(limma)
library(reshape2)
load('/group/stranger-lab/immvar_data/fit.joint.CD14.Robj')
joint.cd14.fit <- data.frame(eb.fit)
joint.cd14.fit$q.value <- p.adjust(joint.cd14.fit$p.value, method="fdr")
joint.cd14.fit$rank <- rank(joint.cd14.fit$p.value) / nrow(joint.cd14.fit)
joint.cd14.fit$chr <- as.character(annots[rownames(joint.cd14.fit),'chr'])
joint.cd14.sig <- subset(joint.cd14.fit, q.value<=0.05)

cd14.disease.gene.counts <- lapply(disease.genes, FUN=function(i) {lapply(i, FUN=function(x, joint.cd14.fit) {sum(as.character(joint.cd14.fit$symbol)%in%x)}, joint.cd14.fit)})
cd14.sig.disease.gene.counts <- lapply(disease.genes, FUN=function(i) {lapply(i, FUN=function(x, joint.cd14.sig) {sum(as.character(joint.cd14.sig$symbol)%in%x)}, joint.cd14.sig)})
sig.5 <- max(joint.cd14.sig$p.value)

cd14.disease.genes <- rep('none', nrow(joint.cd14.fit))
cd14.immune <- joint.cd14.fit$symbol%in%unique(unlist(disease.genes[[1]]))
cd14.metab <- joint.cd14.fit$symbol%in%unique(unlist(disease.genes[[2]]))
cd14.neuro <- joint.cd14.fit$symbol%in%unique(unlist(disease.genes[[3]]))
cd14.cancer <- joint.cd14.fit$symbol%in%unique(unlist(disease.genes[[4]]))
cd14.disease.genes <- replace(cd14.disease.genes, (cd14.immune | cd14.metab | cd14.neuro | cd14.cancer), 'Disease')

pdf(file=paste('/group/stranger-lab/czysz/ImmVar/plots/','gwas_volcano_CD14.pdf',sep=''))
volplot.data <- data.frame(logFC=joint.cd14.fit$coefficients,
        pval=-log10(joint.cd14.fit$p.value),
        chr=as.character(joint.cd14.fit$chr),
        row.names=rownames(joint.cd14.fit),
	disease=cd14.disease.genes)
volplot.data$chr <- as.character(volplot.data$chr)
volplot.data[!(volplot.data$chr=='chrY' | volplot.data$chr=='chrX'),'chr'] <- 'auto'
g <- ggplot(data=volplot.data, aes(x=logFC, y=pval, color=chr))
g <- g + geom_point() + labs(title='ImmVar CD14', x='log2FC', y='-log10(pvalue)') +
        scale_color_manual(values=c('black', 'red','blue'),
                name='Chromosome', breaks=c('auto','chrX', 'chrY'), labels=c('Auto','X', 'Y')) +
        geom_abline(slope=0, intercept=-log10(sig.5)) + geom_point(data=subset(volplot.data, !(disease=="none")), col='green')
g
g <- ggplot(data=volplot.data, aes(x=logFC, y=pval, color=chr))
if (T) {
g <- g + geom_point() + labs(title='ImmVar CD14-Zoomed', x='log2FC', y='-log10(pvalue)') +
        scale_color_manual(values=c('black', 'red','blue'),
                name='Chromosome', breaks=c('auto','chrX', 'chrY'), labels=c('Auto','X', 'Y')) +
        scale_x_continuous(limits = c(-1, 1)) + scale_y_continuous(limits = c(0, 30)) +
        geom_abline(slope=0, intercept=-log10(sig.5))+ geom_point(data=subset(volplot.data, !(disease=='none')), col='green')
}
g
dev.off()

load('/group/stranger-lab/immvar_data/fit.joint.CD4.Robj')
print(topTable(eb.fit))
joint.cd4.fit <- data.frame(eb.fit)
joint.cd4.fit$q.value <- p.adjust(joint.cd4.fit$p.value, method="fdr")
joint.cd4.fit$rank <- rank(joint.cd4.fit$p.value) / nrow(joint.cd4.fit)
joint.cd4.fit$chr <- as.character(annots[rownames(joint.cd4.fit),'chr'])
joint.cd4.sig <- subset(joint.cd4.fit, q.value<=0.05)

cd4.disease.genes <- rep('none', nrow(joint.cd4.fit))
cd4.immune <- joint.cd4.fit$symbol%in%unique(unlist(disease.genes[[1]]))
cd4.metab <- joint.cd4.fit$symbol%in%unique(unlist(disease.genes[[2]]))
cd4.neuro <- joint.cd4.fit$symbol%in%unique(unlist(disease.genes[[3]]))
cd4.cancer <- joint.cd4.fit$symbol%in%unique(unlist(disease.genes[[4]]))
cd4.disease.genes <- replace(cd4.disease.genes, (cd4.immune | cd4.metab | cd4.neuro | cd4.cancer), 'Disease')

cd4.disease.gene.counts <- lapply(disease.genes, FUN=function(i) {lapply(i, FUN=function(x, joint.cd4.fit) {sum(as.character(joint.cd4.fit$symbol)%in%x)}, joint.cd4.fit)})
cd4.sig.disease.gene.counts <- lapply(disease.genes, FUN=function(i) {lapply(i, FUN=function(x, joint.cd4.sig) {sum(as.character(joint.cd4.sig$symbol)%in%x)}, joint.cd4.sig)})
sig.5 <- max(joint.cd14.sig$p.value)

pdf(file=paste('/group/stranger-lab/czysz/ImmVar/plots/','gwas_volcano_CD4.pdf',sep=''))
volplot.data <- data.frame(logFC=joint.cd4.fit$coefficients,
        pval=-log10(joint.cd4.fit$p.value),
        chr=joint.cd4.fit$chr,
        row.names=rownames(joint.cd4.fit),
	disease=cd4.disease.genes)
volplot.data$chr <- as.character(volplot.data$chr)
volplot.data[!(volplot.data$chr=='chrY' | volplot.data$chr=='chrX'),'chr'] <- 'auto'
g <- ggplot(data=volplot.data, aes(x=logFC, y=pval, color=chr))
g <- g + geom_point() + labs(title='ImmVar CD4', x='log2FC', y='-log10(pvalue)') +
        scale_color_manual(values=c('black','red','blue'),
                name='Chromosome', breaks=c('auto','chrX', 'chrY'), labels=c('Auto','X', 'Y')) +
        geom_abline(slope=0, intercept=-log10(sig.5[1])) + geom_point(data=subset(volplot.data, !(disease=='none')), col='green')
g
g <- ggplot(data=volplot.data, aes(x=logFC, y=pval, color=chr))
g <- g + geom_point() + labs(title='ImmVar CD4-Zoomed', xlab='log2FC', ylab='-log10(pvalue)') +
        scale_color_manual(values=c('black','red','blue'),
                name='Chromosome', breaks=c('auto','chrX', 'chrY'), labels=c('Auto','X', 'Y')) +
        scale_x_continuous(limits = c(-1, 1)) + scale_y_continuous(limits = c(0, 30)) +
        geom_abline(slope=0, intercept=-log10(sig.5[1]))+ geom_point(data=subset(volplot.data, !(disease=='none')), col='green')
g
dev.off()

setwd('/group/stranger-lab/immvar/meta/')

cd14.rep.meta <- read.table(file='cd14_meta_weight1.txt', header=T, row.names=1)
cd14.rep.meta <- subset(cd14.rep.meta, Weight==2097)
cd4.rep.meta <- read.table(file='cd4_meta_weight1.txt', header=T, row.names=1)
cd4.rep.meta <- subset(cd4.rep.meta, Weight==719)

setwd('/group/stranger-lab/immvar')

cd14.rep.meta <- cbind(cd14.rep.meta,
        q.value=p.adjust(cd14.rep.meta$P.value, method='fdr'),
        rank=rank(cd14.rep.meta$P.value) / nrow(cd14.rep.meta),
        chr=annots[rownames(cd14.rep.meta), 'chr'],
        gene=annots[rownames(cd14.rep.meta), "symbol_id"])
cd14.rep.meta$chr <- as.character(cd14.rep.meta$chr)
cd14.rep.sig.5 <- subset(cd14.rep.meta, q.value<0.05)


annotateGenes <- function(x) {
	rownames(disease.gene.pairs)[disease.gene.pairs==x]
}

cd14.meta.disease.genes <- rep('none', nrow(cd14.rep.meta))
cd14.meta.immune <- cd14.rep.meta$gene%in%unique(unlist(disease.genes[[1]]))
cd14.meta.metab <- cd14.rep.meta$gene%in%unique(unlist(disease.genes[[2]]))
cd14.meta.neuro <- cd14.rep.meta$gene%in%unique(unlist(disease.genes[[3]]))
cd14.meta.cancer <- cd14.rep.meta$gene%in%unique(unlist(disease.genes[[4]]))
cd14.meta.disease.genes <- replace(cd14.meta.disease.genes, (cd14.meta.immune), 'Inflammatory')
cd14.meta.disease.genes <- replace(cd14.meta.disease.genes, (cd14.meta.metab), 'Metabolism')
cd14.meta.disease.genes <- replace(cd14.meta.disease.genes, (cd14.meta.neuro), 'Neurological')
cd14.meta.disease.genes <- replace(cd14.meta.disease.genes, (cd14.meta.cancer), 'Cancer')

volplot.data = cbind(cd14.rep.meta,disease=cd14.meta.disease.genes)
pdf(file='/group/stranger-lab/czysz/ImmVar/plots/cd14_gwas_meta_volcano.pdf')
g <- ggplot(data=cbind(cd14.rep.meta,disease=cd14.meta.disease.genes), 
	aes(x=-1*Zscore, y=-log10(P.value), 
	color=replace(cd14.rep.meta$chr, !(cd14.rep.meta$chr=='chrX' | cd14.rep.meta$chr=='chrY'), 'auto')))
g <- g + geom_point() + labs(title='CD14 All Meta', x='Zscore', y='Pvalue') +
        scale_color_manual(values=c('black', 'red','blue'),
                name='Chromosome', breaks=c('auto', 'chrX', 'chrY'), labels=c('Auto', 'X', 'Y')) +
        geom_abline(slope=0, intercept=-log10(max(cd14.rep.sig.5$P.value))) + geom_point(data=subset(volplot.data, !(disease=='none')), col='green')
print(g)
g <- ggplot(data=cbind(cd14.rep.meta,disease=cd14.meta.disease.genes), 
	aes(x=Zscore, y=-log10(P.value), 
	color=replace(cd14.rep.meta$chr, !(cd14.rep.meta$chr=='chrX' | cd14.rep.meta$chr=='chrY'), 'auto')))

g <- g + geom_point() + labs(title='CD14 All Meta', x='Zscore', y='Pvalue') +
        scale_color_manual(values=c('black', 'red','blue'),
                name='Chromosome', breaks=c('auto', 'chrX', 'chrY'), labels=c('Auto', 'X', 'Y')) +
        scale_x_continuous(limits = c(-20, 20)) + scale_y_continuous(limits = c(0, 30)) +
        geom_abline(slope=0, intercept=-log10(max(cd14.rep.sig.5$P.value))) + geom_point(data=subset(volplot.data, !(disease=='none')), col='green')
print(g);dev.off()
cd14.meta.gwas <- subset(volplot.data, !(disease=='none'))
print(head(cd14.meta.gwas[order(cd14.meta.gwas$P.value),c(-1:-3,-7,-8)],50))

cd14.pval <- rep(0, nrow(disease.gene.pairs))
cd14.qval <- rep(0, nrow(disease.gene.pairs))
cd14.zscore <- rep(0, nrow(disease.gene.pairs))
for (g in seq(nrow(disease.gene.pairs))) {
	#entry <- cd14.rep.meta[cd14.rep.meta$gene==disease.gene.pairs[g], ]
	entry <- subset(cd14.rep.meta, gene==disease.gene.pairs[g,1])
	if (nrow(entry) > 0) {
	cd14.pval[g] <- signif(entry$P.value, 3)
	cd14.qval[g] <- signif(entry$q.value, 3)
	cd14.zscore[g] <- entry$Zscore } else {
	cd14.pval[g] <- NA
	cd14.qval[g] <- NA
	cd14.zscore[g] <- NA }
}

disease.genes.cd14 <- disease.gene.pairs
disease.genes.cd14 <- cbind(disease.genes.cd14, pval=cd14.pval, qval=cd14.qval, zscore=cd14.zscore)
disease.genes.cd14 <- disease.genes.cd14[order(disease.genes.cd14$pval),]
cd4.rep.meta <- cbind(cd4.rep.meta,
        q.value=p.adjust(cd4.rep.meta$P.value, method='fdr'),
        rank=rank(cd4.rep.meta$P.value) / nrow(cd4.rep.meta),
        chr=annots[rownames(cd4.rep.meta), 'chr'],
        gene=annots[rownames(cd4.rep.meta), "symbol_id"])
cd4.rep.meta$chr <- as.character(cd4.rep.meta$chr)
cd4.rep.sig.5 <- subset(cd4.rep.meta, q.value<0.05)
cd4.meta.disease.genes <- rep('none', nrow(cd4.rep.meta))

cd4.meta.immune <- cd4.rep.meta$gene%in%unique(unlist(disease.genes[[1]]))
cd4.meta.metab <- cd4.rep.meta$gene%in%unique(unlist(disease.genes[[2]]))
cd4.meta.neuro <- cd4.rep.meta$gene%in%unique(unlist(disease.genes[[3]]))
cd4.meta.cancer <- cd4.rep.meta$gene%in%unique(unlist(disease.genes[[4]]))
cd4.meta.disease.genes <- replace(cd4.meta.disease.genes, (cd4.meta.immune), 'Inflammatory')
cd4.meta.disease.genes <- replace(cd4.meta.disease.genes, (cd4.meta.metab), 'Metabolism')
cd4.meta.disease.genes <- replace(cd4.meta.disease.genes, (cd4.meta.neuro), 'Neurological')
cd4.meta.disease.genes <- replace(cd4.meta.disease.genes, (cd4.meta.cancer), 'Cancer')
pdf(file='/group/stranger-lab/czysz/ImmVar/plots/cd4_gwas_meta_volcano.pdf')

volplot.data = cbind(cd4.rep.meta,disease=cd4.meta.disease.genes)
g <- ggplot(data=cbind(cd4.rep.meta,disease=cd4.meta.disease.genes), 
	aes(x=-1*Zscore, y=-log10(P.value), 
	color=replace(cd4.rep.meta$chr, !(cd4.rep.meta$chr=='chrX' | cd4.rep.meta$chr=='chrY'), 'auto')))
g <- g + geom_point() + labs(title='CD4 All Meta', x='Zscore', y='Pvalue') +
        scale_color_manual(values=c('black', 'red','blue'),
                name='Chromosome', breaks=c('auto', 'chrX', 'chrY'), labels=c('Auto', 'X', 'Y')) +
        geom_abline(slope=0, intercept=-log10(max(cd4.rep.sig.5$P.value))) + geom_point(data=subset(volplot.data, !(disease=='none')), col='green')
print(g)
g <- ggplot(data=cbind(cd4.rep.meta,disease=cd4.meta.disease.genes), 
	aes(x=Zscore, y=-log10(P.value), 
	color=replace(cd4.rep.meta$chr, !(cd4.rep.meta$chr=='chrX' | cd4.rep.meta$chr=='chrY'), 'auto')))

g <- g + geom_point() + labs(title='CD4 All Meta', x='Zscore', y='Pvalue') +
        scale_color_manual(values=c('black', 'red','blue'),
                name='Chromosome', breaks=c('auto', 'chrX', 'chrY'), labels=c('Auto', 'X', 'Y')) +
        scale_x_continuous(limits = c(-20, 20)) + scale_y_continuous(limits = c(0, 30)) +
        geom_abline(slope=0, intercept=-log10(max(cd4.rep.sig.5$P.value))) + geom_point(data=subset(volplot.data, !(disease=='none')), col='green')
print(g);dev.off()

cd4.pval <- rep(0, nrow(disease.gene.pairs))
cd4.qval <- rep(0, nrow(disease.gene.pairs))
cd4.zscore <- rep(0, nrow(disease.gene.pairs))
for (g in seq(nrow(disease.gene.pairs))) {
	#entry <- cd14.rep.meta[cd14.rep.meta$gene==disease.gene.pairs[g], ]
	entry <- subset(cd4.rep.meta, gene==disease.gene.pairs[g,1])
	if (nrow(entry) > 0) {
	cd4.pval[g] <- signif(entry$P.value, 3)
	cd4.qval[g] <- signif(entry$q.value, 3)
	cd4.zscore[g] <- entry$Zscore } else {
	cd4.pval[g] <- NA
	cd4.qval[g] <- NA
	cd4.zscore[g] <- NA }
}

disease.genes.cd4 <- disease.gene.pairs
disease.genes.cd4 <- cbind(disease.genes.cd4, pval=cd4.pval, qval=cd4.qval, zscore=cd4.zscore)
disease.genes.cd4 <- disease.genes.cd4[order(disease.genes.cd4$pval),]
cd4.meta.gwas <- subset(volplot.data, !(disease=='none'))
print(head(cd14.meta.gwas[order(cd14.meta.gwas$P.value),c(-1,-2,-3,-7,-8)],50))

shared.disease.genes <- intersect(subset(disease.genes.cd4, qval<0.05)$value, subset(disease.genes.cd14, qval<0.05)$value)
y <- subset(disease.genes.cd14, value%in%shared.disease.genes)
x <- subset(disease.genes.cd4, value%in%shared.disease.genes)
