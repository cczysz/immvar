if (!(require(gwascat))) {
	library(BiocInstaller)
	biocLite('gwascat', suppressUpdates=T)
}

library(gwascat)
library(ggplot2)
load('/group/stranger-lab/immvar_rep/mappings/annot.Robj')

if (!file.exists(file='/group/stranger-lab/czysz/disease_genes.Robj')) {

data(gwrngs38) # load gwas data

#print(sort(unique(getTraits(gwrngs38)))) # print list of all traits in data

inflam.toi <- c("Ankylosing spondylitis", "Crohn's disease", "Ulcerative colitis", "Multiple sclerosis", "Type 1 diabetes", "Rheumatoid arthritis", "Primary biliary cirrhosis", "Systemic lupus erythematosus", "Psoriasis") # save trait of interest
metab.toi <- c('Type 2 diabetes', 'Obesity', 'Weight', "HDL cholesterol", "LDL cholesterol", "Cholesterol, total", "Body mass index")
neuro.toi <- c("Alzheimer's disease", 'Autism', 'Bipolar disorder', "Bipolar disorder and schizophrenia", "Cognitive performance", "Cognitive test performance", "Parkinson's disease", "Schizophrenia")
cancer.toi <- c("Bladder cancer", "Breast cancer", "Colorectal cancer", "Prostate cancer", "Thyroid cancer", "Testicular cancer", "Testicular germ cell cancer", "Testicular germ cell tumor")

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
library(limma)
load('/group/stranger-lab/immvar_data/fit.joint.CD14.Robj')
joint.cd14.fit <- data.frame(eb.fit)
joint.cd14.fit$q.value <- p.adjust(joint.cd14.fit$p.value, method="fdr")
joint.cd14.fit$rank <- rank(joint.cd14.fit$p.value) / nrow(joint.cd14.fit)
joint.cd14.fit$chr <- as.character(annots[rownames(joint.cd14.fit),'chr'])
joint.cd14.sig <- subset(joint.cd14.fit, q.value<=0.05)

cd14.disease.gene.counts <- lapply(y, FUN=function(i) {lapply(i, FUN=function(x, joint.cd14.fit) {sum(as.character(joint.cd14.fit$symbol)%in%x)}, joint.cd14.fit)})
cd14.sig.disease.gene.counts <- lapply(y, FUN=function(i) {lapply(i, FUN=function(x, joint.cd14.sig) {sum(as.character(joint.cd14.sig$symbol)%in%x)}, joint.cd14.sig)})
sig.5 <- max(joint.cd14.sig$p.value)

#pdf(file=paste('/group/stranger-lab/czysz/ImmVar/plots/','volcano_CD14.pdf',sep=''))
volplot.data <- data.frame(logFC=joint.cd14.fit$coefficients,
        pval=-log10(joint.cd14.fit$p.value),
        chr=as.character(joint.cd14.fit$chr),
        row.names=rownames(joint.cd14.fit))
volplot.data$chr <- as.character(volplot.data$chr)
volplot.data[!(volplot.data$chr=='chrY' | volplot.data$chr=='chrX'),'chr'] <- 'auto'
g <- ggplot(data=volplot.data, aes(x=logFC, y=pval, color=chr))
g <- g + geom_point() + labs(title='ImmVar CD14', x='log2FC', y='-log10(pvalue)') +
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
#dev.off()

load('/group/stranger-lab/immvar_data/fit.joint.CD4.Robj')
print(topTable(eb.fit))
joint.cd4.fit <- data.frame(eb.fit)
joint.cd4.fit$q.value <- p.adjust(joint.cd4.fit$p.value, method="fdr")
joint.cd4.fit$rank <- rank(joint.cd4.fit$p.value) / nrow(joint.cd4.fit)
joint.cd4.fit$chr <- as.character(annots[rownames(joint.cd4.fit),'chr'])
joint.cd4.sig <- subset(joint.cd4.fit, q.value<=0.05)

cd4.disease.gene.counts <- lapply(y, FUN=function(i) {lapply(i, FUN=function(x, joint.cd4.fit) {sum(as.character(joint.cd4.fit$symbol)%in%x)}, joint.cd4.fit)})
cd4.sig.disease.gene.counts <- lapply(y, FUN=function(i) {lapply(i, FUN=function(x, joint.cd4.sig) {sum(as.character(joint.cd4.sig$symbol)%in%x)}, joint.cd4.sig)})
sig.5 <- max(joint.cd14.sig$p.value)

#pdf(file=paste('/group/stranger-lab/czysz/ImmVar/plots/','volcano_CD4.pdf',sep=''))
qqplot.data <- data.frame(logFC=joint.cd4.fit$coefficients,
        pval=-log10(joint.cd4.fit$p.value),
        chr=joint.cd4.fit$chr,
        row.names=rownames(joint.cd4.fit))
qqplot.data$chr <- as.character(qqplot.data$chr)
qqplot.data[!(qqplot.data$chr=='chrY' | qqplot.data$chr=='chrX'),'chr'] <- 'auto'
g <- ggplot(data=qqplot.data, aes(x=logFC, y=pval, color=chr))
g <- g + geom_point() + labs(title='ImmVar CD4', x='log2FC', y='-log10(pvalue)') +
        scale_color_manual(values=c('black','red','blue'),
                name='Chromosome', breaks=c('auto','chrX', 'chrY'), labels=c('Auto','X', 'Y')) +
        geom_abline(slope=0, intercept=-log10(sig.5[1]))

g <- g + geom_point() + labs(title='ImmVar CD4-Zoomed', xlab='log2FC', ylab='-log10(pvalue)') +
        scale_color_manual(values=c('black','red','blue'),
                name='Chromosome', breaks=c('auto','chrX', 'chrY'), labels=c('Auto','X', 'Y')) +
        scale_x_continuous(limits = c(-1, 1)) + scale_y_continuous(limits = c(0, 30)) +
        geom_abline(slope=0, intercept=-log10(sig.5[1]))
#dev.off()
