library(lumi)
library(limma)
library(lumiHumanIDMapping)
library(illuminaHumanv3.db)
library(lumiHumanAll.db)
library(sva)
library(massiR)
library(preprocessCore)
library(oligo)
library(annotate)

setwd('/group/stranger-lab/immvar_rep/GenCord')
load('/group/stranger-lab/immvar_rep/mappings/annot.Robj')
importRawData <- function() {
	dat <- lumiR('GSE17080_unnormalized_T-cells.txt', 
	lib.mapping = 'lumiHumanIDMapping',
	columnNameGrepPattern=list(detection=NA, beadNum=NA))
	fdat <<- fData(dat)	
	dat.exp <- exprs(dat)

	if (F) {	
	pdf(file='raw_density.pdf', width=8.5, height=11)
	plot(dat, what='density')
	dev.off()

	pdf(file='raw_density.pdf', width=8.5, height=11)
	plot(dat, what='MAplot', smoothScatter=T)
	plot(example.lumi, what='sampleRelation', method='mds', cex=0.5)
	dev.off()
}


	ids<-c()
	for (id in colnames(dat.exp)) {
		ids<-c(ids, strsplit(id, split="_")[[1]][1])
	}

	# First, quantile normalize between technical replicates, taking the average as the expression for each individual
	replicates <- numeric(ncol(dat))
	for (i in seq(length(unique(ids)))) {
		replicates[which(ids==unique(ids)[i])] <- i
	}
	dat.N <- data.frame(row.names=rownames(dat))
	for (i in seq(length(unique(ids)))) {
		expr <- dat.exp[, replicates==i]
		expr.qn <- normalize.quantiles(expr)
		expr.mean <- rowMeans(expr.qn)
		dat.N <- cbind(dat.N, expr.mean)
	}
	dat.N <- as.matrix(dat.N)
	dat.N.T <- normalize.quantiles(dat.N)
	colnames(dat.N.T) <- unique(ids)
	rownames(dat.N.T) <- rownames(dat.exp)
	save(dat.N.T, file='norm_data.Robj')
	return(dat.N.T)
}

predictSex <- function(temp.dat, plot_massi=T) {
	## Verify massiR results
	y_probes <- data.frame(row.names=c("Ng17goGe3f5AsKKHCo", "oE68.e0uHwT2v8UlJI", "KuyUIoSDxODhLlSLhI", "E35Tklf.OiQtEMSSfk",
		"KLpuVfpdEpAgwfAJ4Y", "TN6X5VSSxLX1wwao50", "NnMe15_xBOv6rgqigI","T7rNJ6HSwRV.fSeqOo",
		"6VE8ScX3J.LOsufT0E", "reXjDTnq9wwlcfUXDc", "WN_KCVC5sX8._VUUV0", "Z_OQlcBHf2HUiiCekk", "QEpD_7OuProlBdrbvo"))
	temp.dat <- as.matrix(selDataMatrix)
	temp.dat["T7rNJ6HSwRV.fSeqOo",] <- -1 * log2(temp.dat["T7rNJ6HSwRV.fSeqOo",]/mean(temp.dat["T7rNJ6HSwRV.fSeqOo",]))
	massi.y.out <- massi_y(data.frame(temp.dat), y_probes)

	massi.select.out <- massi_select(data.frame(temp.dat), y_probes, threshold=1)
	results <- massi_cluster(massi.select.out)
	sample.results <- data.frame(results[[2]])

	if (plot_massi) {
	massi_y_plot(massi.y.out)
	massi_cluster_plot(massi.select.out, results)
	dev.off();dev.off()
	}
	return(sample.results)
}

fitData <- function(exprn, covs) {
	sample.factors <- as.numeric(sample.sex=='male')
	mod = model.matrix(~0+as.factor(sample.factors))
	colnames(mod) <- c('Female', 'Male')
	mod0 = model.matrix(~0+rep(1,ncol(exprn)))

	svobj <- sva(exprn, mod, mod0)
	modSv <- cbind(mod, svobj$sv)

	fit <- lmFit(exprn, modSv)
	contrast.matrix <- c(-1,1,rep(0, svobj$n.sv))
	contrast.fit <- contrasts.fit(fit, contrast.matrix)
	fit <- eBayes(contrast.fit)

	save(fit, file='gencord_fit.Robj')
	return(fit)
}

if (F) {
if (!file.exists('norm_data.Robj')) { dat.N.T <- importRawData() 
} else { load('norm_data.Robj') }
}

dat.N.T <- importRawData()
print(paste("Number of probes (raw data): ", nrow(dat.N.T), sep=''))

fdat <- data.frame(nuID=rownames(fdat), row.names=fdat$ProbeID)
# QC Plots
selDataMatrix <- dat.N.T
geneNames <- getSYMBOL(rownames(selDataMatrix), 'lumiHumanAll.db')
sample.results <- predictSex(as.matrix(selDataMatrix))
probeList <- rownames(selDataMatrix)

# Import probes passing mapping criteria
filt.probe2gene <- read.table(file='/group/stranger-lab/moliva/ImmVar/replication_datasets/mappings/HumanHT-12v3/HumanHT-12v3.probe2gene.filt.tsv', header=F, col.names=c('ProbeID', 'EnsemblID'))

# Add conversion to nuID
nuids <- c()
for (id in filt.probe2gene$ProbeID) {
	nuids<-c(nuids, as.character(fdat[id,]))
}
filt.probe2gene <- cbind(filt.probe2gene, NuID=nuids)
filt.probe2gene <- filt.probe2gene[!(is.na(filt.probe2gene$NuID)),]
save(filt.probe2gene, file='/group/stranger-lab/immvar_rep/mappings/ilmnv3.Robj')
dup.probes <- unique(filt.probe2gene[duplicated(filt.probe2gene$NuID),"NuID"])
# Remove filtered probes from expr matrix
selDataMatrix <- selDataMatrix[rownames(selDataMatrix)%in%filt.probe2gene$NuID,]

pdf('non_summarized_density.pdf')
plot(density(log2(as.numeric(selDataMatrix))), main='Density of Non-Summarized Gencord data')
dev.off()

# Save single-mapping probes
nonDupMat <- selDataMatrix[!(rownames(selDataMatrix)%in%dup.probes),]

ensemblIDs <- c()
for (id in rownames(nonDupMat)) {
	ensemblIDs <- c(ensemblIDs, as.character(filt.probe2gene[filt.probe2gene$NuID==id, "EnsemblID"]))
}

# Some probes map to multiple genes, so need one row per mapping in expr matrix
dupMat <- selDataMatrix[(rownames(selDataMatrix)%in%dup.probes),]

dupMatExp <- NULL
for (id in rownames(dupMat)) {
	ids <- filt.probe2gene[filt.probe2gene$NuID==id, "EnsemblID"]
	for (i in seq(length(ids))) {
		dupMatExp <- rbind(dupMatExp, dupMat[id, ])
	}
	ensemblIDs <- c(ensemblIDs, as.character(ids))
}

selDataMatrix <- rbind(nonDupMat, dupMatExp)
exp.summarized <- summarize(selDataMatrix, probes=ensemblIDs)

if (F) {
pdf('summarized_density.pdf')
plot(density(as.numeric(exp.summarized)), main='Density of Summarized GenCord Data')
dev.off()
}

# Reorder predicted sex to match column order in expression matrix
sample.sex <-c()
for (id in colnames(exp.summarized)) {
	sample.sex<-c(sample.sex, sample.results[sample.results$ID==id, "sex"])	
}
sample.factors <- as.numeric(sample.sex=='male')

# Filter out genes below 10th percentile of expression in both males and females
f.m <- pOverA(0.10, quantile(as.numeric(exp.summarized[, !!sample.factors]), probs=c(0.1)))
f.f <- pOverA(0.10, quantile(as.numeric(exp.summarized[, !sample.factors]), probs=c(0.1)))
ffun <- filterfun(f.m, f.f)
filtered.genes <- genefilter(exp.summarized, ffun)
save(filtered.genes, file='filter_genes.Robj')

print(paste("Genes before expression filtering:", nrow(exp.summarized), sep="\t"))
exp.summarized <- exp.summarized[filtered.genes, ]
print(paste("Genes after expression filtering:", nrow(exp.summarized), sep="\t"))
save(exp.summarized, file="gencord_tcells_summarized.Robj")

fit <- fitData(exp.summarized, sample.factors)
gene.annots <- read.table(file='/group/stranger-lab/moliva/ImmVar/probes_mapping/annotations/gencode.v22.TSS.txt', header=T, row.names=1)
gene.names <- c()
for (id in rownames(fit)) {
	gene.names <- c(gene.names, as.character(gene.annots[id, "symbol_id"]))
}
fit$genes <- gene.names
fit$chr <- annots[rownames(fit), "chr"]

plotVolcano <- function(fit) {
	pdf('volcano.pdf', width=10, height=10)
	volcanoplot(fit, names=fit$genes, cex=0.5, highlight=75)
	dev.off()
}
plotVolcano(fit)
write.table(data.frame(p.val=fit$p.value, p.adj=p.adjust(fit$p.value, method='fdr')), quote=F, row.names=T, col.names=T, file='/group/stranger-lab/immvar_rep/GenCord/fit.txt')
#write.table(topTable(fit, number=Inf, p.value=0.05), file='sig_eb_fit.txt', quote=F)
