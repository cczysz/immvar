library(lumi)
library(limma)
library(lumiHumanIDMapping)
library(illuminaHumanv4.db)
library(lumiHumanAll.db)
library(sva)
library(massiR)
library(oligo)
library(genefilter)

setwd('/group/stranger-lab/immvar_rep/GSE56580')
load('/group/stranger-lab/immvar_rep/mappings/annot.Robj')
if (!file.exists('norm_data.Robj')) {
	dat <- lumiR('non_norm.txt', 
		lib.mapping = 'lumiHumanIDMapping',
		columnNameGrepPattern=list(detection="detectionPval", exprs="intensity"))

	dat.N.T <- lumiExpresso(dat, bg.correct=F, 
		variance.stabilize=F,
		varianceStabilize.param = list(method='log2'), 
		normalize.param=list(method='quantile'))
	save(dat.N.T, file='norm_data.Robj')
} else load('norm_data.Robj')

# QC Plots
if (T) {

	pdf('raw_QC.pdf', width=12, height=8)
	plot(dat, what='density')
	plot(dat, what='sampleRelation', method='mds')
	dev.off()
}

pdf('norm_QC.pdf', width=12, height=8)
plot(dat, what='density')
plot(dat, what='sampleRelation', method='mds')
dev.off()

selDataMatrix <- exprs(dat.N.T)
probeList <- rownames(selDataMatrix)
if (require(lumiHumanAll.db) & require(annotate)) {
	selDataMatrix <- selDataMatrix[!is.na(getSYMBOL(rownames(dat.N.T), 'lumiHumanAll.db')), ]
	geneSymbol <- getSYMBOL(probeList, 'lumiHumanAll.db')
	}

## massiR
# Hard coded Y probes show bimodality between predicted-male and predicted-female samples
predictSex <- function(eset, massi_plot=F) {
	best_y_probes <- data.frame(row.names=c("Ng17goGe3f5AsKKHCo", "oE68.e0uHwT2v8UlJI", "KuyUIoSDxODhLlSLhI", "E35Tklf.OiQtEMSSfk",
		"KLpuVfpdEpAgwfAJ4Y", "TN6X5VSSxLX1wwao50", "NnMe15_xBOv6rgqigI","T7rNJ6HSwRV.fSeqOo", 
		"6VE8ScX3J.LOsufT0E", "reXjDTnq9wwlcfUXDc", "WN_KCVC5sX8._VUUV0", "Z_OQlcBHf2HUiiCekk", "QEpD_7OuProlBdrbvo"))
	# Flip the expression of the XIST probe
	eset["T7rNJ6HSwRV.fSeqOo",] <- -1 * log2(eset["T7rNJ6HSwRV.fSeqOo",]/mean(eset["T7rNJ6HSwRV.fSeqOo",]))
	eset <- data.frame(eset)
	massi.y.out <- massi_y(eset, best_y_probes)

	massi.select.out <- massi_select(eset, best_y_probes, threshold=1)
	results <- massi_cluster(massi.select.out)
	sample.results <- data.frame(results[[2]])

	if (massi_plot) {
	massi_y_plot(massi.y.out)
	massi_cluster_plot(massi.select.out, results)
	dev.off();dev.off()
	}
	return(sample.results)
}

fitData <- function(exprn, covs) {
	mod = model.matrix(~0+as.factor(sample.factors))
	colnames(mod) <- c('Female', 'Male')
	mod0 = model.matrix(~0+rep(1,ncol(dat.N.T)))

	svobj <- sva(exp.summarized, mod, mod0)
	modSv <- cbind(mod, svobj$sv)

	fit <- lmFit(exp.summarized, modSv)
	contrast.matrix <- c(-1,1,rep(0, svobj$n.sv))
	contrast.fit <- contrasts.fit(fit, contrast.matrix)
	fit <- eBayes(contrast.fit)
	return(fit)
}

fdat <- fData(dat.N.T)

sample.results <- predictSex(as.matrix(selDataMatrix))

selDataMatrix <- as.matrix(selDataMatrix)

nuIDs <- rownames(selDataMatrix)
geneNames <- getSYMBOL(rownames(selDataMatrix), 'lumiHumanAll.db')

filt.probes <- read.table(file='/group/stranger-lab/moliva/ImmVar/replication_datasets/mappings/HumanHT-12v4/HumanHT-12v4.probe2gene.filt.tsv', header=F, col.names=c('ProbeID', 'EnsemblID'))
x <- c()
for (id in filt.probes$ProbeID) {
	nuid <- as.character(rownames(fdat)[fdat$ID_REF==id])
	if (length(nuid) == 0) { x <- c(x, NA) 

	} else {x <- c(x, as.character(rownames(fdat)[fdat$ID_REF==id])) }
}
filt.probes$NuID <- x
rm(x)

geneNames <- getSYMBOL(rownames(selDataMatrix), 'lumiHumanAll.db')
#selDataMatrix <- selDataMatrix[rownames(selDataMatrix)%in%filt.probes$NuID,]

tempMat <- NULL
ensIDs <- c()
for (gene in unique(filt.probes$EnsemblID)) {
	info <- filt.probes[filt.probes$EnsemblID==gene,]
	tempMat <- rbind(tempMat, selDataMatrix[rownames(selDataMatrix)%in%info$NuID,])
	ensIDs <- c(ensIDs, rep(gene, sum(rownames(selDataMatrix)%in%info$NuID)))
}
exp.summarized <- summarize(tempMat, probes=ensIDs)
sample.factors <- as.numeric(sample.results$sex=="male")

# Filter lowly expressed genes
f.m <- pOverA(0.10, quantile(as.numeric(exp.summarized[, !!sample.factors]), probs=c(0.1)))
f.f <- pOverA(0.10, quantile(as.numeric(exp.summarized[, !sample.factors]), probs=c(0.1)))
ffun <- filterfun(f.m, f.f)
filtered.genes <- genefilter(exp.summarized, ffun)
save(filtered.genes, file='filter_genes.Robj')
print(paste("Genes before expression filtering:", nrow(exp.summarized), sep="\n"))
exp.summarized <- exp.summarized[filtered.genes, ]
print(paste("Genes after expression filtering:", nrow(exp.summarized), sep="\n"))
save(exp.summarized, file="mesa_tcells_summarized.Robj")


fit <- fitData(exp.summarized, sample.factors)

annot <- read.table(file='/group/stranger-lab/moliva/ImmVar/probes_mapping/annotations/gencode.v22.TSS.txt', header=T, row.names=1)
gene.names <- c()
for (id in rownames(fit)) {
	gene.names <- c(gene.names, as.character(annot[id, "symbol_id"]))
}
fit$genes <- gene.names
save(fit, file='mesa_tcells_fit.Robj')
top_diff <- topTable(fit, number=Inf, p.value=0.05)
top_diff$gene <- geneSymbol[rownames(top_diff)]

plotVolcano <- function(fit) {
	pdf('volcano.pdf', width=10, height=10)
	volcanoplot(fit, names=fit$genes, cex=0.5, highlight=75, main="Sex DE - MESA T-cells")
	dev.off()
}
plotVolcano(fit)
write.table(topTable(fit, number=Inf), file='eb_fit.txt', quote=F)
