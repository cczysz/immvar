library(lumi)
library(limma)
library(lumiHumanIDMapping)
library(illuminaHumanv4.db)
library(lumiHumanAll.db)
library(sva)
library(massiR)
library(oligo)

setwd('/group/stranger-lab/immvar_rep/GSE56045')
# Set to true if Robj needs to be replaced
if (!file.exists('norm_data.Robj')) {
	dat <- lumiR('raw.txt', 
		lib.mapping = 'lumiHumanIDMapping',
		columnNameGrepPattern=list(detection="detectionPval", exprs="intensity"))

	dat.N.T <- lumiExpresso(dat, bg.correct=F,
		varianceStabilize.param = list(method='log2'), 
		normalize.param=list(method='quantile'))
	save(dat.N.T, file='norm_data.Robj', compress="bzip2")
} else load('norm_data.Robj')

# QC Plots
if (F) {
	pdf('raw_density.pdf')
	plot(dat, what='density')
	dev.off()

	pdf('sample_relation.pdf', width=12, height=8)
	plot(dat, what='sampleRelation', cex=0.5)
	plot(dat.N.T, what='sampleRelation', cex=0.5)
	plot(dat, what='sampleRelation', method='mds')
	dev.off()
}

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

sample.results <- predictSex(as.matrix(selDataMatrix))

selDataMatrix <- as.matrix(selDataMatrix)

# Filter lowly expressed genes
f1 <- pOverA(p=0.3, A=quantile(as.numeric(selDataMatrix), probes=c(0.1)))
ffun <- filterfun(f1)
filtered.genes <- genefilter(selDataMatrix, ffun)
selDataMatrix <- selDataMatrix[filtered.genes, ]

geneNames <- getSYMBOL(rownames(selDataMatrix), 'lumiHumanAll.db')
exp.summarized <- basicRMA(selDataMatrix, geneNames, F, F)
sample.factors <- as.numeric(sample.results$sex=="male")

fit <- fitData(exp.summarized, sample.factors)
save(fit, file='gencord_fit.Robj')
top_diff <- topTable(fit, number=Inf, p.value=0.05)
top_diff$gene <- geneSymbol[rownames(top_diff)]

plotVolcano <- function(fit) {
	pdf('volcano.pdf', width=10, height=10)
	volcanoplot(fit, names=rownames(fit), cex=0.5, highlight=75, main="Sex DE - MESA T-cells")
	dev.off()
}
plotVolcano(fit)
write.table(topTable(fit, number=Inf), file='eb_fit.txt', quote=F)
