library(lumi)
library(limma)
library(lumiHumanIDMapping)
library(illuminaHumanv3.db)
library(lumiHumanAll.db)
library(sva)
library(massiR)
library(preprocessCore)
library(oligo)

setwd('/group/stranger-lab/immvar_rep/GenCord')
importRawData <- function() {
	dat <- lumiR('GSE17080_unnormalized_T-cells.txt', 
	lib.mapping = 'lumiHumanIDMapping',
	columnNameGrepPattern=list(detection=NA, beadNum=NA))
	dat.exp <- exprs(dat)

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
	dat.N.T <- log2(normalize.quantiles(dat.N))
	colnames(dat.N.T) <- unique(ids)
	rownames(dat.N.T) <- rownames(dat.exp)
	save(dat.N.T, file='norm_data.Robj')
	return(dat.N.T)
}

predictSex <- function(temp.dat, plot_massi=F) {
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

if (!file.exists('norm_data.Robj')) { dat.N.T <- importRawData() }
else load('norm_data.Robj')

# QC Plots
selDataMatrix <- dat.N.T
probeList <- rownames(selDataMatrix)
if (require(lumiHumanAll.db) & require(annotate)) {
	selDataMatrix <- selDataMatrix[!is.na(getSYMBOL(rownames(selDataMatrix), 'lumiHumanAll.db')), ]
	geneSymbol <- getSYMBOL(probeList, 'lumiHumanAll.db')
}

sample.results <- predictSex(as.matrix(selDataMatrix))
geneNames <- getSYMBOL(rownames(selDataMatrix), 'lumiHumanAll.db')
exp.summarized <- summarize(selDataMatrix, probes=as.character(geneNames))

# Reorder predicted sex to match column order in expression matrix
sample.sex <-c()
for (id in colnames(exp.summarized)) {
	sample.sex<-c(sample.sex, sample.results[sample.results$ID==id, "sex"])	
}

sample.factors <- as.numeric(sample.sex=='male')
fit <- fitData(exp.summarized, sample.factors)

plotVolcano <- function(fit) {
	pdf('volcano.pdf', width=10, height=10)
	volcanoplot(fit, names=rownames(fit), cex=0.5, highlight=75)
	dev.off()
}
plotVolcano(fit)
write.table(topTable(fit, number=Inf, p.value=0.05), file='sig_eb_fit.txt', quote=F)
