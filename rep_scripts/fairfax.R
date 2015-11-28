library(lumi)
library(limma)
library(massiR)
library(oligo)
library(lumiHumanIDMapping)
library(illuminaHumanv4.db)
library(lumiHumanAll.db)
library(sva)
library(genefilter)
library(annotate)

predictSex <- function(temp.dat, plot_massi=F) {
        ## Verify massiR results
        y_probes <- data.frame(row.names=c("Ng17goGe3f5AsKKHCo", "oE68.e0uHwT2v8UlJI", "KuyUIoSDxODhLlSLhI", "E35Tklf.OiQtEMSSfk",
                "KLpuVfpdEpAgwfAJ4Y", "TN6X5VSSxLX1wwao50", "NnMe15_xBOv6rgqigI","T7rNJ6HSwRV.fSeqOo",
                "6VE8ScX3J.LOsufT0E", "reXjDTnq9wwlcfUXDc", "WN_KCVC5sX8._VUUV0", "Z_OQlcBHf2HUiiCekk", "QEpD_7OuProlBdrbvo"))
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

setwd('/group/stranger-lab/immvar_rep/emtab2232/')
if (!file.exists('cd14_norm_exp.Robj')) {
	dat <- lumiR(fileName='mono_raw_286_13_10_2010.txt',
		lib.mapping = 'lumiHumanIDMapping',
		columnNameGrepPattern=list(exprs='AVG_Signal', se.exprs='BEAD_STDEV', detection='Detection', beadNum='Avg_NBEADS'))
	#dat.N <- lumiN(dat, method="quantile")
	dat.N <- lumiExpresso(dat, varianceStabilize.param=list(method='vst'), normalize.param=list(method="quantile"))
	save(dat.N, file="cd14_norm_exp.Robj")
} else load(file='cd14_norm_exp.Robj')

eset <- exprs(dat.N)
fdat <- fData(dat.N)
sample.results <- predictSex(eset)

load('/group/stranger-lab/immvar_rep/mappings/ilmnv4.Robj')
load('/group/stranger-lab/immvar_rep/mappings/annot.Robj')

if (!file.exists('fairfax_monocytes_summarized.Robj')) {
	x <- c()
	for (id in filt.probes$ProbeID) {
		nuid <- as.character(rownames(fdat)[fdat$PROBE_ID==id])
		if (length(nuid) == 0) { x <- c(x, NA)
		} else {x <- c(x, nuid) }
	}
	filt.probes$NuID <- x
	rm(x)

	tempMat <- NULL
	ensIDs <- c()
	for (gene in unique(filt.probes$EnsemblID)) {
		info <- filt.probes[filt.probes$EnsemblID==gene,]
		tempMat <- rbind(tempMat, eset[rownames(eset)%in%info$NuID,])
		ensIDs <- c(ensIDs, rep(gene, sum(rownames(eset)%in%info$NuID)))
	}
	exp.summarized <- summarize(tempMat, probes=ensIDs)
	rm(tempMat)
	save(exp.summarized, file="fairfax_monocytes_summarized.Robj")
} else { load(file="fairfax_monocytes_summarized.Robj") }

q('no')
if (F) {
	massi_y_plot(massi.y.out)
	massi_cluster_plot(massi.select.out, results)
	dev.off();dev.off()
}
#males <- as.numeric(sample.results$sex=='male')

# Reorder predicted sex to match expression colnames
pred.sex <- c()
for (id in colnames(exp.summarized)) {
	pred.sex <- c(pred.sex, sample.results[sample.results$ID==paste("X",id,sep=''), "sex"])
}
males <- as.numeric(pred.sex=='male')

# Filter lowly expressed genes
f.m <- pOverA(0.10, quantile(as.numeric(exp.summarized[, !!males]), probs=c(0.1)))
f.f <- pOverA(0.10, quantile(as.numeric(exp.summarized[, !males]), probs=c(0.1)))
ffun <- filterfun(f.m, f.f)
filtered.genes <- genefilter(exp.summarized, ffun)
q()
mod = model.matrix(~0+as.factor(males))
colnames(mod) <- c('Female', 'Male')
mod0 = model.matrix(~0+rep(1,length(males)))

if (!file.exists('fairfax_sv.Robj')) {
	svobj <- sva(na.omit(exp.summarized), mod, mod0)
	save(svobj, file='fairfax_sv.Robj')
} else load(file='fairfax_sv.Robj')
modSv <- cbind(mod, svobj$sv)

if (!file.exists('fairfax.fit')) {
	fit <- lmFit(na.omit(exp.summarized), modSv)
	contrast.matrix <- c(-1,1,rep(0, svobj$n.sv))
	contrast.fit <- contrasts.fit(fit, contrast.matrix)
	fit <- eBayes(contrast.fit)
	fit$genes <- annots[rownames(fit), "symbol_id"]
	fit$chr <- annots[rownames(fit), "chr"]
	save(fit, file='fairfax.fit')
} else load('fairfax.fit')

pdf(file='volcano.pdf')
volcanoplot(fit, names=fit$genes, highlight=50, cex=0.5)
dev.off()
save.image()
