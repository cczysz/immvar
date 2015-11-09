library(lumi)
library(limma)
library(massiR)
library(lumiHumanIDMapping)
library(illuminaHumanv4.db)
library(lumiHumanAll.db)
library(sva)
library(genefilter)

setwd('/group/stranger-lab/immvar_rep/emtab2232/')
if (!file.exists('norm_cd14.Robj')) {
	proc_dat <- read.table(file='CD14.47231.414.b.txt', header=T, row.names=1)
	save(proc_dat, file='norm_cd14.Robj')
} else load(file='norm_cd14.Robj')

x <- illuminaHumanv4CHR
probes_chr <- mappedkeys(x)
xx <- as.list(x[probes_chr])

y_probes <- data.frame(row.names=names(xx)[xx=="Y"])
massi.y.out <- massi_y(proc_dat, y_probes)
massi.select.out <- massi_select(proc_dat, y_probes, threshold=4)
results <- massi_cluster(massi.select.out)
sample.results <- data.frame(results[[2]])

if (F) {
	massi_y_plot(massi.y.out)
	massi_cluster_plot(massi.select.out, results)
	dev.off();dev.off()
}
males <- as.numeric(sample.results$sex=='male')

# Reorder predicted sex to match expression colnames
pred.sex <- c()
for (id in colnames(proc_dat)) {
	pred.sex <- c(pred.sex, sample.results[sample.results$ID==id, "sex"])
}

#x <- illuminaHumanv4ENSEMBL
#mapped_probes <- mappedkeys(x)
#xx <- as.list(x[mapped_probes])
probeList <- rownames(proc_dat)
if (require(lumiHumanAll.db) & require(annotate)) {
        geneSymbol <- getSYMBOL(probeList, 'illuminaHumanv4.db')
        selDataMatrix <- proc_dat[!is.na(geneSymbol),]
}

selDataMatrix <- as.matrix(selDataMatrix)
# Filter genes by expression
f1 <- pOverA(p=0.3, A=quantile(as.numeric(selDataMatrix), probs=c(0.1)))
ffun <- filterfun(f1)
filtered <- genefilter(selDataMatrix, ffun)
selDataMatrix <- selDataMatrix[filtered, ]

males <- as.numeric(pred.sex=='male')
mod = model.matrix(~0+as.factor(males))
colnames(mod) <- c('Female', 'Male')
mod0 = model.matrix(~0+rep(1,length(males)))
if (file.exists('fairfax_sv.Robj')) {
	#n.sv <- num.sv(selDataMatrix, mod, method="leek")
	svobj <- sva(as.matrix(selDataMatrix), mod, mod0, n.sv=33)
	save(svobj, file='fairfax_sv.Robj')
} else load(file='fairfax_sv.Robj')
modSv <- cbind(mod, svobj$sv)

if (!file.exists('fairfax.fit')) {
	fit <- lmFit(proc_dat, modSv)
	contrast.matrix <- c(-1,1,rep(0, svobj$n.sv))
	contrast.fit <- contrasts.fit(fit, contrast.matrix)
	fit <- eBayes(contrast.fit)
	fit$genes <- geneSymbol
	save(fit, file='fairfax.fit')
} else load('fairfax.fit')

pdf(file='volcano.pdf')
volcanoplot(fit, names=fit$genes, highlight=50, cex=0.5)
dev.off()
