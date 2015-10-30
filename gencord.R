library(BiocInstaller)
library(lumi)
library(limma)
#biocLite('lumiHumanIDMapping')
library(lumiHumanIDMapping)
#biocLite('lumiHumanAll.db')
library(lumiHumanAll.db)
#library(peer)
library(sva)

# Set to true if Robj needs to be replaced
if (T) {
dat <- lumiR('GSE17080_unnormalized_T-cells.txt', 
	lib.mapping = 'lumiHumanIDMapping',
	columnNameGrepPattern=list(detection=NA, beadNum=NA))
dat.N.T <- lumiExpresso(dat, bg.correct = F, varianceStabilize.param = list(method='log2'), normalize.param = list(method='quantile'))
save(dat.N.T, file='norm_data.Robj')
}
load('norm_data.Robj')

# QC Plots
if (T) {
	pdf('raw_density.pdf')
	plot(dat, what='density')
	dev.off()

	#pdf('cov_density.pdf')
	#plot(dat, what='cv')
	#dev.off()

	pdf('sample_relation.pdf', width=12, height=8)
	plot(dat, what='sampleRelation')
	plot(dat.N.T, what='sampleRelation')
	plot(dat, what='sampleRelation', method='mds')
	dev.off()

}
sample_data <- read.table(file='E-GEOD-17080.sdrf.txt', header=T, sep="\t")
sample_data <- sample_data[sample_data$Factor.Value..CELL.TYPE. == 'Primary T-cells', ]
conversions <- read.table(file='ids.txt', header=F)

selDataMatrix <- exprs(dat.N.T)

# Convert between IDs in pData and colNames in raw data file
ids<-c()
for (id in colnames(selDataMatrix)) {
	ids<-c(ids, strsplit(id, split="_")[[1]][1])
}

sample_data[, 1] <- unlist(strsplit(as.character(sample_data[,1]), split=" "))[seq(from=1, to=170, by=2)]

sample_data_ids<-c()
for (id in sample_data[,1]) {
	sample_data_ids<-c(sample_data_ids, as.character(conversions[conversions[,1]==id,2]))
}
sample_data<-cbind(sample_data, ids=sample_data_ids)

probeList <- rownames(selDataMatrix)
if (require(lumiHumanAll.db) & require(annotate)) {
	selDataMatrix <- selDataMatrix[!is.na(getSYMBOL(rownames(selDataMatrix), 'lumiHumanAll.db')), ]
	geneSymbol <- getSYMBOL(probeList, 'lumiHumanAll.db')
}

sampleFactors <- c()
for (id in ids) {
	sampleFactors<-rbind(sampleFactors, as.character(sample_data[sample_data$ids==id, "Characteristics..gender."]))
}

sampleFactors <- as.numeric(sampleFactors=="Male")

mod = model.matrix(~0+as.factor(sampleFactors))
colnames(mod) <- c('Female', 'Male')

mod0 = model.matrix(~1, data=pData(dat.N.T))

#n.sv <- num.sv(selDataMatrix, mod, method="leek")
svobj <- sva(selDataMatrix, mod, mod0)

modSv <- cbind(mod, svobj$sv)

fit <- lmFit(selDataMatrix, modSv)
contrast.matrix <- c(-1,1,rep(0, svobj$n.sv))
contrast.fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(contrast.fit)

fit$probeID <- fData(dat.N.T)[rownames(fit),]
fit$genes <- geneSymbol[rownames(fit)]
top_diff <- topTable(fit, number=Inf, p.value=0.05)
top_diff$gene <- geneSymbol[rownames(top_diff)]

pdf('volcano.pdf')
volcanoplot(fit, names=fit$genes, cex=0.5, highlight=50)
dev.off()
write.table(topTable(fit, number=Inf), file='eb_fit.txt', quote=F)
