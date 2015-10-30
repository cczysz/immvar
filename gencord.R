library(BiocInstaller)
library(lumi)
library(limma)
#biocLite('lumiHumanIDMapping')
library(lumiHumanIDMapping)
#biocLite('lumiHumanAll.db')
library(lumiHumanAll.db)
#library(peer)
library(sva)


RunPeer <- function(expression, k=20, covs) {
	if (F) {
        # Expression must be in NxG. N number of samples, G number of genes
        model = PEER()
        PEER_setPhenoMean(model,t(expression))
        PEER_setNk(model,k)
	if (!is.null(covs)) {PEER_setCovariates(model,as.matrix(covs))}
        PEER_update(model)
        peer.factors = PEER_getX(model)
        return(peer.factors)
	}
	mod = model.matrix(~0+as.factor(covs))
	mod0 = model.matrix(~1)
	n.sv <- num.sv(expression, mod, method="leek")
	svobj <- sva(expression, mod, mod0, n.sv=n.sv)

	modSv <- cbind(mod, svobj$sv)
	fit <- lmFit(expression, modSv)
	fit <- eBayes(fit)
	return(fit)
	#return(svobj)
}

MakeResiduals <- function(input.row,peer.factors) {
        fit <- lm(input.row ~ peer.factors[, -1] )
        residuals(fit)
}

PerformDE <- function(eset, covs) {
	design <- model.matrix(~ 0 + as.factor(covs))
	colnames(design) <- c("Female", "Male")
	fit <- lmFit(eset, design)
	fit$genes <- data.frame(ID= probeList, geneSymbol=geneSymbol)
	efit <- eBayes(fit, robust=T)
	return(list(fit=fit, efit=efit))
}

#dat <- readBeadSummaryData('GSE17080_unnormalized_T-cells.txt', skip=0, columns=list(exprs="AVG_Signal", nObservations = "Avg_NBEADS"))
#dat <- lumiR('GSE17080_unnormalized_T-cells.txt', lib.mapping = 'lumiHumanIDMapping', columnNameGrepPattern=list(detection=NA, beadNum=NA))
#dat <- lumiR('GSE17080_unnormalized_T-cells.txt', lib.mapping = 'lumiHumanIDMapping', columnNameGrepPattern=list(exprs="AVG_Signal", beadNum="Avg_NBEADS", detection=NA))
#eset <- log2(exprs(dat))
#dat.N.T <- lumiExpresso(dat, bg.correct = F, varianceStabilize.param = list(method='log2'), normalize.param = list(method='quantile'))
#rm(dat)
#save(dat.N.T, file='norm_data.Robj')
load('norm_data.Robj')
#norm.eset <- normaliseIllumina(exprs(dat), method="quantile", transform='log2')

#eset <- exprs(dat)
#log2.eset <- lumiT(eset, method="log2")
#norm.eset <- lumiN(log2.eset)

sample_data <- read.table(file='E-GEOD-17080.sdrf.txt', header=T, sep="\t")
sample_data <- sample_data[sample_data$Factor.Value..CELL.TYPE. == 'Primary T-cells', ]
conversions <- read.table(file='ids.txt', header=F)

selDataMatrix <- exprs(dat.N.T)

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

#sampleFactors <- sample_data[, "Factor.Value..GENDER."]
sampleFactors <- as.numeric(sampleFactors=="Male")

if (F) {
design <- model.matrix(~0+factor(sampleFactors))
colnames(design) <- c("Female", "Male")
fit <- lmFit(selDataMatrix, design)
fit <- eBayes(fit)
fit$genes <- data.frame(ID= probeList, geneSymbol=geneSymbol)
 }

mod = model.matrix(~0+as.factor(sampleFactors))
colnames(mod) <- c('Female', 'Male')
mod0 = model.matrix(~1, data=pData(dat.N.T))
#n.sv <- num.sv(selDataMatrix, mod, method="leek")
svobj <- sva(selDataMatrix, mod, mod0)

modSv <- cbind(mod, svobj$sv)
fit <- lmFit(selDataMatrix, modSv)
contrast.matrix <- c(-1,1,rep(0, svobj$n.sv))
#contrast.matrix <- makeContrasts(Male-Female, levels=modSv[,1:2])
contrast.fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(contrast.fit)
#fit <- RunPeer(selDataMatrix, 20, covs=sampleFactors)

top_diff <- topTable(fit, number=Inf, p.value=0.05)
top_diff$gene <- geneSymbol[rownames(top_diff)]
#peer.factors <- RunPeer(selDataMatrix, 20, covs=sampleFactors)

#eset.res <- apply(as.matrix(tcell.eset), 1, MakeResiduals, peer.factors=peer.factors)
#eset.res <- t(eset.res)
#de.out <- PerformDE(eset.res, covs=as.numeric(bcell.sex=='Male'))
