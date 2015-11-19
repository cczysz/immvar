library(limma)
library(stringr)
library(gplots)

setwd('/group/stranger-lab/moliva/ImmVar/Robjects')
meta=read.table('../data/pheno/metadata.csv',header=T,sep=',',comment.char = '~')

load('CD4.joint.norm.fit.Robj')
load('phen.Robj')

phen.cd4 <- phen[phen$CellType=="CD4TNve",]
phen.cd14 <- phen[phen$CellType=="CD14+16-Mono",]
cd4.svs <- cd4.joint.norm.fit$design[,-1:-2]


load('CD14.joint.norm.fit.Robj')
cd14.svs <- cd14.joint.norm.fit$design[,-1:-2]

#cd4.immvar.ids=str_extract(rownames(cd4.joint.norm.fit$design),'IGTB[0-9]+');
cd4.immvar.ids=as.character(phen.cd4$ImmVarID2)
cd14.immvar.ids=as.character(phen.cd14$ImmVarID2)
#cd14.immvar.ids=str_extract(rownames(cd14.joint.norm.fit$design),'IGTB[0-9]+');

meta.cd4 <- meta[meta$Study.ID%in%cd4.immvar.ids,]
meta.cd14 <- meta[meta$Study.ID%in%cd14.immvar.ids,]

if (T) {
meta.cd4.tmp <- NULL
for (id in cd4.immvar.ids) {
	meta.cd4.tmp <- rbind(meta.cd4.tmp, meta.cd4[meta.cd4$Study.ID==id,][1,])
}
}

if (T) {
meta.cd14.tmp <- NULL
for (id in cd14.immvar.ids) {
	meta.cd14.tmp <- rbind(meta.cd14.tmp, meta.cd14[meta.cd14$Study.ID==id,][1,])
}
}
#covs <- c(7,8,9,10,11,12,13)
covs <- c(7,10,11,12,13,18,19)

bp <- as.character(meta.cd4.tmp$Blood.Pressure)
#bp[325] <- "NA/NA"
bp[131] <- "NA/NA"
new.bp <- c()
for (i in bp) {
	new.bp<-c(new.bp, unlist(strsplit(i, split="[/]")))
}

meta.cd4.tmp$BP.sys <- new.bp[seq(1, length(new.bp), by=2)]
meta.cd4.tmp$BP.dia <- new.bp[seq(2, length(new.bp), by=2)]

covs <- c(7,10,11,12,13,18,19)
cov.mat <- matrix(nrow=ncol(cd4.svs), ncol=length(covs))
for (cov in seq(length(covs))) {
	for (sv in seq(ncol(cd4.svs))) {
		x <- lm(cd4.svs[, sv] ~ meta.cd4.tmp[, covs[cov]])
		cov.mat[sv, cov] <- summary.lm(x)$r.squared
		
}
}
colnames(cov.mat) <- names(meta.cd4.tmp)[covs]

pop <- phen.cd4$Race
race.covs <- c()
for (sv in seq(ncol(cd4.svs))) {
	x <- lm(cd4.svs[,sv] ~ pop)
	race.covs <- c(race.covs, summary.lm(x)$r.squared)
}
cov.mat <- cbind(cov.mat, pop=race.covs)

batch <- phen.cd4$Batch
batch.covs <- c()
for (sv in seq(ncol(cd4.svs))) {
	x <- lm(cd4.svs[,sv] ~ batch)
	batch.covs <- c(batch.covs, summary.lm(x)$r.squared)
}
cov.mat <- cbind(cov.mat, batch=batch.covs)

frozen <- as.numeric(phen.cd4$Batch)
frozen.covs <- c()
for (sv in seq(ncol(cd4.svs))) {
	x <- lm(cd4.svs[,sv] ~ frozen)
	frozen.covs <- c(frozen.covs, summary.lm(x)$r.squared)
}
cov.mat <- cbind(cov.mat, frozen=frozen.covs)
pdf(file='/scratch/t.cczysz/cd4_heatmap.pdf')
heatmap(cov.mat[,-4], col=cm.colors(256))
heatmap.2(cov.mat[,-4],trace='none')
dev.off()

bp <- as.character(meta.cd14.tmp$Blood.Pressure)
#bp[325] <- "NA/NA"
#bp[131] <- "NA/NA"
new.bp <- c()
for (i in bp) {
	new.bp<-c(new.bp, unlist(strsplit(i, split="[/]")))
}

meta.cd14.tmp$BP.sys <- new.bp[seq(1, length(new.bp), by=2)]
meta.cd14.tmp$BP.dia <- new.bp[seq(2, length(new.bp), by=2)]

covs <- c(7,10,11,12,13,18,19)
cov.mat <- matrix(nrow=ncol(cd14.svs), ncol=length(covs))
for (cov in seq(length(covs))) {
	for (sv in seq(ncol(cd14.svs))) {
		x <- lm(cd14.svs[, sv] ~ meta.cd14.tmp[, covs[cov]])
		cov.mat[sv, cov] <- summary.lm(x)$r.squared
		
}
}
colnames(cov.mat) <- names(meta.cd14.tmp)[covs]

pop <- phen.cd14$Race
race.covs <- c()
for (sv in seq(ncol(cd14.svs))) {
	x <- lm(cd14.svs[,sv] ~ pop)
	race.covs <- c(race.covs, summary.lm(x)$r.squared)
}
cov.mat <- cbind(cov.mat, pop=race.covs)

batch <- phen.cd14$Batch
batch.covs <- c()
for (sv in seq(ncol(cd14.svs))) {
	x <- lm(cd14.svs[,sv] ~ batch)
	batch.covs <- c(batch.covs, summary.lm(x)$r.squared)
}
cov.mat <- cbind(cov.mat, batch=batch.covs)

frozen <- as.numeric(phen.cd14$Batch)
frozen.covs <- c()
for (sv in seq(ncol(cd14.svs))) {
	x <- lm(cd14.svs[,sv] ~ frozen)
	frozen.covs <- c(frozen.covs, summary.lm(x)$r.squared)
}
cov.mat <- cbind(cov.mat, frozen=frozen.covs)
pdf(file='/scratch/t.cczysz/cd14_heatmap.pdf')
heatmap(cov.mat[,-4], col=cm.colors(256))
heatmap.2(cov.mat[,-4],trace='none')
dev.off()
