library(limma)
library(stringr)
library(gplots)
library(gtools)

setwd('/group/stranger-lab/moliva/ImmVar/Robjects')
meta=read.table('../data/pheno/metadata.csv',header=T,sep=',',comment.char = '~')

#load('CD4.joint.norm.fit.Robj')
load('phen.Robj')

phen.cd4 <- phen[phen$CellType=="CD4TNve",]
phen.cd4 <- phen.cd4[mixedorder(as.character(phen.cd4$ImmVarID2)),]
phen.cd14 <- phen[phen$CellType=="CD14+16-Mono",]
phen.cd14 <- phen.cd14[mixedorder(as.character(phen.cd14$ImmVarID2)),]

#x<-prcomp(t(exp_genes.cd4.joint.norm))
#pc1 <- x$x[,1]

#pc1
load('CD14.joint.norm.fit.Robj')
load('CD4.joint.norm.fit.Robj')
load('CD4.joint.norm.exp_genes.Robj')
load('CD14.joint.norm.exp_genes.Robj')
cd4.names <- apply(as.matrix(colnames(exp_genes.cd4.joint.norm), ncol=1), 1, function(x) {unlist(strsplit(x, split=':'))[3]})
cd14.names <- apply(as.matrix(colnames(exp_genes.cd14.joint.norm), ncol=1), 1, function(x) {unlist(strsplit(x, split=':'))[3]})

load('/group/stranger-lab/immvar_data/sv.CD4.Robj')
cd4.svs <- modSv[mixedorder(cd4.names),-1:-2]
load('/group/stranger-lab/immvar_data/sv.CD14.Robj')
cd14.svs <- modSv[mixedorder(cd14.names),-1:-2]
#cd14.svs <- cd14.joint.norm.fit$design[,-1:-2]

cd4.immvar.ids=as.character(phen.cd4$ImmVarID2)
cd14.immvar.ids=as.character(phen.cd14$ImmVarID2)

meta.cd4 <- meta[meta$Study.ID%in%cd4.immvar.ids,]
meta.cd14 <- meta[meta$Study.ID%in%cd14.immvar.ids,]

cd4.exp <- exp_genes.cd4.joint.norm[, mixedorder(cd4.names)]
cd14.exp <- exp_genes.cd14.joint.norm[, mixedorder(cd14.names)]

#cd4.pc <- prcomp(t(cd4.exp))
#cd4.pc1 <- cd4.pc$x[,1]

#cd14.pc <- prcomp(t(cd14.exp))
#cd14.pc1 <- cd14.pc$x[,1]

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
#bp[131] <- "NA/NA"
bp[244] <- "NA/NA"
new.bp <- c()
for (i in bp) {
	new.bp<-c(new.bp, unlist(strsplit(i, split="[/]")))
}

meta.cd4.tmp$BP.sys <- new.bp[seq(1, length(new.bp), by=2)]
meta.cd4.tmp$BP.dia <- new.bp[seq(2, length(new.bp), by=2)]

covs <- c(7,10,11,12,13,18,19)
cov.mat <- matrix(nrow=ncol(cd4.svs), ncol=length(covs))
pc.cov.mat <- matrix(nrow=ncol(cd4.svs), ncol=length(covs))
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
#pdf(file='/scratch/t.cczysz/cd4_heatmap.pdf')
pdf(file='/group/stranger-lab/czysz/ImmVar/plots/cd4_svavscovs.pdf')
#heatmap(cov.mat[,-4], col=cm.colors(256))
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
pval.mat <- matrix(nrow=ncol(cd14.svs), ncol=length(covs))
for (cov in seq(length(covs))) {
	for (sv in seq(ncol(cd14.svs))) {
		x <- lm(cd14.svs[, sv] ~ meta.cd14.tmp[, covs[cov]])
		cov.mat[sv, cov] <- summary.lm(x)$r.squared
		pval.mat[sv, cov] <- summary(x)$coefficients[2,4]
		
}
}
colnames(cov.mat) <- names(meta.cd14.tmp)[covs]

pop <- phen.cd14$Race
race.covs <- c()
pval.race <- c()
for (sv in seq(ncol(cd14.svs))) {
	x <- lm(cd14.svs[,sv] ~ pop)
	race.covs <- c(race.covs, summary.lm(x)$r.squared)
	pval.race <- c(pval.race, summary(x)$coefficients[2,4])
}
cov.mat <- cbind(cov.mat, pop=race.covs)
pval.mat <- cbind(pval.mat, pop=pval.race)

batch <- phen.cd14$Batch
batch.covs <- c()
batch.pval <- c()
for (sv in seq(ncol(cd14.svs))) {
	x <- lm(cd14.svs[,sv] ~ batch)
	batch.covs <- c(batch.covs, summary.lm(x)$r.squared)
	batch.pval <- c(batch.pval, summary(x)$coefficients[2,4])
}
cov.mat <- cbind(cov.mat, batch=batch.covs)
pval.mat <- cbind(pval.mat, batch=batch.pval)

frozen <- as.numeric(phen.cd14$Frozen)
frozen.covs <- c()
frozen.pval <- c()
for (sv in seq(ncol(cd14.svs))) {
	x <- lm(cd14.svs[,sv] ~ frozen)
	frozen.covs <- c(frozen.covs, summary.lm(x)$r.squared)
	frozen.pval <- c(frozen.pval, summary(x)$coefficients[2,4])
}
cov.mat <- cbind(cov.mat, frozen=frozen.covs)
pval.mat <- cbind(pval.mat, frozen=frozen.pval)
pdf(file='/group/stranger-lab/czysz/ImmVar/plots/cd14_svavscovs.pdf', width=10, height=10)
#heatmap(cov.mat[,-4], col=cm.colors(256))
heatmap.2(cov.mat[,-4],cellnote=signif(pval.mat[,-4],2), notecex=1, notecol='black',col=heat.colors(20)[seq(20,1)],trace='none', Rowv=NULL, key.xlab='R^2', margins=c(7,7))
dev.off()

#x<-prcomp(t(exp_genes.cd4.joint.norm))
