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

load('CD14.joint.norm.fit.Robj')
load('CD4.joint.norm.fit.Robj')
load('CD4.joint.norm.exp_genes.Robj')
load('CD14.joint.norm.exp_genes.Robj')
cd4.names <- apply(as.matrix(colnames(exp_genes.cd4.joint.norm), ncol=1), 1, function(x) {unlist(strsplit(x, split=':'))[3]})
cd14.names <- apply(as.matrix(colnames(exp_genes.cd14.joint.norm), ncol=1), 1, function(x) {unlist(strsplit(x, split=':'))[3]})

cd4.pcs<-prcomp(t(exp_genes.cd4.joint.norm))
cd4.pcs <- cd4.pcs$x[,1:20]

cd14.pcs<-prcomp(t(exp_genes.cd14.joint.norm))
cd14.pcs <- cd14.pcs$x[,1:20]

load('/group/stranger-lab/immvar_data/sv.CD4.Robj')
cd4.svs <- modSv[mixedorder(cd4.names),-1:-2]
load('/group/stranger-lab/immvar_data/sv.CD14.Robj')
cd14.svs <- modSv[mixedorder(cd14.names),-1:-2]

cd4.immvar.ids=as.character(phen.cd4$ImmVarID2)
cd14.immvar.ids=as.character(phen.cd14$ImmVarID2)

meta.cd4 <- meta[meta$Study.ID%in%cd4.immvar.ids,]
meta.cd14 <- meta[meta$Study.ID%in%cd14.immvar.ids,]

cd4.exp <- exp_genes.cd4.joint.norm[, mixedorder(cd4.names)]
cd14.exp <- exp_genes.cd14.joint.norm[, mixedorder(cd14.names)]

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
covs <- c(7,10,11,12,13,18,19)

bp <- as.character(meta.cd4.tmp$Blood.Pressure)
bp[244] <- "NA/NA"
new.bp <- c()
for (i in bp) {
	new.bp<-c(new.bp, unlist(strsplit(i, split="[/]")))
}

meta.cd4.tmp$BP.sys <- new.bp[seq(1, length(new.bp), by=2)]
meta.cd4.tmp$BP.dia <- new.bp[seq(2, length(new.bp), by=2)]

buildCor <- function(cov, mat) {
	cov.vec <- c()
	p.vec <- c()
	for (vec in seq(ncol(mat))){
		cov.vec <- c(cov.vec,summary.lm(lm(mat[,vec] ~ cov))$r.squared)
}
	return(cov.vec)
}
covs <- c(7,10,11,12,13,18,19)
cov.mat <- matrix(nrow=ncol(cd4.svs), ncol=length(covs))
colnames(cov.mat) <- names(meta.cd4.tmp)[covs]

# Build matrix of individual x covariate
#meta.cd4.tmp <- cbind(meta.cd4.tmp[,covs], pop=phen.cd4$Race, batch=phen.cd4$Batch, frozen=as.numeric(phen.cd4$Frozen))
pc.cov.mat <- matrix(nrow=20, ncol=length(covs))
colnames(pc.cov.mat) <- names(meta.cd4.tmp)[covs]

if (T) {
for (cov in seq(length(covs))) {
	cov.vec <- buildCor(meta.cd4.tmp[, covs[cov]], cd4.svs)
	cov.mat[, cov] <- cov.vec	
	cov.vec <- buildCor(meta.cd4.tmp[, covs[cov]], cd4.pcs)
	pc.cov.mat[, cov] <- cov.vec
	if (F) {
	for (sv in seq(ncol(cd4.svs))) {
		x <- lm(cd4.svs[, sv] ~ meta.cd4.tmp[, covs[cov]])
		cov.mat[sv, cov] <- summary.lm(x)$r.squared
		
		}
	}
}
}


pop <- phen.cd4$Race
race.covs <- c()
for (sv in seq(ncol(cd4.svs))) {
	x <- lm(cd4.svs[,sv] ~ pop)
	race.covs <- c(race.covs, summary.lm(x)$r.squared)
}
race.covs <- buildCor(pop, cd4.svs)
cov.mat <- cbind(cov.mat, pop=race.covs)
race.covs <- buildCor(pop, cd4.pcs)
pc.cov.mat <- cbind(pc.cov.mat, pop=race.covs)

batch <- phen.cd4$Batch
batch.covs <- c()
for (sv in seq(ncol(cd4.svs))) {
	x <- lm(cd4.svs[,sv] ~ batch)
	batch.covs <- c(batch.covs, summary.lm(x)$r.squared)
}
batch.covs <- buildCor(batch, cd4.svs)
cov.mat <- cbind(cov.mat, batch=batch.covs)
batch.covs <- buildCor(batch, cd4.pcs)
cov.mat <- cbind(pc.cov.mat, batch=batch.covs)

frozen <- as.numeric(phen.cd4$Frozen)
frozen.covs <- c()
for (sv in seq(ncol(cd4.svs))) {
	x <- lm(cd4.svs[,sv] ~ frozen)
	frozen.covs <- c(frozen.covs, summary.lm(x)$r.squared)
}
frozen.covs <- buildCor(frozen, cd4.svs)
cov.mat <- cbind(cov.mat, frozen=frozen.covs)
frozen.covs <- buildCor(frozen, cd4.pcs)
pc.cov.mat <- cbind(pc.cov.mat, frozen=frozen.covs)
#pdf(file='/scratch/t.cczysz/cd4_heatmap.pdf')
pdf(file='/group/stranger-lab/czysz/ImmVar/plots/cd4_pcavscovs.pdf')
#heatmap(cov.mat[,-4], col=cm.colors(256))
heatmap.2(pc.cov.mat[,-4],trace='none')
dev.off()
q()

cov.mat <- matrix(nrow=ncol(cd4.svs), ncol=length(covs))
pval.mat <- matrix(nrow=ncol(cd4.svs), ncol=length(covs))
for (cov in seq(length(covs))) {
	for (sv in seq(ncol(cd4.svs))) {
		x <- lm(cd4.svs[, sv] ~ meta.cd4.tmp[, covs[cov]])
		cov.mat[sv, cov] <- summary.lm(x)$r.squared
		pval.mat[sv, cov] <- summary(x)$coefficients[2,4]
		
}
}
colnames(cov.mat) <- names(meta.cd4.tmp)[covs]

pop <- phen.cd4$Race
race.covs <- c()
pval.race <- c()
for (sv in seq(ncol(cd4.svs))) {
	x <- lm(cd4.svs[,sv] ~ pop)
	race.covs <- c(race.covs, summary.lm(x)$r.squared)
	pval.race <- c(pval.race, summary(x)$coefficients[2,4])
}
cov.mat <- cbind(cov.mat, pop=race.covs)
pval.mat <- cbind(pval.mat, pop=pval.race)
race.cov <- buildCor(pop, cd4.pcs)
pc.cov.mat <- cbind(pc.cov.mat, pop=race.covs)

batch <- phen.cd4$Batch
batch.covs <- c()
batch.pval <- c()
for (sv in seq(ncol(cd4.svs))) {
	x <- lm(cd4.svs[,sv] ~ batch)
	batch.covs <- c(batch.covs, summary.lm(x)$r.squared)
	batch.pval <- c(batch.pval, summary(x)$coefficients[2,4])
}
cov.mat <- cbind(cov.mat, batch=batch.covs)
pval.mat <- cbind(pval.mat, batch=batch.pval)

frozen <- as.numeric(phen.cd4$Frozen)
frozen.covs <- c()
frozen.pval <- c()
for (sv in seq(ncol(cd4.svs))) {
	x <- lm(cd4.svs[,sv] ~ frozen)
	frozen.covs <- c(frozen.covs, summary.lm(x)$r.squared)
	frozen.pval <- c(frozen.pval, summary(x)$coefficients[2,4])
}
cov.mat <- cbind(cov.mat, frozen=frozen.covs)
pval.mat <- cbind(pval.mat, frozen=frozen.pval)
pdf(file='/group/stranger-lab/czysz/ImmVar/plots/cd4_svavscovs.pdf', width=10, height=10)
#heatmap(cov.mat[,-4], col=cm.colors(256))
heatmap.2(cov.mat[,-4],cellnote=signif(pval.mat[,-4],2), notecex=1, notecol='black',col=heat.colors(20)[seq(20,1)],trace='none', Rowv=NULL, key.xlab='R^2', margins=c(7,7))
dev.off()

# CD14
#meta.cd14.tmp <- cbind(meta.cd14.tmp[,covs], pop=phen.cd14$Race, batch=phen.cd14$Batch, frozen=as.numeric(phen.cd14$Frozen))
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
pc.cov.mat <- matrix(nrow=20, ncol=length(covs))
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
race.cov <- buildCor(pop, cd14.pcs)
pc.cov.mat <- cbind(pc.cov.mat, pop=race.covs)

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
