library(limma)
library(stringr)
library(gplots)
library(gtools)

setwd('/group/stranger-lab/moliva/ImmVar/Robjects')
meta=read.table('../data/pheno/metadata.csv',header=T,sep=',',comment.char = '~')

load('phen.Robj')

# Load phenotypes and order by ImmVar ID
phen.cd4 <- phen[phen$CellType=="CD4TNve",]
phen.cd4 <- phen.cd4[mixedorder(as.character(phen.cd4$ImmVarID2)),]
phen.cd14 <- phen[phen$CellType=="CD14+16-Mono",]
phen.cd14 <- phen.cd14[mixedorder(as.character(phen.cd14$ImmVarID2)),]

# Load fit objects
#load('CD14.joint.norm.fit.Robj')
#load('CD4.joint.norm.fit.Robj')

load('CD4.joint.norm.exp_genes.Robj')
load('CD14.joint.norm.exp_genes.Robj')
# Rename column names to ImmVar ID
cd4.names <- apply(as.matrix(colnames(exp_genes.cd4.joint.norm), ncol=1), 1, function(x) {unlist(strsplit(x, split=':'))[3]})
cd14.names <- apply(as.matrix(colnames(exp_genes.cd14.joint.norm), ncol=1), 1, function(x) {unlist(strsplit(x, split=':'))[3]})

# Reorder expression by name
cd4.exp <- exp_genes.cd4.joint.norm[, mixedorder(cd4.names)]
cd14.exp <- exp_genes.cd14.joint.norm[, mixedorder(cd14.names)]

# Calculate PCs of expression data
cd4.pcs<-prcomp(t(exp_genes.cd4.joint.norm))
cd4.pcs <- cd4.pcs$x[,1:20]

cd14.pcs<-prcomp(t(exp_genes.cd14.joint.norm))
cd14.pcs <- cd14.pcs$x[,1:20]

# Save surrogate variables
load('/group/stranger-lab/immvar_data/sv.CD4.Robj')
cd4.svs <- modSv[mixedorder(cd4.names),-1:-2]
load('/group/stranger-lab/immvar_data/sv.CD14.Robj')
cd14.svs <- modSv[mixedorder(cd14.names),-1:-2]

cd4.immvar.ids=as.character(phen.cd4$ImmVarID2)
cd14.immvar.ids=as.character(phen.cd14$ImmVarID2)

# Pull meta data separately for each celltype for IDs present in expression
meta.cd4 <- meta[meta$Study.ID%in%cd4.immvar.ids,]
meta.cd14 <- meta[meta$Study.ID%in%cd14.immvar.ids,]

# Remove duplicate entries
meta.cd4.tmp <- NULL
for (id in cd4.immvar.ids) {
	meta.cd4.tmp <- rbind(meta.cd4.tmp, meta.cd4[meta.cd4$Study.ID==id,][1,])
}

meta.cd14.tmp <- NULL
for (id in cd14.immvar.ids) {
	meta.cd14.tmp <- rbind(meta.cd14.tmp, meta.cd14[meta.cd14$Study.ID==id,][1,])
}

bp <- as.character(meta.cd4.tmp$Blood.Pressure)
bp[244] <- "NA/NA"
new.bp <- c()
for (i in bp) {
	new.bp<-c(new.bp, unlist(strsplit(i, split="[/]")))
}
meta.cd4.tmp$BP.sys <- as.numeric(new.bp[seq(1, length(new.bp), by=2)])
meta.cd4.tmp$BP.dia <- as.numeric(new.bp[seq(2, length(new.bp), by=2)])

bp <- as.character(meta.cd14.tmp$Blood.Pressure)
new.bp <- c()
for (i in bp) {
	new.bp<-c(new.bp, unlist(strsplit(i, split="[/]")))
}

meta.cd14.tmp$BP.sys <- as.numeric(new.bp[seq(1, length(new.bp), by=2)])
meta.cd14.tmp$BP.dia <- as.numeric(new.bp[seq(2, length(new.bp), by=2)])

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
cov.pval.mat <- matrix(nrow=ncol(cd4.svs), ncol=length(covs))
colnames(cov.pval.mat) <- names(meta.cd4.tmp)[covs]

# Build matrix of individual x covariate
pc.cov.mat <- matrix(nrow=20, ncol=length(covs))
colnames(pc.cov.mat) <- names(meta.cd4.tmp)[covs]
pc.pval.mat <- matrix(nrow=20, ncol=length(covs))
colnames(pc.pval.mat) <- names(meta.cd4.tmp)[covs]

if (T) {
for (cov in seq(length(covs))) {
	for (sv in seq(ncol(cd4.svs))) {
		x <- lm(cd4.svs[, sv] ~ meta.cd4.tmp[, covs[cov]])
		cov.mat[sv, cov] <- summary.lm(x)$r.squared
		cov.pval.mat[sv, cov] <- summary.lm(x)$coefficients[2,4]	
		}
	for (pc in seq(ncol(cd4.pcs))) {
		x <- lm(cd4.pcs[, pc] ~ meta.cd4.tmp[, covs[cov]])
		pc.cov.mat[pc, cov] <- summary.lm(x)$r.squared
		pc.pval.mat[pc, cov] <- summary.lm(x)$coefficients[2,4]	
		}
}
}

pop <- as.character(phen.cd4$Race)
pop <- replace(pop, pop=='Caucasian', 0)
pop <- replace(pop, pop=='African-American', 1)
pop <- replace(pop, pop=='Asian', 2)
pop <- as.numeric(pop)
race.covs <- c()
race.pval <- c()
for (sv in seq(ncol(cd4.svs))) {
	x <- lm(cd4.svs[,sv] ~ pop)
	race.covs <- c(race.covs, summary.lm(x)$r.squared)
	race.pval <- c(race.pval, summary.lm(x)$coefficients[2,4])
}
cov.mat <- cbind(cov.mat, pop=race.covs)
cov.pval.mat <- cbind(cov.pval.mat, pop=race.pval)

race.covs <- c()
race.pval <- c()
for (pc in seq(ncol(cd4.pcs))) {
	x <- lm(cd4.pcs[,pc] ~ pop)
	race.covs <- c(race.covs, summary.lm(x)$r.squared)
	race.pval <- c(race.pval, summary.lm(x)$coefficients[2,4])
}
pc.cov.mat <- cbind(pc.cov.mat, pop=race.covs)
pc.pval.mat <- cbind(pc.pval.mat, pop=race.pval)

batch <- phen.cd4$Batch
batch.covs <- c()
batch.pval <- c()
for (sv in seq(ncol(cd4.svs))) {
	x <- lm(cd4.svs[,sv] ~ batch)
	batch.covs <- c(batch.covs, summary.lm(x)$r.squared)
	batch.pval <- c(batch.pval, summary.lm(x)$coefficients[2,4])
}
cov.mat <- cbind(cov.mat, batch=batch.covs)
cov.pval.mat <- cbind(cov.pval.mat, batch=batch.pval)

batch.covs <- c()
batch.pval <- c()
for (pc in seq(ncol(cd4.pcs))) {
	x <- lm(cd4.pcs[,pc] ~ batch)
	batch.covs <- c(batch.covs, summary.lm(x)$r.squared)
	batch.pval <- c(batch.pval, summary.lm(x)$coefficients[2,4])
}
pc.cov.mat <- cbind(pc.cov.mat, batch=batch.covs)
pc.pval.mat <- cbind(pc.pval.mat, batch=batch.pval)

frozen <- as.numeric(phen.cd4$Frozen)
frozen.covs <- c()
frozen.pval <- c()
for (sv in seq(ncol(cd4.svs))) {
	x <- lm(cd4.svs[,sv] ~ frozen)
	frozen.covs <- c(frozen.covs, summary.lm(x)$r.squared)
	frozen.pval <- c(frozen.pval, summary.lm(x)$coefficients[2,4])
}
cov.mat <- cbind(cov.mat, frozen=frozen.covs)
cov.pval.mat <- cbind(cov.pval.mat, frozen=frozen.pval)

frozen.covs <- c()
frozen.pval <- c()
for (pc in seq(ncol(cd4.pcs))) {
	x <- lm(cd4.pcs[,pc] ~ frozen)
	frozen.covs <- c(frozen.covs, summary.lm(x)$r.squared)
	frozen.pval <- c(frozen.pval, summary.lm(x)$coefficients[2,4])
}
pc.cov.mat <- cbind(pc.cov.mat, frozen=frozen.covs)
pc.pval.mat <- cbind(pc.pval.mat, frozen=frozen.covs)
#pdf(file='/scratch/t.cczysz/cd4_heatmap.pdf')
pdf(file='/group/stranger-lab/czysz/ImmVar/plots/cd4_pcavscovs.pdf', width=9, height=9)
#heatmap.2(pc.cov.mat[,-4],trace='none')
heatmap.2(cov.mat[,-4],cellnote=signif(cov.pval.mat[,-4],2), notecex=1, notecol='black',col=heat.colors(20)[seq(20,1)],trace='none', Rowv=NULL, key.xlab='R^2', margins=c(7,7))#, main='CD4+ SV vs Covariates')
heatmap.2(pc.cov.mat[,-4],cellnote=signif(pc.pval.mat[,-4],2), notecex=1, notecol='black',col=heat.colors(20)[seq(20,1)],trace='none', Rowv=NULL, key.xlab='R^2', margins=c(7,7))#, main='CD4+ PCs vs Covariates')
dev.off()

rm(cov.mat)
rm(cov.pval.mat)
rm(pc.cov.mat)
rm(pc.pval.mat)

covs <- c(7,10,11,12,13,18,19)
cov.mat <- matrix(nrow=ncol(cd14.svs), ncol=length(covs))
colnames(cov.mat) <- names(meta.cd14.tmp)[covs]
cov.pval.mat <- matrix(nrow=ncol(cd14.svs), ncol=length(covs))
colnames(cov.pval.mat) <- names(meta.cd14.tmp)[covs]

# Build matrix of individual x covariate
pc.cov.mat <- matrix(nrow=20, ncol=length(covs))
colnames(pc.cov.mat) <- names(meta.cd14.tmp)[covs]
pc.pval.mat <- matrix(nrow=20, ncol=length(covs))
colnames(pc.pval.mat) <- names(meta.cd14.tmp)[covs]

if (T) {
for (cov in seq(length(covs))) {
	for (sv in seq(ncol(cd14.svs))) {
		x <- lm(cd14.svs[, sv] ~ meta.cd14.tmp[, covs[cov]])
		cov.mat[sv, cov] <- summary.lm(x)$r.squared
		cov.pval.mat[sv, cov] <- summary.lm(x)$coefficients[2,4]	
		}
	for (pc in seq(ncol(cd14.pcs))) {
		x <- lm(cd14.pcs[, pc] ~ meta.cd14.tmp[, covs[cov]])
		pc.cov.mat[pc, cov] <- summary.lm(x)$r.squared
		pc.pval.mat[pc, cov] <- summary.lm(x)$coefficients[2,4]	
		}
}
}

pop <- as.character(phen.cd14$Race)
pop <- replace(pop, pop=='Caucasian', 0)
pop <- replace(pop, pop=='African-American', 1)
pop <- replace(pop, pop=='Asian', 2)
pop <- as.numeric(pop)
race.covs <- c()
race.pval <- c()
for (sv in seq(ncol(cd14.svs))) {
	x <- lm(cd14.svs[,sv] ~ pop)
	race.covs <- c(race.covs, summary.lm(x)$r.squared)
	race.pval <- c(race.pval, summary.lm(x)$coefficients[2,4])
}
cov.mat <- cbind(cov.mat, pop=race.covs)
cov.pval.mat <- cbind(cov.pval.mat, pop=race.pval)

race.covs <- c()
race.pval <- c()
for (pc in seq(ncol(cd14.pcs))) {
	x <- lm(cd14.pcs[,pc] ~ pop)
	race.covs <- c(race.covs, summary.lm(x)$r.squared)
	race.pval <- c(race.pval, summary.lm(x)$coefficients[2,4])
}
pc.cov.mat <- cbind(pc.cov.mat, pop=race.covs)
pc.pval.mat <- cbind(pc.pval.mat, pop=race.pval)

batch <- phen.cd14$Batch
batch.covs <- c()
batch.pval <- c()
for (sv in seq(ncol(cd14.svs))) {
	x <- lm(cd14.svs[,sv] ~ batch)
	batch.covs <- c(batch.covs, summary.lm(x)$r.squared)
	batch.pval <- c(batch.pval, summary.lm(x)$coefficients[2,4])
}
cov.mat <- cbind(cov.mat, batch=batch.covs)
cov.pval.mat <- cbind(cov.pval.mat, batch=batch.pval)

batch.covs <- c()
batch.pval <- c()
for (pc in seq(ncol(cd14.pcs))) {
	x <- lm(cd14.pcs[,pc] ~ batch)
	batch.covs <- c(batch.covs, summary.lm(x)$r.squared)
	batch.pval <- c(batch.pval, summary.lm(x)$coefficients[2,4])
}
pc.cov.mat <- cbind(pc.cov.mat, batch=batch.covs)
pc.pval.mat <- cbind(pc.pval.mat, batch=batch.pval)

frozen <- as.numeric(phen.cd14$Frozen)
frozen.covs <- c()
frozen.pval <- c()
for (sv in seq(ncol(cd14.svs))) {
	x <- lm(cd14.svs[,sv] ~ frozen)
	frozen.covs <- c(frozen.covs, summary.lm(x)$r.squared)
	frozen.pval <- c(frozen.pval, summary.lm(x)$coefficients[2,4])
}
cov.mat <- cbind(cov.mat, frozen=frozen.covs)
cov.pval.mat <- cbind(cov.pval.mat, frozen=frozen.pval)

frozen.covs <- c()
frozen.pval <- c()
for (pc in seq(ncol(cd14.pcs))) {
	x <- lm(cd14.pcs[,pc] ~ frozen)
	frozen.covs <- c(frozen.covs, summary.lm(x)$r.squared)
	frozen.pval <- c(frozen.pval, summary.lm(x)$coefficients[2,4])
}
pc.cov.mat <- cbind(pc.cov.mat, frozen=frozen.covs)
pc.pval.mat <- cbind(pc.pval.mat, frozen=frozen.covs)

sex <- as.character(phen.cd14$Sex)
sex <- as.numeric(sex=='Male')
sex.covs <- c()
sex.pval <- c()
for (sv in seq(ncol(cd14.svs))) {
	x <- lm(cd14.svs[,sv] ~ sex)
	sex.covs <- c(sex.covs, summary.lm(x)$r.squared)
	sex.pval <- c(sex.pval, summary.lm(x)$coefficients[2,4])
}
cov.mat <- cbind(cov.mat, sex=sex.covs)
cov.pval.mat <- cbind(cov.pval.mat, sex=sex.pval)

sex.covs <- c()
sex.pval <- c()
for (pc in seq(ncol(cd14.pcs))) {
	x <- lm(cd14.pcs[,pc] ~ sex)
	sex.covs <- c(sex.covs, summary.lm(x)$r.squared)
	sex.pval <- c(sex.pval, summary.lm(x)$coefficients[2,4])
}
pc.cov.mat <- cbind(pc.cov.mat, sex=sex.covs)
pc.pval.mat <- cbind(pc.pval.mat, sex=sex.covs)
#pdf(file='/scratch/t.cczysz/cd14_heatmap.pdf')
pdf(file='/group/stranger-lab/czysz/ImmVar/plots/cd14_pcavscovs.pdf', width=9, height=9)
#heatmap.2(pc.cov.mat[,-4],trace='none')
heatmap.2(cov.mat[,-4],cellnote=signif(cov.pval.mat[,-4],2), notecex=1, notecol='black',col=heat.colors(20)[seq(20,1)],trace='none', Rowv=NULL, key.xlab='R^2', margins=c(7,7))#, main='CD14+ SV vs Covariates')
heatmap.2(pc.cov.mat[,-4],cellnote=signif(pc.pval.mat[,-4],2), notecex=1, notecol='black',col=heat.colors(20)[seq(20,1)],trace='none', Rowv=NULL, key.xlab='R^2', margins=c(7,7))#, main='CD14+ PCs vs Covariates')
dev.off()
