library(oligo)
library(limma)
library(ggplot2)

population <<- "Caucasian"
cell.type <<- "CD14"
load('/group/stranger-lab/moliva/ImmVar/Robjects/phen.Robj')
phen <- phen[phen$CellType == 'CD14+16-Mono', ]

# Populations: Caucasian, African-American, Asian
# Cell types: CD14, CD4
# Imported as "exp_genes": GxN matrix of normalized expression values
LoadData <- function(population=population,cell.type=cell.type) {
	data.dir <- "/group/stranger-lab/moliva/ImmVar/Robjects/"
	file.path <- paste(data.dir,"exp_genes.",cell.type,".",population,".Robj",sep="")
	load(file=file.path)
}

MakeResiduals <- function(input.row,peer.factors) {
	residuals(lm(input.row ~ 0 + peer.factors[, -1]))
}

PerformDEAnalysis <- function(expr,samples) {
	design <- model.matrix(~0 + samples)
	colnames(design) <- c("Female","Male")

	fit <- lmFit(expr, design)

	contrast.matrix <- makeContrasts(mf = Male - Female, levels=design)
	contrast.fit <- contrasts.fit(fit, contrast.matrix)
	eb.fit <- eBayes(contrast.fit, robust=F)
}

AnalyzeFit <- function(eb.fit) {
	top.de.genes <- row.names(topTable(eb.fit, number=100))

	DE_probesets <- unique(merge_probes_DF[merge_probes_DF$gene_ensembl%in%top.de.genes, ]$fsetid)

	g = ggplot(data=data.frame(eb.fit),aes(x=coefficients,y=lods)) 
	pdf('/home/t.cri.cczysz/volcano.pdf')
	g + geom_point() + xlab("fold change") + ylab("log odds")
	dev.off() 

	pdf(file='/home/t.cri.cczysz/de_probesets.pdf')
	for (set in top.de.genes) {
		plot(density(expr.residuals[set, !!sex]),col='blue',xlim=c(-1,1),ylim=c(0,1))
		lines(density(expr.residuals[set, !sex]),col='red')
	}
	dev.off()
}

#LoadData()
data.dir <- "/group/stranger-lab/moliva/ImmVar/Robjects/"
file.path <- paste(data.dir,"exp_genes.",cell.type,".",population,".Robj",sep="")
load(file=file.path)

source('/home/t.cri.cczysz/thesis/scripts/peer.R')
sex <- as.numeric(phen[phen$Race == 'Caucasian', ]$Sex == 'Male')
peer.factors <- RunPeer(exp_genes,k=20,sex)

expr.residuals <- apply(exp_genes, 1, MakeResiduals, peer.factors=peer.factors)
expr.residuals <- t(expr.residuals)

eb.fit <- PerformDEAnalysis(expr.residuals,as.factor(sex))

AnalyzeFit(eb.fit)
