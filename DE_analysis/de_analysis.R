library(oligo)
library(limma)
library(ggplot2)

# Populations: Caucasian, African-American, Asian
# Cell types: CD14, CD4
# Imported as "exp_genes": GxN matrix of normalized expression values
LoadData <- function(population=population,cell.type=cell.type) {
	data.dir <- "/group/stranger-lab/moliva/ImmVar/Robjects/"
	file.path <- paste(data.dir,"exp_genes.",cell.type,".",population,".Robj",sep="")
	load(file=file.path)
}

MakeResiduals <- function(input.row,peer.factors) {
	fit <- lm(input.row ~ peer.factors[, -1:-2] - 1)
	residuals(fit)
}

PerformDEAnalysis <- function(expr,samples) {
	design <- model.matrix(~ as.factor(samples)-1)
	colnames(design) <- c("Female","Male")

	fit <- lmFit(expr, design)

	contrast.matrix <- makeContrasts(Male - Female, levels=design)
	contrast.fit <- contrasts.fit(fit, contrast.matrix)
	eb.fit <- eBayes(contrast.fit, robust=T)
}

AnalyzeFit <- function(eb.fit, expr.residuals, sex) {
	top.de.genes <- row.names(topTable(eb.fit, number=100))
	# load('/group/stranger-lab/moliva/ImmVar/probes_mapping/Robjects/merge_probes_DF.Robj')

	# DE_probesets <- unique(merge_probes_DF[merge_probes_DF$gene_ensembl%in%top.de.genes, ]$fsetid)

	#g = ggplot(data=data.frame(eb.fit),aes(x=coefficients,y=lods)) 
	pdf('/home/t.cri.cczysz/volcano.pdf')
	volcanoplot(eb.fit)
	#g + geom_point() + xlab("fold change") + ylab("log odds")
	dev.off() 
	if (F) {

# Plot p.value here
	pdf(file = '/home/t.cri.cczysz/qqplot.pdf')
	x <- hist(seq(1000),plot=F,freq=F)

	dev.off()
	}

	pdf(file='/home/t.cri.cczysz/de_probesets.pdf')
	for (set in top.de.genes) {
		plot(density(expr.residuals[set, !!sex]),col='blue',xlim=c(-1,1),ylim=c(0,3),main=set)
		lines(density(expr.residuals[set, !sex]),col='red')
	}
	dev.off()
}

#LoadData()
population <<- "Caucasian"
cell.type <<- "CD14"
load('/group/stranger-lab/moliva/ImmVar/Robjects/phen.Robj')
phen <- phen[phen$CellType == 'CD14+16-Mono', ]
data.dir <- "/group/stranger-lab/moliva/ImmVar/Robjects/"
file.path <- paste(data.dir,"exp_genes.",cell.type,".",population,".Robj",sep="")
load(file=file.path)

source('/home/t.cri.cczysz/thesis/scripts/peer.R')
sex <- as.numeric(phen[phen$Race == 'Caucasian', ]$Sex == 'Male')
peer.factors <- RunPeer(exp_genes,k=20,sex)

expr.residuals <- apply(exp_genes, 1, MakeResiduals, peer.factors=peer.factors)
expr.residuals <- t(expr.residuals)

eb.fit <- PerformDEAnalysis(expr.residuals, sex)

AnalyzeFit(eb.fit, expr.residuals, sex)

print(topTable(eb.fit, number=10))
save(expr.residuals,file='/group/stranger-lab/immvar_data/CD14.Caucasian.Residuals.Robj')
