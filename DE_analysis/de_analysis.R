library(oligo)
library(limma)
library(ggplot2)
library(peer)

# Populations: Caucasian, African-American, Asian
# Cell types: CD14, CD4
# Imported as "exp_genes": GxN matrix of normalized expression values
LoadData <- function(population=population,cell.type=cell.type) {
	data.dir <- "/group/stranger-lab/moliva/ImmVar/Robjects/"
	file.path <- paste(data.dir,"exp_genes.",cell.type,".",population,".Robj",sep="")
	load(file=file.path)
}

RunPeer <- function(expression, k=20, covs) {
	model = PEER()
	
	# Expression must be in NxG. N number of samples, G number of genes
	PEER_setPhenoMean(model,t(expression))
	PEER_setNk(model,k)
	#PEER_setAdd_mean(model, TRUE)
	PEER_setCovariates(model,as.matrix(covs)+1)
	PEER_update(model)
	peer.factors = PEER_getX(model)
	return(peer.factors)
}

MakeResiduals <- function(input.row,peer.factors) {
	fit <- lm(input.row ~ peer.factors[, -1] - 1)
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

	save.path <- "/group/stranger-lab/czysz/ImmVar"
	save.file.name <- paste("de_genes",cell.type,population,"txt",sep=".")
	write.table(topTable(eb.fit,number=100),file=save.file.name)
	#write.fit(eb.fit,file=paste(save.path,save.file.name,sep="/"),adjust="BH")

	# DE_probesets <- unique(merge_probes_DF[merge_probes_DF$gene_ensembl%in%top.de.genes, ]$fsetid)

	volcano.file <- paste("volcano",population,cell.type,"pdf",sep=".")
	pdf(file=paste(save.path,volcano.file,sep="/"))
	volcanoplot(eb.fit)
	#g = ggplot(data=data.frame(eb.fit),aes(x=coefficients,y=lods)) 
	#g + geom_point() + xlab("fold change") + ylab("log odds")
	dev.off() 
	
	if (F) {

	# Plot p.value here
	pdf(file = '/home/t.cri.cczysz/qqplot.pdf')
	x <- hist(seq(1000),plot=F,freq=F)
	dev.off()
	}

	de.expr.file <- paste("de_expression",population,cell.type,"pdf",sep=".")
	pdf(file=paste(save.path,de.expr.file,sep="/"))
	for (set in top.de.genes) {
		plot(density(expr.residuals[set, !!sex]),col='blue',xlim=c(-10,10),ylim=c(0,2),main=set)
		lines(density(expr.residuals[set, !sex]),col='red')
	}
	dev.off()
}

for (population in c("Caucasian","African-American","Asian")) {
for (cell.type in c("CD14","CD4")) {
	population <<- population
	cell.type <<- cell.type
	load('/group/stranger-lab/moliva/ImmVar/Robjects/phen.Robj')
	
	if (cell.type == "CD14") phen.cell.type <- "CD14+16-Mono"
	else phen.cell.type <- "CD4TNve"

	phen <- phen[phen$CellType == phen.cell.type, ]

	data.dir <- "/group/stranger-lab/moliva/ImmVar/Robjects/"
	file.path <- paste(data.dir,"exp_genes.",cell.type,".",population,".Robj",sep="")
	load(file=file.path)

	sex <- as.numeric(phen[phen$Race == population, ]$Sex == 'Male')
	peer.factors <- RunPeer(exp_genes,k=20,sex)

	expr.residuals <- apply(exp_genes, 1, MakeResiduals, peer.factors=peer.factors)
	expr.residuals <- t(expr.residuals)

	save.file.name <- paste("de_genes",cell.type,population,"txt",sep=".")
	# save(expr.residuals, file = paste("/group/stranger-lab/immvar_data/",save.file.name,sep=""))

	eb.fit <- PerformDEAnalysis(expr.residuals, sex)

	AnalyzeFit(eb.fit, expr.residuals, sex)

	print(topTable(eb.fit, number=10))
	}
}
