library(oligo)
library(limma)
#library(ggplot2)
library(peer)

# Populations: Caucasian, African-American, Asian
# Cell types: CD14, CD4
# Imported as "exp_genes": GxN matrix of normalized expression values

RunPeer <- function(expression, k=20, covs) {
	model = PEER()
	
	# Expression must be in NxG. N number of samples, G number of genes
	PEER_setPhenoMean(model,t(expression))
	PEER_setNk(model,k)
	#PEER_setCovariates(model,as.matrix(covs)+1)
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
	colnames(design) <- c("CD4","CD14")
	print(design)

	fit <- lmFit(expr, design)

	contrast.matrix <- makeContrasts(CD4 - CD14, levels=design)
	contrast.fit <- contrasts.fit(fit, contrast.matrix)
	eb.fit <- eBayes(contrast.fit, robust=T)
}

AnalyzeFit <- function(eb.fit, expr.residuals, type, cell) {
	top.de.genes <- row.names(topTable(eb.fit, number=Inf,p.value=0.05))

	save.path <- "/group/stranger-lab/czysz/ImmVar"
	#save.file.name <- paste("de_genes_bt", population,cell, "txt", sep=".")
	#write.table(topTable(eb.fit,number=Inf),file=paste(save.path,save.file.name,sep="/"))

	de.expr.file <- paste("de_expression_bt", population, cell, "pdf", sep=".")
	pdf(file=paste(save.path,de.expr.file,sep="/"))
	for (set in top.de.genes) {

		d.m <- density(expr.residuals[,!!type])
		d.f <- density(expr.residuals[,!type])
		xmin <- floor(min(min(d.m$x),min(d.f$x)))
		xmax <- ceiling(max(max(d.m$x),max(d.f$x)))
		ymax <- ceiling(max(max(d.m$y),max(d.f$y)))

		pval <- wilcox.test(expr.residuals[set,!!type], expr.residuals[set,!type])$p.value
		plot(dm,col='blue',main=paste(set,"\n","Wilcox Test: ",pval,sep=""),xlim=c(xmin,xmax),ylim=c(0,ymax))
		lines(df,col='red')
	}
	dev.off()
}
FTest <- function(fit, expr.residuals, sex) {
	ttable <- topTable(fit,number=Inf)
	ftest.results <- data.frame()
	for (gene in rownames(ttable)) {
		gene <- gene
		# Filter for no expression in either males or females
			# Get list from Meri
		f.test <- var.test(expr.residuals[gene,!!sex],expr.residuals[gene,!sex])
		#f.pval <- f.test$p.value
		dat.f <- data.frame(f=f.test$statistic, 
			p.val=f.test$p.value, 
			ratio=f.test$estimate[[1]],
			row.names=f.test$data.name)
		rownames(dat.f) <- gene
		ftest.results <- rbind(ftest.results, dat.f)
	}
	f.name <- paste("ftest",cell.type,population,"Robj",sep=".")
	save(ftest.results,file=paste('/group/stranger-lab/immvar_data/',f.name,sep=""))
}

for (population in c("Caucasian")){ #,"African-American","Asian")) {
	population <<- population
	load('/group/stranger-lab/moliva/ImmVar/Robjects/phen.Robj')
	
	cd14 <- "CD14+16-Mono"
	cd4 <- "CD4TNve"
	
	cd14.phen <- phen[phen$CellType == cd14 & phen$Race==population, ]
	cd4.phen <- phen[phen$CellType == cd4 & phen$Race==population, ]

	data.dir <- "/group/stranger-lab/immvar_data/"
	cd4.exp.file <- paste("exp_genes_intersect",'CD4',population,"Robj",sep='.')
	cd14.exp.file <- paste("exp_genes_intersect",'CD14',population,"Robj",sep='.')
	
	load(file=paste(data.dir,cd4.exp.file,sep=''))
	cd4.exp <- exp_genes
	load(file=paste(data.dir,cd14.exp.file,sep=''))
	cd14.exp <- exp_genes

	shared_genes <- intersect(rownames(cd4.exp),rownames(cd14.exp))
	cd4.exp <- cd4.exp[rownames(cd4.exp)%in%shared_genes,]
	cd14.exp <- cd14.exp[rownames(cd14.exp)%in%shared_genes,]

	cd4.sex <- cd4.phen[cd4.phen$ImmVarID2%in%colnames(cd4.exp),"Sex"]	
	cd14.sex <- cd14.phen[cd14.phen$ImmVarID2%in%colnames(cd14.exp),"Sex"]	
	
	cell.type <- c(rep(0, ncol(cd4.exp)), rep(1, ncol(cd14.exp)))

	cd4.peer.factors <- RunPeer(cd4.exp, k=20, cov=NULL)
	cd14.peer.factors <- RunPeer(cd14.exp, k=20, cov=NULL)
	
	cd4.exp.residuals <- apply(as.matrix(cd4.exp), 1, MakeResiduals, peer.factors=cd4.peer.factors)
	cd4.exp.residuals <- t(cd4.exp.residuals)
	
	cd14.exp.residuals <- apply(as.matrix(cd14.exp), 1, MakeResiduals, peer.factors=cd14.peer.factors)
	cd14.exp.residuals <- t(cd14.exp.residuals)
	
	cd4.eb.fit <- PerformDEAnalysis(cd4.exp.residuals, cd4.sex)
	cd14.eb.fit <- PerformDEAnalysis(cd14.exp.residuals, cd14.sex)

	AnalyzeFit(cd4.eb.fit, cd4.exp.residuals, cd4.sex=='Male',"CD4")
	AnalyzeFit(cd14.eb.fit, cd14.exp.residuals, cd14.sex=='Male', "CD14")
	
	#FTest(eb.fit, all.exp, cell.type)
}
