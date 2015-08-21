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
	colnames(design) <- c("CD14","CD4")

	fit <- lmFit(expr, design)

	contrast.matrix <- makeContrasts(Male - Female, levels=design)
	contrast.fit <- contrasts.fit(fit, contrast.matrix)
	eb.fit <- eBayes(contrast.fit, robust=T)
}

AnalyzeFit <- function(eb.fit, expr.residuals, type) {
	top.de.genes <- row.names(topTable(eb.fit, number=Inf,p.value=0.05))
	# load('/group/stranger-lab/moliva/ImmVar/probes_mapping/Robjects/merge_probes_DF.Robj')

	save.path <- "/group/stranger-lab/czysz/ImmVar"
	save.file.name <- paste("de_genes_bt", population, "txt", sep=".")
	write.table(topTable(eb.fit,number=Inf,p.value=0.05),file=paste(save.path,save.file.name,sep="/"))

	de.expr.file <- paste("de_expression_bt", population, "pdf", sep=".")
	pdf(file=paste(save.path,de.expr.file,sep="/"))
	for (set in top.de.genes) {
		#plot(density(expr.residuals[set, !!sex]),col='blue',xlim=c(-10,10),ylim=c(0,2),main=set)
		cd14.exp <- expr.residuals[set, !!type]
		cd4.exp <- expr.residuals[set, !type]

		d.14 <- density(cd14.exp)
		d.4 <- density(cd4.exp)
		xmin <- floor(min(min(d.14$x),min(d.4$x)))
		xmax <- ceiling(max(max(d.14$x),max(d.4$x)))
		ymax <- ceiling(max(max(d.14$y),max(d.4$y)))

		pval <- wilcox.test(cd14.exp, cd4.exp)$p.value
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

for (population in c("Caucasian","African-American","Asian")) {
	population <<- population
	load('/group/stranger-lab/moliva/ImmVar/Robjects/phen.Robj')
	
	cd14 <- "CD14+16-Mono"
	cd4 <- "CD4TNve"
	
	cd14.phen <- phen[phen$CellType == cd14, ]
	cd4.phen <- phen[phen$CellType == cd4, ]

	data.dir <- "/group/stranger-lab/immvar_data/"
	cd4.res.file.name <- paste("residuals",cd4,population,"Robj",sep=".")
	load(file = paste(data.dir,res.file.name,sep=""))
	cd4.res <- expr.residuals

	cd14.res.file.name <- paste("residuals",cd14,population,"Robj",sep=".")
	load(file = paste(data.dir,cd14.res.file.name,sep=""))
	cd14.res <- expr.residuals

	cell.type <- c(rep(0, nrow(cd4.res)), rep(1, nrow(cd14.res)))

	shared.genes <- intersect(rownames(cd4.res), rownames(cd14.res))
	
	all.exp <- numeric()
	for (gene in shared.genes){
		all.exp <- rbind(all.exp, cbind(cd4.res[gene, ], cd14.res[gene, ]))
	}
	eb.fit <- PerformDEAnalysis(all.exp, cell.type)

	fit.save.name <- paste("bt_cell_fit",population,"Robj",sep=".")
	save(eb.fit, file=paste("/group/stranger-lab/immvar_data/",fit.save.name,sep=""))
	
	AnalyzeFit(eb.fit, all.exp, cell.type)
	
	#FTest(eb.fit, all.exp, cell.type)
}
