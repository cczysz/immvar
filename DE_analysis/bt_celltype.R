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
	colnames(design) <- c("Female","Male")

	fit <- lmFit(expr, design)

	contrast.matrix <- makeContrasts(Male - Female, levels=design)
	contrast.fit <- contrasts.fit(fit, contrast.matrix)
	eb.fit <- eBayes(contrast.fit, robust=T)
}

AnalyzeFit <- function(eb.fit, expr.residuals, sex) {
	top.de.genes <- row.names(topTable(eb.fit, number=Inf,p.value=0.05))
	# load('/group/stranger-lab/moliva/ImmVar/probes_mapping/Robjects/merge_probes_DF.Robj')

	save.path <- "/group/stranger-lab/czysz/ImmVar"
	save.file.name <- paste("de_genes",cell.type,population,"txt",sep=".")
	write.table(topTable(eb.fit,number=Inf,p.value=0.05),file=paste(save.path,save.file.name,sep="/"))

	de.expr.file <- paste("de_expression",population,cell.type,"pdf",sep=".")
	pdf(file=paste(save.path,de.expr.file,sep="/"))
	for (set in top.de.genes) {
		#plot(density(expr.residuals[set, !!sex]),col='blue',xlim=c(-10,10),ylim=c(0,2),main=set)
		male <- expr.residuals[set, !!sex]
		female <- expr.residuals[set, !sex]

		dm <- density(male)
		df <- density(female)
		xmin <- floor(min(min(dm$x),min(df$x)))
		xmax <- ceiling(max(max(dm$x),max(df$x)))
		ymax <- ceiling(max(max(dm$y),max(df$y)))

		pval <- wilcox.test(male,female)$p.value
		plot(dm,col='blue',main=paste(set,"\n","Wilcox Test: ",pval,sep=""),xlim=c(xmin,xmax),ylim=c(0,ymax))
		lines(df,col='red')
	}
	dev.off()

	ttable <- topTable(eb.fit,number=Inf)
	ftest.results <- c()
	for (gene in rownames(ttable)) {
		f.test <- var.test(expr.residuals[gene,!!sex],expr.residuals[gene,!sex])
		f.pval <- f.test$p.value
		ftest.results <- rbind(ftest.results, c(gene,f.test,f.pval))
	}
	write.table(ftest.results,file="/group/stranger-lab/czysz/ImmVar/ftestresults.txt")
	pdf(file=paste(save.path,"plots",paste("ftest",population,cell.type,"pdf",sep="."),sep="/"))
	plot(density(f.pval),main=paste("F test",cell.type,population,))
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
for (cell.type in c("CD14","CD4")) {
	population <<- population
	cell.type <<- cell.type
	load('/group/stranger-lab/moliva/ImmVar/Robjects/phen.Robj')
	
	if (cell.type == "CD14") phen.cell.type <- "CD14+16-Mono"
	else phen.cell.type <- "CD4TNve"

	phen <- phen[phen$CellType == phen.cell.type, ]

	#data.dir <- "/group/stranger-lab/moliva/ImmVar/Robjects/"
	data.dir <- "/group/stranger-lab/immvar_data/"
	res.file.name <- paste("residuals",cell.type,population,"Robj",sep=".")
	load(file = paste(data.dir,res.file.name,sep=""))
	fit.file.name <- paste("fit",cell.type,population,"Robj",sep=".")
	load(file = paste(data.dir,fit.file.name,sep=""))
	sex <- as.numeric(phen[phen$Race == population, ]$Sex == 'Male')

	if (F) {
		file.path <- paste(data.dir,"exp_genes.",cell.type,".",population,".Robj",sep="")
		load(file=file.path)

		sex <- as.numeric(phen[phen$Race == population, ]$Sex == 'Male')
		peer.factors <- RunPeer(exp_genes,k=20,sex)

		expr.residuals <- apply(exp_genes, 1, MakeResiduals, peer.factors=peer.factors)
		expr.residuals <- t(expr.residuals)

		save.file.name <- paste("de_genes",cell.type,population,"txt",sep=".")
		# save(expr.residuals, file = paste("/group/stranger-lab/immvar_data/",save.file.name,sep=""))
	

	eb.fit <- PerformDEAnalysis(expr.residuals, sex)

	fit.save.name <- paste("fit",cell.type,population,"Robj",sep=".")
	save(eb.fit, file=paste("/group/stranger-lab/immvar_data/",fit.save.name,sep=""))
	}

	#AnalyzeFit(eb.fit, expr.residuals, sex)
	
	FTest(eb.fit,expr.residuals, sex)

	}
}
