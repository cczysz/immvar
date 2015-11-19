library(oligo)
library(limma)
#library(ggplot2)
library(sva)
library(genefilter)
library(gtools)
library(ggplot2)

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
	fit <- lm(input.row ~ peer.factors[, -1] )
	residuals(fit)
}

PerformDEAnalysis <- function(expr,samples, cau, afr) {
	mod = model.matrix(~0+as.factor(samples) + cau + afr)
        colnames(mod) <- c('Female', 'Male', 'Cau', 'Afr')
        #mod0 = model.matrix(~0+rep(1,ncol(expr)))
        mod0 = model.matrix(~0+ cau + afr)

        svobj <- sva(expr, mod, mod0)
        modSv <- cbind(mod, svobj$sv)

        fit <- lmFit(expr, modSv)
        contrast.matrix <- c(-1,1,0,0,rep(0, svobj$n.sv))
        contrast.fit <- contrasts.fit(fit, contrast.matrix)
        fit <- eBayes(contrast.fit)
	return(fit)
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

}
FTest <- function(fit, exp_genes, expr.residuals, sex) {
	#ttable <- topTable(fit,number=Inf)
	males <- !!sex
	females <- !sex

	# Expression-Based Filtering
	females_filt=kOverA(A=quantile(exp_genes,probs = 0.1),k = round(length(!sex))/3)
	males_filt=kOverA(A=quantile(exp_genes,probs = 0.1),k = round(length(!!sex))/3)
	
	male_filt_genes <- genefilter(exp_genes[,males],filterfun(males_filt))
	female_filt_genes <- genefilter(exp_genes[,females],filterfun(females_filt))

	filt_genes_AND <- intersect(male_filt_genes[!male_filt_genes],female_filt_genes[!female_filt_genes])
	filt_genes_OR <- union(male_filt_genes[!male_filt_genes],female_filt_genes[!female_filt_genes])

	# Variance-based filtering
	var.filt.males <- varFilter(exp_genes[,males],var.cutoff=0.1) # Default cutoff of 0.5 removes 50% of genes
	var.filt.females <- varFilter(exp_genes[,females],var.cutoff=0.1)
	var.filt <- intersect(rownames(var.filt.males), rownames(var.filt.females))

	genes_to_filter <- filt_genes_OR[!(filt_genes_OR%in%var.filt)]
	residuals.filtered <- expr.residuals[!(rownames(expr.residuals)%in%genes_to_filter),]

	# Filtering
	#filter.val <- quantile(as.numeric(expr.residuals),probs=c(0.1))
	#f1 <- pOverA(p=0.1, A=filter.val)
	#f1 <- kOverA(k=10, A=filter.val)
	#ffun <- filterfun(f1)
	#filt.exp <- genefilter(expr.residuals, ffun)

	ftest.results <- data.frame()
	for (gene in rownames(residuals.filtered)) {
		gene.exp <- residuals.filtered[gene, ]
		# Filter for no expression in either males or females
		f.test <- var.test(gene.exp[!!sex], gene.exp[!sex])
		#f.pval <- f.test$p.value
		dat.f <- data.frame(f=f.test$statistic, 
			p.val=f.test$p.value, 
			ratio=f.test$estimate[[1]],
			row.names=f.test$data.name)
		rownames(dat.f) <- gene
		ftest.results <- rbind(ftest.results, dat.f)
	}
	p.adj <- p.adjust(ftest.results$p.val, method="fdr")
	ftest.results <- cbind(ftest.results, p.adj)
	sig.results <- ftest.results[ftest.results$p.adj < 0.05, ]
	f.name <- paste("ftest",cell.type,population,"Robj",sep=".")
	save(ftest.results,sig.results,file=paste('/group/stranger-lab/immvar_data/',f.name,sep=""))

	f.name <- paste("sig.ftest",cell.type,population,"pdf",sep=".")
	#sig.results <- ftest.results[ftest.results$p.adj < 0.05, ]
	pdf(file=paste('/group/stranger-lab/czysz/ImmVar/plots/',f.name,sep=""))
	for (gene in rownames(sig.results)) {
		male_exp <- density(expr.residuals[gene, males])
		female_exp <- density(expr.residuals[gene, females])
		x_min <- min(min(male_exp$x), min(female_exp$x))
		x_max <- max(max(male_exp$x), max(female_exp$x))
		ylim <- max(max(male_exp$y), max(female_exp$y))
		plot(density(expr.residuals[gene,males]),col="blue",xlim=c(x_min,x_max),ylim=c(0,ylim),main=gene)
		lines(density(expr.residuals[gene,females]),col="red")
	}
	dev.off()
	
}

#for (cell.type in c("CD4","CD14")) {
for (cell.type in c("CD14")) {
	#population <<- population
	cell.type <<- cell.type
	load('/group/stranger-lab/moliva/ImmVar/Robjects/phen.Robj')
	load('/group/stranger-lab/immvar_rep/mappings/annot.Robj')
	#load('/group/stranger-lab/moliva/ImmVar/probes_mapping/Robjects/merge_probes_DF.Robj')
	#annots <- read.table(file='/group/stranger-lab/moliva/ImmVar/probes_mapping/annotations/gene_annot', header=F, row.names=1)
	#colnames(annots) <- c('symbol', 'chr', 'start', 'stop', 'strand')
	
	if (cell.type == "CD14") phen.cell.type <- "CD14+16-Mono"
	else phen.cell.type <- "CD4TNve"

	phen <- phen[phen$CellType == phen.cell.type, ]

	#data.dir <- "/group/stranger-lab/moliva/ImmVar/Robjects/"
	#load(file=paste(data.dir,paste("exp_genes",cell.type,population,"Robj",sep="."),sep=""))
	load(file=paste('/group/stranger-lab/moliva/ImmVar/Robjects/', paste(cell.type,'.joint.norm.exp_genes.Robj',sep=''),sep=''))
	if (cell.type == "CD4") { exp_genes <- exp_genes.cd4.joint.norm 
	new_ids <- c()
	for (id in colnames(exp_genes)) {
		new_ids <- c(new_ids, unlist(strsplit(id, split=":"))[3])	
	}	
	colnames(exp_genes) <- new_ids
		} else { exp_genes <- exp_genes.cd14.joint.norm 
	new_ids <- c()
	for (id in colnames(exp_genes)) {
		new_ids <- c(new_ids, unlist(strsplit(id, split=":"))[3])	
	}	
	colnames(exp_genes) <- new_ids
	}


	rownames(phen) <- phen$ImmVarID2

	sex <- phen[colnames(exp_genes),"Sex"]
	
	cau <- as.numeric(phen[colnames(exp_genes), "Race"]=="Caucasian")
	afr <- as.numeric(phen[colnames(exp_genes), "Race"]=="African-American")
	eb.fit <- PerformDEAnalysis(exp_genes, as.numeric(sex=="Male"), cau, afr)
	eb.fit$N <- as.numeric(rep(ncol(exp_genes), nrow(eb.fit)))
	eb.fit$chr <- as.character(annots[rownames(eb.fit), "chr"])
	eb.fit$symbol <- as.character(annots[rownames(eb.fit), "symbol_id"])
	fit.save.name <- paste("fit.joint",cell.type,"Robj",sep=".")
	#fit.save.name <- paste("fit",cell.type,population,"Robj",sep=".")
	save(eb.fit, file=paste("/group/stranger-lab/immvar_data/",fit.save.name,sep=""))
	
	pdf(file=paste('/group/stranger-lab/czysz/ImmVar/plots/', paste(cell.type,'_volcano.pdf',sep=''),sep=''))
	title <- paste(cell.type, sep=" ")
	volcanoplot(eb.fit, cex=0.5, names=eb.fit$symbol, highlight=50, main=paste(title, "all genes"))
	#g <- ggplot(data=eb.fit, aes(x=coefficients , y=p.value))
	#fit.sex <- eb.fit[(eb.fit$chr == 'chrY' | eb.fit$chr=='chrX'),] 
	#fit.auto <- eb.fit[!(eb.fit$chr == 'chrY' | eb.fit$chr=='chrX'),] 
	#volcanoplot(fit.sex, cex=0.5, names=fit.sex$symbol, highlight=50, main=paste(title, "X and Y genes"))
	#volcanoplot(fit.auto, cex=0.5, names=fit.auto$symbol, highlight=50, main=paste(title, "Autosomal genes"))
	dev.off()
	
	#AnalyzeFit(eb.fit, exp_genes, sex)
	#FTest(eb.fit,exp_genes,expr.residuals, sex)

	}
#}
