library(oligo)
setwd('/group/stranger-lab/immvar_data/')

# Set to true or false depending on if these data are already loaded or not
if (T) {
############# CD14 Data
load('/group/stranger-lab/immvar_rep/GSE56034.Robj')
cd14.eset <- log2(exprs(gse56034[[1]]))
cd14.pdat <- pData(gse56034[[1]]) # Population not given. Must overlap with phen

EU_PC.cd14 <- as.matrix(read.table(file='/group/stranger-lab/immvar_rep/GSE56034_GSM.ImmVarCD14.EU.PC20.txt', header=T, row.names=1))
EA_PC.cd14 <- as.matrix(read.table(file='/group/stranger-lab/immvar_rep/GSE56034_GSM.ImmVarCD14.EA.PC10.txt', header=T, row.names=1))
AA_PC.cd14 <- as.matrix(read.table(file='/group/stranger-lab/immvar_rep/GSE56034_GSM.ImmVarCD14.AA.PC14.txt', header=T, row.names=1))

cau_ids.cd14 <- colnames(EU_PC.cd14)
asn_ids.cd14 <- colnames(EA_PC.cd14)
afr_ids.cd14 <- colnames(AA_PC.cd14)

############## CD4 Data
load('/group/stranger-lab/immvar_rep/GSE56033.Robj')
cd4.eset <- log2(exprs(gse56033[[1]]))
cd4.pdat <- pData(gse56033[[1]])

EU_PC.cd4 <- as.matrix(read.table(file='/group/stranger-lab/immvar_rep/GSE56033_GSM.ImmVarCD4.EU.PC20.txt', header=T, row.names=1))
AA_PC.cd4 <- as.matrix(read.table(file='/group/stranger-lab/immvar_rep/GSE56033_GSM.ImmVarCD4.AA.PC12.txt', header=T, row.names=1))
EA_PC.cd4 <- as.matrix(read.table(file='/group/stranger-lab/immvar_rep/GSE56033_GSM.ImmVarCD4.EA.PC12.txt', header=T, row.names=1))

cau_ids.cd4 <- colnames(EU_PC.cd4)
asn_ids.cd4 <- colnames(EA_PC.cd4)
afr_ids.cd4 <- colnames(AA_PC.cd4)

#load('/group/stranger-lab/czysz/ImmVar/cd4.probe2ens.Robj')
#load('/group/stranger-lab/czysz/ImmVar/cd14.probe2ens.Robj')

load('/group/stranger-lab/immvar_data/stats.Robj')
#rownames(cd14.eset) <- cd14.raj.ensembl
#rownames(cd4.eset) <- cd4.raj.ensembl

orig_ids <- rownames(cd4.eset)
rownames(cd14.eset) <- stats.all$ensembl_gene
rownames(cd4.eset) <- stats.all$ensembl_gene
}

cd14.ids <- list(cau_ids.cd14, afr_ids.cd14, asn_ids.cd14)
cd4.ids <- list(cau_ids.cd4, afr_ids.cd4, asn_ids.cd4)
cd4.pcs <- list(EU_PC.cd4, AA_PC.cd4, EA_PC.cd4)
cd14.pcs <- list(EU_PC.cd14, AA_PC.cd14, EA_PC.cd14)
pop <- c('Caucasian', 'African-American', 'Asian')

if (F) {
pdf(file='/group/stranger-lab/czysz/ImmVar/plots/all.comparison.pdf')
plot(0,0,xlim=c(-5,5), ylim=c(0,2))
cols <- rainbow(3)
for (i in seq(3)) {
	#lines(density(na.omit(as.numeric(cd4.eset[, colnames(cd4.eset)%in%cd4.ids[[i]]]))), col=cols[i], lty=1, xlim=c(0,20))
	#lines(density(na.omit(as.numeric(cd14.eset[, colnames(cd14.eset)%in%cd14.ids[[i]]]))), col=cols[i], lty=2)
	lines(density(na.omit(as.numeric(cd4.pcs[[i]][, colnames(cd4.pcs[[i]])%in%cd4.ids[[i]]]))), col=cols[i], lty=1, xlim=c(0,20))
	lines(density(na.omit(as.numeric(cd14.pcs[[i]][, colnames(cd14.pcs[[i]])%in%cd14.ids[[i]]]))), col=cols[i], lty=2, xlim=c(0,20))
}
dev.off()
}

pdf(file='/group/stranger-lab/czysz/ImmVar/plots/cd4.comparison.pdf')
for (i in seq(1)) {
	cd4.raj.mean <- apply(na.omit(cd4.eset[!duplicated(rownames(cd4.eset)), colnames(cd4.eset)%in%cd4.ids[[i]]]), 1 ,mean)
	cd4.raj.median <- apply(na.omit(cd4.eset[!duplicated(rownames(cd4.eset)), colnames(cd4.eset)%in%cd4.ids[[i]]]), 1 ,median)
	cd4.raj.var <- apply(na.omit(cd4.eset[!duplicated(rownames(cd4.eset)), colnames(cd4.eset)%in%cd4.ids[[i]]]), 1 , var)
	load(file=paste('/group/stranger-lab/immvar_data/', paste('exp_genes.CD4',pop[i],'Robj',sep='.'),sep=''))
	
	cd4.cc.mean <- apply(exp_genes, 1 ,mean)	
	cd4.cc.median <- apply(exp_genes, 1 ,median)	
	cd4.cc.var <- apply(exp_genes, 1, var)
	shared.genes <- intersect(names(cd4.cc.mean), names(cd4.raj.mean))

	# Use only genes in both datasets	
	cc.mean.shared <- cd4.cc.mean[names(cd4.cc.mean)%in%shared.genes]
	raj.mean.shared <- cd4.raj.mean[names(cd4.raj.mean)%in%shared.genes]
	
	# Change to apply eventually
	x<-c()
	for (gene in names(raj.mean.shared)) { x <- rbind(x, matrix(c(cc.mean.shared[gene], raj.mean.shared[gene]), ncol=2)) }
	rownames(x) <- names(raj.mean.shared)
	colnames(x) <- c("Renormalized", "Downloaded")
	cd4.diff.x <- abs(x[,1] - x[,2])
	top.cd4.diff.x <- cd4.diff.x[cd4.diff.x >= quantile(cd4.diff.x, probs=c(0.99))]
	
	# Plot the distributions of the top few differential (mean-wise) genes
	#for (gene in names(tail(sort(top.cd4.diff.x)))) {
		#plot(density(exp_genes[as.character(gene),]), col='red', main=as.character(gene),xlim=c(0,15))
		#lines(density(cd4.eset[as.character(gene), colnames(cd4.eset)%in%cd4.ids[[i]]]), col='blue')
	#}

	# Use only genes in both datasets	
	cc.median.shared <- cd4.cc.median[names(cd4.cc.median)%in%shared.genes]
	raj.median.shared <- cd4.raj.median[names(cd4.raj.median)%in%shared.genes]
	
	# Change to apply eventually
	y<-c()
	for (gene in names(raj.median.shared)) { y <- rbind(y, matrix(c(cc.median.shared[gene], raj.median.shared[gene]), ncol=2)) }
	rownames(y) <- names(raj.median.shared)
	colnames(y) <- c("Renormalized", "Downloaded")
	cd4.diff.y <- abs(y[,1] - y[,2])
	top.cd4.diff.y <- cd4.diff.y[cd4.diff.y >= quantile(cd4.diff.y, probs=c(0.99))]

	# Use only genes in both datasets	
	cc.var.shared <- cd4.cc.var[names(cd4.cc.var)%in%shared.genes]
	raj.var.shared <- cd4.raj.var[names(cd4.raj.var)%in%shared.genes]
	
	# Change to apply eventually
	z<-c()
	for (gene in names(raj.var.shared)) { z <- rbind(z, matrix(c(cc.var.shared[gene], raj.var.shared[gene]), ncol=2)) }
	rownames(z) <- names(raj.var.shared)
	colnames(z) <- c("Renormalized", "Downloaded")
	cd4.diff.z <- abs(z[,1] - z[,2])
	top.cd4.diff.z <- cd4.diff.z[cd4.diff.z >= quantile(cd4.diff.z, probs=c(0.99))]

	print(cor(x))
	print(cor(y))
	print(cor(z))
	
	plot(x, main=paste(paste(pop[i], "Means",sep=' '),nrow(x),sep='\n'), sub=cor(x[,1], x[,2]))
	#points(x[intersect(intersect(names(top.cd4.diff.x), names(top.cd4.diff.y)), names(top.cd4.diff.z)),], col='red')
	points(x[cd4.diff.x > quantile(cd4.diff.x, probs=c(0.99)),], col='red')
	
	plot(y, main=paste(pop[i], "Medians",sep=' '), sub=cor(y[,1], y[,2]))
	#points(y[intersect(intersect(names(top.cd4.diff.x), names(top.cd4.diff.y)), names(top.cd4.diff.z)),], col='red')
	points(y[cd4.diff.y > quantile(cd4.diff.y, probs=c(0.99)),], col='red')
	
	plot(z, main=paste(pop[i], "Variance",sep=' '), sub=cor(z[,1], z[,2]), xlim=c(0,2))
	#points(z[intersect(intersect(names(top.cd4.diff.x), names(top.cd4.diff.y)), names(top.cd4.diff.z)),], col='red')
	points(z[cd4.diff.z > quantile(cd4.diff.z, probs=c(0.99)),], col='red')
	cd4.mean <- x
	cd4.med <- y
	cd4.var <- z	
	
	save(cd4.mean, cd4.med, cd4.var, file=paste('/group/stranger-lab/immvar_data/cd4.',pop[i],'.diff.Robj',sep=''))
}
dev.off()

pdf(file='/group/stranger-lab/czysz/ImmVar/plots/cd14.comparison.pdf')
for (i in seq(1)) {
	cd14.raj.mean <- apply(na.omit(cd14.eset[!duplicated(rownames(cd14.eset)), colnames(cd14.eset)%in%cd14.ids[[i]]]), 1 ,mean)
	cd14.raj.median <- apply(na.omit(cd14.eset[!duplicated(rownames(cd14.eset)), colnames(cd14.eset)%in%cd14.ids[[i]]]), 1 ,median)
	cd14.raj.var <- apply(na.omit(cd14.eset[!duplicated(rownames(cd14.eset)), colnames(cd14.eset)%in%cd14.ids[[i]]]), 1 , var)
	load(file=paste('/group/stranger-lab/immvar_data/', paste('exp_genes.CD14',pop[i],'Robj',sep='.'),sep=''))
	
	cd14.cc.mean <- apply(exp_genes, 1 ,mean)	
	cd14.cc.median <- apply(exp_genes, 1 ,median)	
	cd14.cc.var <- apply(exp_genes, 1, var)
	shared.genes <- intersect(names(cd14.cc.mean), names(cd14.raj.mean))
	
	cc.mean.shared <- cd14.cc.mean[names(cd14.cc.mean)%in%shared.genes]
	raj.mean.shared <- cd14.raj.mean[names(cd14.raj.mean)%in%shared.genes]
	
	x<-c()
	for (gene in names(raj.mean.shared)) { x <- rbind(x, matrix(c(cc.mean.shared[gene], raj.mean.shared[gene]), ncol=2)) }
	rownames(x) <- names(raj.mean.shared)
	colnames(x) <- c("Renormalized", "Downloaded")
	cd14.diff.x <- abs(x[,1] - x[,2])
	top.cd14.diff.x <- cd14.diff.x[cd14.diff.x >= quantile(cd14.diff.x, probs=c(0.99))]
	
	if (F) {
	for (gene in names(tail(sort(top.cd14.diff.x)))) {
		plot(density(exp_genes[as.character(gene),]), col='red', main=as.character(gene),xlim=c(0,15))
		lines(density(cd14.eset[as.character(gene), colnames(cd14.eset)%in%cd14.ids[[i]]]), col='blue')
	}
	#write.table(cd14.diff.x[cd14.diff.x > quantile(cd14.diff.x, probs=c(0.99))], file='/home/t.cri.cczysz/cd14.cau.mean.outlier.txt')
	}
	cc.median.shared <- cd14.cc.median[names(cd14.cc.median)%in%shared.genes]
	raj.median.shared <- cd14.raj.median[names(cd14.raj.median)%in%shared.genes]
	y<-c()
	for (gene in names(raj.median.shared)) { y <- rbind(y, matrix(c(cc.median.shared[gene], raj.median.shared[gene]), ncol=2)) }
	rownames(y) <- names(raj.median.shared)
	colnames(y) <- c("Renormalized", "Downloaded")
	cd14.diff.y <- abs(y[,1] - y[,2])
	top.cd14.diff.y <- cd14.diff.y[cd14.diff.y >= quantile(cd14.diff.y, probs=c(0.99))]

	cc.var.shared <- cd14.cc.var[names(cd14.cc.var)%in%shared.genes]
	raj.var.shared <- cd14.raj.var[names(cd14.raj.var)%in%shared.genes]
	z<-c()
	for (gene in names(raj.var.shared)) { z <- rbind(z, matrix(c(cc.var.shared[gene], raj.var.shared[gene]), ncol=2)) }
	rownames(z) <- names(raj.var.shared)
	colnames(z) <- c("Renormalized", "Downloaded")
	cd14.diff.z <- abs(z[,1] - z[,2])
	top.cd14.diff.z <- cd14.diff.z[cd14.diff.z >= quantile(cd14.diff.z, probs=c(0.99))]

	print(cor(x))
	print(cor(y))
	print(cor(z))
	
	plot(x, main=paste(paste(pop[i], "Means",sep=' '),nrow(x),sep='\n'), sub=cor(x[,1], x[,2]))
	#points(x[intersect(intersect(names(top.cd14.diff.x), names(top.cd14.diff.y)), names(top.cd14.diff.z)),], col='red')
	points(x[cd14.diff.x > quantile(cd14.diff.x, probs=c(0.99)),], col='red')
	
	plot(y, main=paste(pop[i], "Medians",sep=' '), sub=cor(y[,1], y[,2]))
	#points(y[intersect(intersect(names(top.cd14.diff.x), names(top.cd14.diff.y)), names(top.cd14.diff.z)),], col='red')
	points(y[cd14.diff.y > quantile(cd14.diff.y, probs=c(0.99)),], col='red')
	
	plot(z, main=paste(pop[i], "Variance",sep=' '), sub=cor(z[,1], z[,2]), xlim=c(0,2))
	#points(z[intersect(intersect(names(top.cd14.diff.x), names(top.cd14.diff.y)), names(top.cd14.diff.z)),], col='red')
	points(z[cd14.diff.z > quantile(cd14.diff.z, probs=c(0.99)),], col='red')
	
	cd14.mean <- x
	cd14.med <- y
	cd14.var <- z	
	save(cd14.mean, cd14.med, cd14.var, file=paste('/group/stranger-lab/immvar_data/cd14.',pop[i],'.diff.Robj',sep=''))
}
dev.off()

