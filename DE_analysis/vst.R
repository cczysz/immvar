importExp <- function(cell) {
	exp.list <- list()
	for (pop in c('Caucasian','African-American','Asian')) {
		load(file=paste('exp_genes',cell,pop,'Robj',sep='.'))
		exp.list <- unlist(list(exp.list, list(exp_genes)), recursive=FALSE)
	}
	return(exp.list)
}

calculateVst <- function(pop1.exp, pop2.exp) {

	n1 <- length(pop1.exp)
	n2 <- length(pop2.exp)
	v1 <- var(pop1.exp)
	v2 <- var(pop2.exp)
	Vs <- (v1*n1 + v2*n2)/(n1+n2)

	Vt <- var(c(pop1.exp, pop2.exp))
	Vst <- (Vt - Vs)/Vt
	return(Vst)
}

pairwiseVst <- function(exp1, exp2) {
	x <- numeric()
	shared.genes <- intersect(rownames(exp1), rownames(exp2))
	for (gene in shared.genes) {
		x <- c(x,calculateVst(exp1[gene, ], exp2[gene, ]))
	}
	names(x) <- shared.genes
	return(x)
}

setwd('/group/stranger-lab/immvar_data/')

exp.cd14 <- importExp("CD14")
exp.cd4 <- importExp("CD4")

vst.list <- list()
for (i in seq(2)) {
for (j in (i+1):3) {
	vst <- pairwiseVst(exp.cd14[[i]], exp.cd14[[j]])
	vst.list <- unlist(list(vst.list, list(vst)), recursive=FALSE)
}
}

shared.genes <- intersect(intersect(names(vst.list[[1]]), names(vst.list[[2]])), names(vst.list[[3]]))

for (i in seq(3)) { vst.list[[i]] <- vst.list[[i]][shared.genes] }

vst.list[[1]] <- vst.list[[1]]^2
vst.list[[2]] <- vst.list[[2]]^2
vst.list[[3]] <- vst.list[[3]]^2

top.diff.genes <- c()
for ( i in seq(3)) {
	top.diff <- vst.list[[i]][vst.list[[i]] > quantile(vst.list[[i]], probs=c(0.99))]
	top.diff.genes <- c(top.diff.genes, names(top.diff))	
}

#top.diff.genes <- names(top.diff.genes[table(top.diff.genes) == 3])
