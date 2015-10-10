importExp <- function(cell) {
	exp.list <- list()
	for (pop in c('Caucasian','African-American','Asian')) {
		file.n=paste('exp_genes',cell,pop,'Robj',sep='.')
		load(file=paste('/group/stranger-lab/moliva/ImmVar/Robjects/',file.n,sep=''))
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
		x <- c(x,calculateVst(na.omit(exp1[gene, ]), na.omit(exp2[gene, ])))
	}
	names(x) <- shared.genes
	return(x)
}

library(oligo)
library(GEOquery)

setwd('/group/stranger-lab/immvar_data/')

exp.cd14 <- importExp("CD14")
exp.cd4 <- importExp("CD4")

load('/group/stranger-lab/moliva/ImmVar/Robjects/phen.Robj')
load('/group/stranger-lab/moliva/ImmVar/probes_mapping/Robjects/merge_probes_DF.Robj')

phen.cd14 <<- phen[phen$PhenotypeMarkers=="CD3- CD14+ CD16-", ]
phen.cd4 <<- phen[phen$PhenotypeMarkers=="CD3+ CD14- CD8- CD4+ CD62L+ CD25-", ]

pops <- c("Caucasian", "African-American", "Asian")
#if (!file.exists('/group/stranger-lab/czysz/ImmVar/vst.sex.cd14.Robj')) {
cd14.vst.list <- list()
for (i in seq(3)) {
	tmp.phen <- phen.cd14[phen.cd14$Race==pops[i], ]
	vst <- pairwiseVst(exp.cd14[[i]][, (tmp.phen$Sex=="Male")], exp.cd14[[i]][, (tmp.phen$Sex=="Female")])
	cd14.vst.list <- unlist(list(cd14.vst.list, list(vst)), recursive=FALSE)
	} 
save(cd14.vst.list,file='/group/stranger-lab/czysz/ImmVar/vst.sex.cd14.Robj')
#} else load('/group/stranger-lab/czysz/ImmVar/vst.cd14.Robj')

cd14.shared.genes <- intersect(intersect(names(cd14.vst.list[[1]]), names(cd14.vst.list[[2]])), names(cd14.vst.list[[3]]))

for (i in seq(3)) { cd14.vst.list[[i]] <- cd14.vst.list[[i]][cd14.shared.genes] }

cd14.top.diff.genes <- c()
for ( i in seq(3)) {
	top.diff <- cd14.vst.list[[i]][cd14.vst.list[[i]] > quantile(cd14.vst.list[[i]], probs=c(0.99),na.rm=T)]
	cd14.top.diff.genes <- c(cd14.top.diff.genes, names(top.diff))	
}

#if (!file.exists('/group/stranger-lab/czysz/ImmVar/vst.sex.cd4.Robj')) {
cd4.vst.list <- list()
for (i in seq(3)) {
	tmp.phen <- phen.cd4[phen.cd4$Race==pops[i], ]
	vst <- pairwiseVst(exp.cd4[[i]][, (tmp.phen$Sex=="Male")], exp.cd4[[i]][, (tmp.phen$Sex=="Female")])
	cd4.vst.list <- unlist(list(cd4.vst.list, list(vst)), recursive=FALSE)
	} 
	save(cd4.vst.list,file='/group/stranger-lab/czysz/ImmVar/vst.sex.cd4.Robj')
#} else load('/group/stranger-lab/czysz/ImmVar/vst.cd4.Robj')

cd4.shared.genes <- intersect(intersect(names(cd4.vst.list[[1]]), names(cd4.vst.list[[2]])), names(cd4.vst.list[[3]]))

for (i in seq(3)) { cd4.vst.list[[i]] <- cd4.vst.list[[i]][cd4.shared.genes] }

cd4.top.diff.genes <- c()
for ( i in seq(3)) {
	top.diff <- cd4.vst.list[[i]][cd4.vst.list[[i]] > quantile(cd4.vst.list[[i]], probs=c(0.99),na.rm=T)]
	cd4.top.diff.genes <- c(cd4.top.diff.genes, names(top.diff))	
}

annot <- read.table(file='/group/stranger-lab/moliva/ImmVar/probes_mapping/annotations/gencode.v22.TSS.txt',header=T,row.names=1,as.is=T)

cd14.rep <- c("UTS2","FLNB","SPTBN1","SMAGP","EMP1","TTC39C","PSPH","SPATA20","PPIL3","F2RL1","LRRC6","RFX2","LMNA","P2RX5","PRH1","SIGLEC14","GATM","LILRA3")
#rep_geneid <- annot[annot$symbol_id%in%cd14.rep,]
cd14.rep_geneid <- unique(merge_probes_DF[merge_probes_DF$gene_symbols%in%cd14.rep,"gene_ensembl"])

cd4.rep <- c("UTS2","CRIP2","NR1D1","C11orf21","VIM","TRPM2","TMEM14C","RPL36AL","PTCH1","PSPH","HOXB2","HEBP2","GPR137B","AFAP1","CCDC144A","FHIT","PPFIBP2","GSTM4")
cd4.rep_geneid <- unique(merge_probes_DF[merge_probes_DF$gene_symbols%in%cd4.rep,"gene_ensembl"])

cd14.all <- cbind(cd14.vst.list[[1]][names(cd14.vst.list[[1]])%in%cd14.rep_geneid],
	cd14.vst.list[[2]][names(cd14.vst.list[[2]])%in%cd14.rep_geneid],
	cd14.vst.list[[3]][names(cd14.vst.list[[3]])%in%cd14.rep_geneid])
colnames(cd14.all) <- c("EU-AA","EU-EA","AA-EA")

cd4.all <- cbind(cd4.vst.list[[1]][names(cd4.vst.list[[1]])%in%cd4.rep_geneid],
	cd4.vst.list[[2]][names(cd4.vst.list[[2]])%in%cd4.rep_geneid],
	cd4.vst.list[[3]][names(cd4.vst.list[[3]])%in%cd4.rep_geneid])
colnames(cd4.all) <- c("EU-AA","EU-EA","AA-EA")
