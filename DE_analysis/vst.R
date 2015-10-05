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

library(oligo)
library(GEOquery)

setwd('/group/stranger-lab/immvar_data/')

#exp.cd14 <- importExp("CD14")
#exp.cd4 <- importExp("CD4")

load('/group/stranger-lab/immvar_rep/GSE56034.Robj')
eset <- exprs(gse56034[[1]])
eset <- log2(eset)
pdat <- pData(gse56034[[1]]) # Population not given. Must overlap with phen

EU_PC <- as.matrix(read.table(file='/group/stranger-lab/immvar_data/GSE56034_GSM.ImmVarCD14.EU.PC20.txt', header=T, row.names=1))
EA_PC <- as.matrix(read.table(file='/group/stranger-lab/immvar_data/GSE56034_GSM.ImmVarCD14.EA.PC10.txt', header=T, row.names=1))
AA_PC <- as.matrix(read.table(file='/group/stranger-lab/immvar_data/GSE56034_GSM.ImmVarCD14.AA.PC14.txt', header=T, row.names=1))

cau_ids <- colnames(EU_PC)
asn_ids <- colnames(EA_PC)
afr_ids <- colnames(AA_PC)

load('/group/stranger-lab/moliva/ImmVar/Robjects/phen.Robj')
load('/group/stranger-lab/moliva/ImmVar/probes_mapping/Robjects/merge_probes_DF.Robj')

phen.cd14 <<- phen[phen$PhenotypeMarkers=="CD3- CD14+ CD16-",]

geo_ids <- as.character(pdat$title)
split_ids <- strsplit(geo_ids,split="[.]")
geo_ids <- c()
for (i in seq(length(split_ids))) {
	id <- split_ids[[i]][1]
	geo_ids <- c(geo_ids,id)
}

convertIDs <- function(id) { 
	race <- as.character(phen.cd14[phen.cd14$ImmVarID2==id,]$Race)
	return(race)
}
race <- apply(as.matrix(geo_ids),1,convertIDs)
race <- as.character(race)

pops <- list(cau_ids, afr_ids, asn_ids)
if (F) {
vst.list <- list()
for (i in seq(2)) {
for (j in (i+1):3) {
	#vst <- pairwiseVst(exp.cd14[[i]], exp.cd14[[j]])
	vst <- pairwiseVst(eset[, colnames(eset)%in%pops[[i]]],  eset[, colnames(eset)%in%pops[[j]]])
	vst.list <- unlist(list(vst.list, list(vst)), recursive=FALSE)
}
}}

save(vst.list,file='/group/stranger-lab/czysz/ImmVar/pc.vst.towfique.Robj')

shared.genes <- intersect(intersect(names(vst.list[[1]]), names(vst.list[[2]])), names(vst.list[[3]]))

for (i in seq(3)) { vst.list[[i]] <- vst.list[[i]][shared.genes] }

top.diff.genes <- c()
for ( i in seq(3)) {
	top.diff <- vst.list[[i]][vst.list[[i]] > quantile(vst.list[[i]], probs=c(0.99),na.rm=T)]
	top.diff.genes <- c(top.diff.genes, names(top.diff))	
}

annot <- read.table(file='/group/stranger-lab/moliva/ImmVar/probes_mapping/annotations/gencode.v22.TSS.txt',header=T,row.names=1,as.is=T)
cd14.rep <- c("UTS2","FLNB","SPTBN1","SMAGP","EMP1","TTC39C","PSPH","SPATA20","PPIL3","F2RL1","LRRC6","RFX2","LMNA","P2RX5","PRH1","SIGLEC14","GATM","LILRA3")
#rep_geneid <- annot[annot$symbol_id%in%cd14.rep,]
rep_geneid <- merge_probes_DF[merge_probes_DF$gene_symbols%in%cd14.rep,"transcript_cluster_id"]
rep_geneid <- unique(rep_geneid)

cd4.rep <- c("UTS2","CRIP2","NR1D1","C11orf21","VIM","TRPM2","TMEM14C","RPL36AL","PTCH1","HOXB2","HEBP2","GPR137B","AFAP1","CCDC144A","FHIT","PPFIBP2","GSTM4")

all <- cbind(vst.list[[1]][names(vst.list[[1]])%in%rep_geneid],vst.list[[2]][names(vst.list[[2]])%in%rep_geneid],vst.list[[3]][names(vst.list[[3]])%in%rep_geneid])
#all <- cbind(vst.list[[1]][names(vst.list[[1]])%in%(rep_geneid)],vst.list[[2]][names(vst.list[[2]])%in%rownames(rep_geneid)],vst.list[[3]][names(vst.list[[3]])%in%rownames(rep_geneid)])
colnames(all) <- c("EU-AA","EU-EA","AA-EA")

write.table(cbind(rep_geneid$symbol_id,all),file='/home/t.cri.cczysz/vst.rep.txt',quote=F)
