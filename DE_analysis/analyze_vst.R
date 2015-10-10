load('/group/stranger-lab/czysz/ImmVar/vst.cd14.Robj')
cd14.vst.cc <- cd14.vst.list;rm(cd14.vst.list)
load('/group/stranger-lab/czysz/ImmVar/vst.cd4.Robj')
cd4.vst.cc <- cd4.vst.list;rm(cd4.vst.list)
load('/group/stranger-lab/czysz/ImmVar/cd14.vst.towfique.Robj')
cd14.vst.raj <- cd14.vst.list;rm(cd14.vst.list)
load('/group/stranger-lab/czysz/ImmVar/cd4.vst.towfique.Robj')
cd4.vst.raj <- cd4.vst.list;rm(cd4.vst.list)

load('/group/stranger-lab/moliva/ImmVar/probes_mapping/Robjects/merge_probes_DF.Robj')

cd14.shared.ensembl <- intersect(intersect(names(cd14.vst.cc[[1]]), names(cd14.vst.cc[[2]])), names(cd14.vst.cc[[3]]))
cd4.shared.ensembl <- intersect(intersect(names(cd4.vst.cc[[1]]), names(cd4.vst.cc[[2]])), names(cd4.vst.cc[[3]]))

cd14.shared.probes <- intersect(intersect(names(cd14.vst.raj[[1]]), names(cd14.vst.raj[[2]])), names(cd14.vst.raj[[3]]))
cd4.shared.probes <- intersect(intersect(names(cd4.vst.raj[[1]]), names(cd4.vst.raj[[2]])), names(cd4.vst.raj[[3]]))

convertEnsemblID <- function(id, merge_probes_DF) {
	e_id <- as.character(merge_probes_DF[merge_probes_DF$gene_ensembl==id, "transcript_cluster_id"][1])	
	e_id <- unlist(strsplit(e_id, split=":"))[1]
	return(e_id)
}

convertprobeID <- function(id, merge_probes_DF) {
	e_id <- as.character(merge_probes_DF[merge_probes_DF$transcript_cluster_id==id, "gene_ensembl"][1])	
	e_id <- unlist(strsplit(e_id, split=":"))[1]
	return(e_id)
}
#cd14.raj.ensembl <- apply(as.matrix(names(cd14.vst.raj[[1]])),1,convertprobeID, merge_probes_DF=merge_probes_DF)
#cd4.raj.ensembl <- apply(as.matrix(names(cd4.vst.raj[[1]])),1,convertprobeID, merge_probes_DF=merge_probes_DF)

#save(cd14.raj.ensembl,file='/group/stranger-lab/czysz/ImmVar/cd14.probe2ens.Robj')
#save(cd4.raj.ensembl,file='/group/stranger-lab/czysz/ImmVar/cd4.probe2ens.Robj')

load('/group/stranger-lab/czysz/ImmVar/cd4.probe2ens.Robj')
load('/group/stranger-lab/czysz/ImmVar/cd14.probe2ens.Robj')

tests <- c('EUvsAA', 'EUvsEA', 'AAvsEA')
for (i in seq(3)) { 
	names(cd4.vst.raj[[i]]) <- cd4.raj.ensembl 
	shared.genes <- intersect(names(cd4.vst.raj[[i]]), names(cd4.vst.cc[[i]]))

	x <- NULL
	for (gene in shared.genes) {
		x <- rbind(x, matrix(c(cd4.vst.raj[[i]][gene], cd4.vst.cc[[i]][gene]),ncol=2))
	}
	rownames(x) <- shared.genes
	raj.sig <- quantile(x[,1],probs=c(0.99)); cc.sig <- quantile(x[,2],probs=c(0.99))
	file.n <- paste('cd4','cor',tests[i],'pdf',sep='.')
	pdf(file=paste('/group/stranger-lab/czysz/ImmVar/plots/',file.n,sep=''))
	plot(x, main=paste(tests[i], cor(x[,1],x[,2]), sep=' '), xlab='Raj', ylab='CC')
	abline(0,1,col='black'); abline(v=raj.sig,col='red'); abline(h=cc.sig, col='blue')
	dev.off()
}

for (i in seq(3)) { 
	names(cd14.vst.raj[[i]]) <- cd14.raj.ensembl 
	shared.genes <- intersect(names(cd14.vst.raj[[i]]), names(cd14.vst.cc[[i]]))

	x <- NULL
	for (gene in shared.genes) {
		x <- rbind(x, matrix(c(cd14.vst.raj[[i]][gene], cd14.vst.cc[[i]][gene]),ncol=2))
	}
	rownames(x) <- shared.genes
	raj.sig <- quantile(x[,1],probs=c(0.99)); cc.sig <- quantile(x[,2],probs=c(0.99))
	file.n <- paste('cd14','cor',tests[i],'pdf',sep='.')
	pdf(file=paste('/group/stranger-lab/czysz/ImmVar/plots/',file.n,sep=''))
	plot(x, main=paste(tests[i],cor(x[,1],x[,2]),sep=' '), xlab='Raj', ylab='CC')
	abline(0,1,col='black'); abline(v=raj.sig,col='red'); abline(h=cc.sig, col='blue')
	dev.off()
}

if (F) {
cd14.raj.ensembl <- c()
for (id in cd14.shared.probes) {
	e_id <- as.character(merge_probes_DF[merge_probes_DF$transcript_cluster_id==id, "gene_ensembl"][1])	
	e_id <- unlist(strsplit(e_id, split=":"))[1]
	cd14.raj.ensembl <- c(cd14.raj.ensembl, e_id)
}
names(cd14.vst.raj) <- cd14.raj.ensembl

cd4.raj.ensembl <- c()
for (id in cd4.shared.probes) {
	e_id <- as.character(merge_probes_DF[merge_probes_DF$transcript_cluster_id==id, "gene_ensembl"][1])	
	e_id <- unlist(strsplit(e_id, split=":"))[1]
	cd4.raj.ensembl <- c(cd4.raj.ensembl, e_id)
}
names(cd4.vst.raj) <- cd4.raj.ensembl
}
