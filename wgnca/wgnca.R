library(WGCNA)
options(stringsAsFactors = FALSE);
enableWGCNAThreads(nThreads=4)

loadExpression <- function(celltype) {
	setwd('/group/stranger-lab/moliva/ImmVar/Robjects')
	#for (cell in c('CD14', 'CD4')) {

	load('phen.Robj')
	if (celltype == 'CD14') {phen <-subset(phen, CellType=="CD14+16-Mono")} else {phen <- subset(phen, CellType=="CD4TNve") }

	load(file=paste(celltype, '.joint.norm.exp_genes.Robj', sep=''))
	if (celltype == 'CD14') { expr <- exp_genes.cd14.joint.norm } else { expr <- exp_genes.cd4.joint.norm }

	expr <- list(male=t(expr[, phen$Sex=='Male']), female=t(expr[, phen$Sex=='Female']))
	return(expr)
}

createTOM <- function(expr) {
	powers = c(c(1:10), seq(from = 12, to=20, by=2))
	sft = pickSoftThreshold(expr, powerVector = powers, verbose = 5)
	softPower = which(-sign(sft$fitIndices[,3])*sft$fitIndices[,2]>0.8)[1];
	adjacency = adjacency(expr, power = softPower);
	dissTOM = 1-TOMsimilarity(adjacency)
	return(dissTOM)
}

identifyModules <- function(geneTree, dissTOM, minModuleSize=30) {

	dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
		deepSplit = 2, pamRespectsDendro = FALSE,
		minClusterSize = minModuleSize);
}

performGOAnalysis <- function(expr, moduleColors) {

	geneNames.ens <- colnames(expr)
	geneNames.sym <- as.character(annots[geneNames.ens, 'symbol_id'])
	ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")
	biomart.results=getBM(ensembl,attributes=c("ensembl_gene_id","hgnc_symbol","entrezgene"),filters="ensembl_gene_id",values=substring(geneNames.ens,1,15))
	uni=biomart.results[!duplicated(biomart.results[,1]),]

	matchings <- (substring(geneNames.ens,1,15)%in%uni[,1])
	geneNames.ens <- geneNames.ens[matchings]
	moduleColors.go <- moduleColors[matchings]


	GOenr <- GOenrichmentAnalysis(moduleColors.go, uni[, 3])
}

setwd('/group/stranger-lab/immvar/wgnca')

for (celltype in c('CD14', 'CD4')) {
	exp.list <- loadExpression(celltype)
for (i in seq(length(exp.list))) {

	pdf.prefix <- paste(celltype, names(exp.list)[i], sep='_')
	print(pdf.prefix)
	dissTOM <- createTOM(exp.list[[i]])
	geneTree = hclust(as.dist(dissTOM), method = "average");

	dynamicMods <- identifyModules(geneTree, dissTOM, 30)
	dynamicColors = labels2colors(dynamicMods)

	pdf(file=paste(pdf.prefix, 'dendogram.pdf', sep='.'), width=8, height=6)
	plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
	dendroLabels = FALSE, hang = 0.03,
	addGuide = TRUE, guideHang = 0.05,
	main = "Gene dendrogram and module colors")
	dev.off()
	
	MEList = moduleEigengenes(expr, colors = dynamicColors)
	MEs = MEList$eigengenes
	# Calculate dissimilarity of module eigengenes
	MEDiss = 1-cor(MEs);
	# Cluster module eigengenes
	METree = hclust(as.dist(MEDiss), method = "average");

	pdf(file=paste(pdf.prefix, 'eigengene_tree.pdf', sep='.'), width=7, height=6)
	plot(METree, main = "Clustering of module eigengenes",
	xlab = "", sub = "")
	MEDissThres = 0.25
	# Plot the cut line into the dendrogram
	abline(h=MEDissThres, col = "red")
	dev.off()

	merge = mergeCloseModules(expr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
	# The merged module colors
	mergedColors = merge$colors;
	# Eigengenes of the new merged modules:
	mergedMEs = merge$newMEs;
	moduleColors = mergedColors
	# Construct numerical labels corresponding to the colors
	colorOrder = c("grey", standardColors(50));
	moduleLabels = match(moduleColors, colorOrder)-1;
	MEs = mergedMEs;

	pdf(file=paste(pdf.prefix, 'dendogram.pdf', sep='.'), width=8, height=6)
	plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
		c("Dynamic Tree Cut", "Merged dynamic"),
		dendroLabels = FALSE, hang = 0.03,
		addGuide = TRUE, guideHang = 0.05)
	dev.off()
	
	plotTOM <- dissTOM^7
	diag(plotTOM) = NA
	pdf(file=paste(pdf.prefix, 'TOM.cormat.pdf', sep='_'), width=9, height=9)
	TOMplot(plotTOM, geneTree, moduleColors)
	dev.off()

	save.image(file=paste(pdf.prefix, 'Robj', sep='.'))
}
}
