load.cel.files <- function(population) {

	library(oligo)
	load("Robjects/phen.Robj")

if (F) {
	if ( cell_type == "CD14" ) {
		
		## Load EUR CD14 CEL files
		cell_type_markers="CD14+16-Mono";
	} 
	else if ( cell_type == "CD4" ) {
		## Load EUR CD4 CEL files
		cell_type_markers="CD4TNve";
	} 
	else {
		stop("Cell type needs to be CD14 or CD4")
	}
	if ( !population%in%c("Caucasian","African-American","Asian") ) {
		stop("Population needs to be Caucasian African-American or Asian")
	}
}
	cell_type_markers="CD14+16-Mono";
	files_sufix=as.character(phen[phen$Race == population & phen$CellType == cell_type_markers,"FileName"])
	filenames.cd14=unlist(lapply(files_sufix,function(x) paste(c("data/mRNAexp/CEL_files","CD14",x),collapse="/")));

	cell_type_markers="CD4TNve";
	files_sufix=as.character(phen[phen$Race == population & phen$CellType == cell_type_markers,"FileName"])
	filenames.cd4=unlist(lapply(files_sufix,function(x) paste(c("data/mRNAexp/CEL_files","CD4",x),collapse="/")));

	#raw_oligo=read.celfiles(filenames=c(filenames.cd14,filenames.cd4),phenoData=AnnotatedDataFrame(phen[phen$Race == population ,]))
	raw_oligo=read.celfiles(filenames=c(filenames.cd14,filenames.cd4))
	return(raw_oligo)
}

backcorrect.normalize.probe.level <- function(ExpressionFeatureSet) {

	library(oligo)
	bgCorrected <- backgroundCorrect(ExpressionFeatureSet);
	normalized2 <- normalize(bgCorrected, method="quantile");
	return(normalized2);
}
