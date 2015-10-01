load.cel.files <- function(population,cell_type) {

	library(oligo)
	load("Robjects/phen.Robj")

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
	shared_ids <- names(table(phen$ImmVarID2)[table(phen$ImmVarID2)>1])
	#files_sufix=as.character(phen[phen$Race == population & phen$CellType == cell_type_markers & phen$ImmVarID2%in%shared_ids,"FileName"])
	files_sufix=as.character(phen[phen$Race == population & phen$CellType == cell_type_markers,"FileName"])
	filenames=unlist(lapply(files_sufix,function(x) paste(c("data/mRNAexp/CEL_files",cell_type,x),collapse="/")));
	raw_oligo=read.celfiles(filenames=filenames,phenoData=AnnotatedDataFrame(phen[phen$Race == population & phen$CellType == cell_type_markers,]))
	#raw_oligo=read.celfiles(filenames=filenames,phenoData=AnnotatedDataFrame(phen[phen$Race == population & phen$CellType == cell_type_markers & phen$ImmVarID2%in%shared_ids,]))
	return(raw_oligo)
}

backcorrect.normalize.probe.level <- function(ExpressionFeatureSet) {

	library(oligo)
	bgCorrected <- backgroundCorrect(ExpressionFeatureSet);
	normalized2 <- normalize(bgCorrected, method="quantile");
	return(normalized2);
}
