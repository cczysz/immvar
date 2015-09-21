# #!/home/t.cri.cczysz/bin/Rscript
#!/apps/compilers/R/3.1.0/bin/Rscript
### load libraries

load.cel.files <- function(population) {

	load("Robjects/phen.Robj")
	shared_ids <- names(table(phen$ImmVarID2)[table(phen$ImmVarID2)>1])

	cell_type_markers="CD14+16-Mono";
	#immvar_ids <- phen[phen$Race == population, "ImmVarID2"]
	#shared_samples <- unique(immvar_ids[duplicated(immvar_ids)])
	files_sufix=as.character(phen[phen$Race == population & phen$CellType == cell_type_markers & phen$ImmVarID2%in%shared_ids,"FileName"])
	#files_sufix=phen[phen$Race == population & phen$CellType == cell_type_markers & phen$ImmVarID2%in%shared_ids,"FileName"]
	filenames.cd14=unlist(lapply(files_sufix,function(x) paste(c("data/mRNAexp/CEL_files","CD14",x),collapse="/")));

	cell_type_markers="CD4TNve";
	#files_sufix=phen[phen$Race == population & phen$CellType == cell_type_markers & phen$ImmVarID2%in%shared_ids,"FileName"]
	files_sufix=as.character(phen[phen$Race == population & phen$CellType == cell_type_markers & phen$ImmVarID2%in%shared_ids,"FileName"])
	filenames.cd4=unlist(lapply(files_sufix,function(x) paste(c("data/mRNAexp/CEL_files","CD4",x),collapse="/")));

	#raw_oligo=read.celfiles(filenames=c(filenames.cd14,filenames.cd4),phenoData=AnnotatedDataFrame(phen[phen$Race == population ,]))
	raw_oligo=read.celfiles(filenames=c(filenames.cd14,filenames.cd4))
	return(raw_oligo)
}

backcorrect.normalize.probe.level <- function(ExpressionFeatureSet) {
	bgCorrected <- backgroundCorrect(ExpressionFeatureSet);
	normalized2 <- normalize(bgCorrected, method="quantile");
	return(normalized2);
}

library(genefilter)
library(argparse)
library(oligo)
library(pd.hugene.1.0.st.v1)

### Import functions

#sapply(list.files(pattern="[.]R$", path="/group/stranger-lab/moliva/ImmVar/Rfuncs/", full.names=TRUE), source);
#source('/group/stranger-lab/immvar/DE_analysis/generate_exp_object.R')
### MANAGE ARGUMENTS

parser=ArgumentParser()
parser$add_argument("-p", "--population", type="character", default="Caucassian",
    help="Samples' population to filter [default %(default)s], options: Caucassian,African-American,Asian ",
    metavar="character")
#parser$add_argument("-c", "--cell_type", type="character", default="CD14",
    #help="Cell type to filter [default %(default)s], options: CD14, CD4",
    #metavar="character")
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false",
    dest="verbose", help="Print NO output")
args <- parser$parse_args()
args <- list(population="Caucasian", verbose=F)
if( !args$population%in%c("Caucasian","African-American","Asian")) {
	stop("Population needs to be Caucasian African-American or Asian")
}
#if( !args$cell_type%in%c("CD4","CD14")) {
	#stop("Cell type need to be either CD14 or CD4")
#}

###### MAPPING BASED FILTERING
setwd('/group/stranger-lab/moliva/ImmVar')
par.genes<<-as.character(read.table("probes_mapping/annotations/par.genes.txt")[,1])
Y.degenerate<<-as.character(read.table("probes_mapping/annotations/chrY.degenerate.txt")[,1])
Y.degenerate_Xhomolog<<-as.character(read.table("probes_mapping/annotations/chrY.degenerate_Xhomolog.txt")[,1])
Y.all<<-as.character(read.table("probes_mapping/annotations/chrY.genes.txt")[,2])


file="probes_mapping/Robjects/merge_probes_DF_lite.Robj";
if (file.exists(file)) {
        if ( args$verbose ) { print(paste(c("Load file",file),collapse="")) }
        load(file)
} else {
	merge_probes_DF=unique(merge_probes_DF[,c("probeId","gene_ensembl","overlapping_snps","crosshyb_type")])

	## Different cross-hybridization inexes for some probes. To pick a single value per probe, assign the most stringent (highest multimapping index) one:

	dup=as.character(merge_probes_DF$probeId[duplicated(merge_probes_DF$probeId)])
	index_rows=seq(nrow(merge_probes_DF))

	print("Pick one cross-hybridization value per probe")
	for (d in seq(length(dup))) { 
		row=index_rows[merge_probes_DF$probeId%in%dup[d]]; 
		cross=as.character(merge_probes_DF[row,"crosshyb_type"]); 
		if (sum(as.numeric(cross%in%c("1","3"))) == 2) { merge_probes_DF[row[1],"crosshyb_type"] = "3"} 
		else if (sum(as.numeric(cross%in%c("1","0"))) == 2) { merge_probes_DF[row[1],"crosshyb_type"] = "0";} 
		else if (sum(as.numeric(cross%in%c("3","0"))) == 2) { merge_probes_DF[row[1],"crosshyb_type"] = "0";} 
		else {break}; 
		merge_probes_DF=merge_probes_DF[-row[2],]; 
		length(index_rows) <- length(index_rows)-1;
	}
	save(merge_probes_DF,file="probes_mapping/Robjects/merge_probes_DF_lite.Robj")
}

## Counts after pre-filtering

if ( args$verbose ) { print(c("Num genes before filtering:",length(unique(merge_probes_DF$gene_ensembl)))) }
if ( args$verbose ) { print(c("Num probes before filtering:",length(unique(merge_probes_DF$probeId)))) }

print("Filter probes overlapping snps MAF>0.1 in either EUR, EAS, AFR")
probes_snps=unique(merge_probes_DF[merge_probes_DF$overlapping_snps >0,"probeId"])
if ( args$verbose ) { print(c("Num probes overlapping snps with MAF>0.1% in either Cucassian, Asian or African-American populations:",length(probes_snps))) }
# [1] 39846

### Filter multimapping probes

print("Filter multimapping probes")
probes_multimapping=unique(merge_probes_DF[grep(":",merge_probes_DF$gene_ensembl),"probeId"])
if ( args$verbose ) { print(c("Num probes mapping to >1 gene:", length(probes_multimapping))) }
#[1] 73674

### Filtering of probes belonging to probesets with cross-hybridization, except PAR

print("Filter probes belonging to probesets with cross-hybridization, except PAR")
probes_cross=unique(merge_probes_DF[merge_probes_DF$crosshyb_type != 1  & !merge_probes_DF$gene_ensembl%in%par.genes,"probeId"])
if ( args$verbose ) { print(c("Num probes belonging to probesets with cross-hybridization, except PAR:",length(probes_cross))) }
# [1] 47699

filt_probes=unique(c(probes_snps,probes_multimapping,probes_cross))
merge_probes_DF_filt=merge_probes_DF[!merge_probes_DF$probeId%in%filt_probes,]
if ( args$verbose ) { print(c("Num probes filtered:",length(filt_probes))) }
#[1] 122062

## Counts after 1st filtering

if ( args$verbose ) { print(c("Num genes after annotation based filtering:",length(unique(merge_probes_DF_filt$gene_ensembl)))) }
if ( args$verbose ) { print(c("Num probes after annotation based filtering:",length(unique(merge_probes_DF_filt$probeId)))) }

###### EXPRESSION BASED FILTERING
### Background-correction, normalization of probe-level expression values.
### Step already done. Load Robj

if (F) {
file=paste(c("Robjects/oligo.bgSubstractedNormalized",args$cell_type,args$population,"Robj"),collapse=".");
if (file.exists(file)) {
	if ( args$verbose ) { print(paste(c("Load file",file),collapse="")) }
	load(file)
} else {
	if ( args$verbose ) { print(paste(c("Create probe-level, background-corrected and normalized expression values. Store in file:",file),collapse="")) }
	oligo_raw=load.cel.files(args$population,args$cell_type)
	normalized2=backcorrect.normalize.probe.level(oligo_raw)
	save(normalized2,file=file)
}
}

oligo_raw=load.cel.files(args$population)
normalized2=backcorrect.normalize.probe.level(oligo_raw)
# save(normalized2,file=file)
### Eliminate filtered probes

samples=as.character(pData(normalized2)$ImmVarID2)
normalized2=exprs(normalized2)[as.character(rownames(exprs(normalized2))[rownames(exprs(normalized2))%in%merge_probes_DF$probeId]),]
colnames(normalized2)=samples
normalized2=normalized2[as.character(rownames(normalized2)[!rownames(normalized2)%in%filt_probes]),]
normalized2=normalized2[as.character(merge_probes_DF_filt$probeId),]

### Calculate densities gene expression. Identify global minimum, which will be the gene expression threshold.
### Previously, store ImmVarID2 for males and females on vectors of the same name

load("Robjects/phen.Robj")
phen <- phen[phen$Race == args$population,]
cell.cd14 <- c(rep(1,ncol(normalized2)/2),rep(0,ncol(normalized2)/2))
cell.cd4 <- c(rep(0,ncol(normalized2)/2),rep(1,ncol(normalized2)/2))

dCD14 <- density(log2(normalized2[,cell.cd14]))
min_cd14=optimize(approxfun(dCD14$x,dCD14$y),interval=c(3,4))$minimum

dCD4 <- density(log2(normalized2[,cell.cd4]))
min_cd4=optimize(approxfun(dCD4$x,dCD4$y),interval=c(3,4))$minimum

###  Filter Y-linked non-PAR probes for which median expression in females is above threshold and distribution of expressions in male and females are equal (wilxon test pval > 10-5) in any cell type in any population

#file=paste(c("probes_mapping/annotations/",args$cell_type,args$population,"FiltProbeExp.txt"),collapse=".")

file2=paste(c("probes_mapping/annotations/","all","FiltProbeExp.txt"),collapse=".")
filt_Y_probe_exp=as.data.frame(read.table(file2))[,1];

FilterYProbes <- function(exp, sex, merge_probes_DF_filt) {
	filt_Y_probe=vector()
	d.exp <- density(log2(exp[,cell]))
	min_exp=optimize(approxfun(d.exp$x,d.exp$y),interval=c(3,4))$minimum

	median_genexp_males=apply(log2(exp[,!!sex]),1,median)
	median_genexp_females=apply(log2(exp[,!sex]),1,median)

	for (gene in unique(as.character(merge_probes_DF_filt[merge_probes_DF_filt$gene_ensembl%in%Y.all[!Y.all%in%par.genes],"gene_ensembl"]))) {
		probes=as.character(merge_probes_DF_filt[grep(gene,merge_probes_DF_filt$gene_ensembl),"probeId"]);
		for (probe in probes) {
			pval=wilcox.test(log2(exp[probe,!!sex]),log2(exp[probe,!sex]))$p.value
			if (pval > 0.00001 && median_genexp_females_probes[probe] > min_females) {
			filt_Y_probe=c(filt_Y_probe,probe)
			}
		}
	}

	filt_Y_probe=unique(filt_Y_probe)
	return(filt_Y_probe)
}

#filt_Y_probes_CD4 <- FilterYProbes(exp, sex, merge_probes_DF_filt)
#filt_probes=c(filt_Y_probe_exp,filt_probe_min_exp)
filt_probes=as.character(c(filt_Y_probe_exp))

merge_probes_DF_filt=merge_probes_DF_filt[!as.character(merge_probes_DF_filt$probeId)%in%filt_probes,]
normalized2=normalized2[as.character(rownames(normalized2)[!as.character(rownames(normalized2))%in%filt_probes]),]

## Counts after 2nd filtering

if ( args$verbose ) { print(c("Num genes after Y-linked non-PAR genes probes expression based filtering:",length(unique(merge_probes_DF_filt$gene_ensembl)))) }
if ( args$verbose ) { print(c("Num probes after  Y-linked non-PAR genes probes expression based filtering:",length(unique(merge_probes_DF_filt$probeId))))}

exp_genes <- basicRMA(normalized2, as.character(merge_probes_DF_filt$gene_ensembl), normalize=F, background=F)
probes_per_gene=as.data.frame(unlist(table(as.character(merge_probes_DF_filt$gene_ensembl))))

#pdf("Probes_per_gene.pdf")
#hist(probes_per_gene$Freq,breaks=seq(max(probes_per_gene$Freq)))
#abline(v=8)
#dev.off()

### Filter for genes with < 5 probes

genes=as.character(probes_per_gene$Var1[probes_per_gene$Freq > 4])
if ( args$verbose ) { print(c("Num genes after min number of probes (>=5) filtering:",length(genes))) }
if ( args$verbose ) { print(c("Num probes supporting genes with min number of probes (>=7):",length(merge_probes_DF_filt$probeId[merge_probes_DF_filt$gene_ensembl%in%genes]))) }
exp_genes=exp_genes[genes,];


### Filter genes expressed below 10% quantile in >2/3 of males and >2/3 of females
if (F) {
females_filt=kOverA(A=quantile(exp_genes,probs = 0.1),k = round(length(females)/3))
males_filt=kOverA(A=quantile(exp_genes,probs = 0.1),k = round(length(males)/3))
genes_above_threshold_index=genefilter(exp_genes[,males],filterfun(males_filt)) & genefilter(exp_genes[,females],filterfun(females_filt))
#print(dim(exp_genes))
exp_genes=exp_genes[genes_above_threshold_index,]
females_filt=kOverA(A=quantile(exp_genes,probs = 0.1),k = round(length(females))/3)
males_filt=kOverA(A=quantile(exp_genes,probs = 0.1),k = round(length(males))/3)

male_filt_genes <- genefilter(exp_genes[,males],filterfun(males_filt))
female_filt_genes <- genefilter(exp_genes[,females],filterfun(females_filt))

filt_genes_AND <- intersect(male_filt_genes[!male_filt_genes], female_filt_genes[!female_filt_genes])
filt_genes_OR <- union(male_filt_genes[!male_filt_genes], female_filt_genes[!female_filt_genes])
exp_genes=exp_genes[(male_filt_genes | female_filt_genes),]
}
#var.filt.males <- varFilter(exp_genes[,males],var.cutoff=0.25) # Default cutoff of 0.5 removes 50% of genes
#var.filt.females <- varFilter(exp_genes[,females],var.cutoff=0.25)
#var.filt <- intersect(rownames(var.filt.males), rownames(var.filt.females))
#exp_genes <- varFilter(exp_genes, var.cutoff=0.25)
#genes_to_filter <- filt_genes_OR[!(filt_genes_OR%in%var.filt)]
#print(dim(exp_genes))
#exp_genes <- exp_genes[var.filt%in%rownames(exp_genes),]
#exp_genes <- exp_genes[var.filt,]
#print(dim(exp_genes))

#if ( args$verbose ) { print(c("Num genes expressed above 10% quantile in >=1/3 of males and >=1/3 of females:",nrow(exp_genes))) }
#if ( args$verbose ) { print(c("Num of probes supporting genes expressed above 10% quantile in >=1/3 of males and =>=1/3 of females:",length(merge_probes_DF_filt$probeId[merge_probes_DF_filt$gene_ensembl%in%rownames(exp_genes)]))) }

file=paste(c("exp_genes_bt_cell",args$population,"Robj"),collapse=".")
save(exp_genes,file=paste('/group/stranger-lab/immvar_data/',file,sep=''))
