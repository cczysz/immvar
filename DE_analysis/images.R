library(limma)
library(ggplot2)

data.dir <- "/group/stranger-lab/immvar_data"
plot.dir <<- "/group/stranger-lab/czysz/ImmVar/plots"

cell.types <- c("CD14","CD4")
populations <- c("Caucasian","African-American","Asian")

MakeVolcano <- function(fit) {

}

MakeQQ <- function(fit,cell.type,population) {

	ttable <- topTable(fit,number=Inf)
	file.name <- paste("qq_plot",cell.type,population,"pdf",sep=".")
	pdf(file=paste(plot.dir,file.name,sep="/"))
	qqplot(runif(1000),ttable$P.Value,main=paste("QQ Plot ",cell.type,population,sep="."))
	abline(0,1)
	dev.off()	
	
	file.name <- paste("density",cell.type,population,"pdf",sep=".")
	pdf(file=paste(plot.dir,file.name,sep="/"))
	plot(density(ttable$P.Value),main=paste(cell.type,population))
	dev.off()
}

MakeExprs <- function(fit,residuals) {

}

for (population in populations) {
for (cell.type in cell.types) {
	residual.file <- paste("residuals",cell.type,population,"Robj",sep=".")
	fit.file <- paste("fit",cell.type,population,"Robj",sep=".")
	load(file=paste(data.dir,residual.file,sep="/"))
	load(file=paste(data.dir,fit.file,sep="/"))

	#MakeVolcano(eb.fit)
	MakeQQ(eb.fit,cell.type,population)
	#MakeExprs(eb.fit, expr.residuals)

}
}
