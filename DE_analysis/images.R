library(limma)
library(ggplot2)

data.dir <- "/group/stranger-lab/immvar_data"
plot.dir <<- "/group/stranger-lab/czysz/ImmVar/plots"

cell.types <- c("CD14","CD4")
populations <- c("Caucasian","African-American","Asian")

MakeVolcano <- function(fit,cell.type,population) {
	file.n <- paste("volcano",cell.type,population,"pdf",sep=".")

	ttable <- topTable(fit, number=Inf)
	chrs <- c()
	for (gene in rownames(ttable)) {
		dat <- annots[annots$gene_id==gene, ]
		chrs <- c(chrs,dat$chr)
	}	
	g <- ggplot(data=data.frame(ttable), aes(x=logFC,y=-log10(P.Value),color=as.factor(chrs),label=chrs))
	g + geom_text(size=3) + theme(legend.position="none")
	ggsave(filename=paste("/group/stranger-lab/czysz/ImmVar/plots/",file.n,sep=""))
}

MakeQQ <- function(fit,cell.type,population) {

	ttable <- topTable(fit,number=Inf)
	file.name <- paste("qq_plot",cell.type,population,"pdf",sep=".")
	o <- -log10(sort(fit$p.value,decreasing=F))
	e <- -log10(seq(length(o))/length(o))
	N <- length(o)
	c95 <- rep(0,N)
	c05 <- rep(0,N)
	 
	for(i in 1:N){
	c95[i] <- qbeta(0.95,i,N-i+1)
	c05[i] <- qbeta(0.05,i,N-i+1)
	}
	MAX <- max(c(o,e)) 
	pdf(file=paste(plot.dir,file.name,sep="/"))
	qqplot(e,o,ylim=c(0,MAX), main=paste(cell.type,population))
	par(new=T)
	plot(e, -log(c95,10), ylim=c(0,MAX), type="l", 
		axes=FALSE, xlab="", ylab="")
	par(new=T)
	plot(e, -log(c05,10), ylim=c(0,MAX), type="l", 
		axes=FALSE, xlab="", ylab="")
	## add the diagonal
	abline(0,1,col="red")
	if (F) {
	plot(e,o,pch=19,cex=0.25,
		xlab=expression(Expected~~-log[10](italic(p))),
		ylab=expression(Observed~~-log[10](italic(p))),
		xlim=c(0,max(e)),
		ylim=c(0,max(e)))
	}
	dev.off()	
	
	file.name <- paste("density",cell.type,population,"pdf",sep=".")
	pdf(file=paste(plot.dir,file.name,sep="/"))
	plot(density(ttable$P.Value),main=paste(cell.type,population))
	dev.off()
}

MakeExprs <- function(fit,residuals) {

}

annots <<- read.table(file="/group/stranger-lab/moliva/ImmVar/probes_mapping/annotations/gencode.v22.TSS.txt",header=F,as.is=T)
colnames(annots) <- c("gene_id","gene_symbol","chr","start","stop","strand")

for (population in populations) {
for (cell.type in cell.types) {
	residual.file <- paste("residuals",cell.type,population,"Robj",sep=".")
	fit.file <- paste("fit",cell.type,population,"Robj",sep=".")
	load(file=paste(data.dir,residual.file,sep="/"))
	load(file=paste(data.dir,fit.file,sep="/"))

	MakeVolcano(eb.fit,cell.type,population)
	#MakeQQ(eb.fit,cell.type,population)
	#MakeExprs(eb.fit, expr.residuals)

}
}
