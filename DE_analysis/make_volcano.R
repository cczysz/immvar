library(ggplot2)
library(limma)

makeVolcanoPlot <- function(fit, annot, file.out) {

	if (!("symbol"%in%names(fit))) {
		fit$symbol <- annot[rownames(fit), "symbol_id"]	
		fit$chr <- annot[rownames(fit), "chr"]
	}

	fit <- data.frame(fit)
	fit$p.value <- -log10(fit$p.value)
	pdf(file=file.out)
	g <- ggplot(data=fit, aes(x=coefficients, y=p.value)) + geom_point() + 
		geom_point(data=fit[fit$chr=='chrY',],col='red') + 
		geom_point(data=fit[fit$chr=='chrX',], col='blue')
	return(g)
}
