library(limma)

wd='/group/stranger-lab/immvar_data/'
setwd(wd)

for (cell in c('CD14', 'CD4')) {

f_in <- paste('fit.joint', cell, 'Robj', sep='.')
load(file=f_in)

if (cell == 'CD14') { 
	cd14.fit <- eb.fit
	ttable<-topTable(eb.fit, number=Inf)
	sig<-topTable(eb.fit, number=Inf, p.value=0.05)
} else {
	cd4.fit <- eb.fit
	ttable<-topTable(eb.fit, number=Inf) 
	sig<-topTable(eb.fit, number=Inf, p.value=0.05)
}
}

files.in <- c('emtab2232/fairfax_fit.Robj', 'GenCord/gencord_fit.Robj', 'GSE56580/mesa_tcells_fit.Robj', 'GSE56045/gencord_fit.Robj')
studies <- c('Fairfax', 'Gencord', 'MesaT', 'MesaM')
sample_size <- as.numeric(c(432, 85, 227, 1264))
# fit,fit
fit_dir <- "/group/stranger-lab/immvar_rep/"

rep.fits <- list()
for (i in seq(length(files.in))) {
	load(paste(fit_dir, files.in[i], sep=''))
	rep.fits[[i]] <- fit
}

pdf(file='/group/stranger-lab/czysz/ImmVar/plots/fc_plots.pdf', width=8, height=8)
for (i in seq(length(files.in)-1)) {
for (j in seq(i+1, length(files.in))) {
	x <- data.frame(gene=rownames(rep.fits[[i]]$coefficient), logFC=rep.fits[[i]]$coefficient)
	y <- data.frame(gene=rownames(rep.fits[[j]]$coefficient), logFC=rep.fits[[j]]$coefficient)
	joint <- merge(x, y, by="gene")
	title <- paste(studies[i], studies[j], 'logFC', sep=' ')
	smoothScatter(joint[, -1], main=title, xlab=studies[i], ylab=studies[j], sub=cor(joint[, 2], joint[, 3]))
	abline(0,1)

	# Compare rep to immvar results	
	y <- data.frame(gene=rownames(cd14.fit$coefficient), logFC=cd14.fit$coefficient)
	joint <- merge(x, y, by="gene")
	title <- paste(studies[i], 'ImmVar CD14', 'logFC', sep=' ')
	smoothScatter(joint[, -1], main=title, xlab=studies[i], ylab='ImmVar CD14', sub=cor(joint[, 2], joint[, 3]))
	abline(0,1)

	y <- data.frame(gene=rownames(cd4.fit$coefficient), logFC=cd4.fit$coefficient)
	joint <- merge(x, y, by="gene")
	title <- paste(studies[i], 'ImmVar CD4', 'logFC', sep=' ')
	smoothScatter(joint[, -1], main=title, xlab=studies[i], ylab='ImmVar CD4', sub=cor(joint[, 2], joint[, 3]))
	abline(0,1)
}
}
dev.off()
