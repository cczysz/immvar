cell.types <- c('CD14', 'CD4')
load('/group/stranger-lab/immvar_rep/mappings/annot.Robj')

meta.files <- list(CD14='/group/stranger-lab/immvar/meta/cd14_meta_all1.txt', CD4='/group/stranger-lab/immvar/meta/cd4_meta_all1.txt')
fit.files <- list(
	CD14=list(MESA='/group/stranger-lab/immvar_rep/GSE56045/gencord_fit.Robj', ImmVar='/group/stranger-lab/immvar_data/fit.joint.CD14.Robj', Fairfax='/group/stranger-lab/immvar_rep/emtab2232/fairfax_fit.Robj'),
	CD4=list(MESA='/group/stranger-lab/immvar_rep/GSE56580/mesa_tcells_fit.Robj', ImmVar='/group/stranger-lab/immvar_data/fit.joint.CD4.Robj', Gencord='/group/stranger-lab/immvar_rep/GenCord/gencord_fit.Robj'))

fits <- list()

for (cell in cell.types) {
	cell <- 'CD4'
	meta.all <- read.table(file=meta.files[[cell]], header=T)
	meta.annots <- annots[meta.all$MarkerName,]

	meta.all$q.value <- p.adjust(meta.all$P.value, method='fdr')
	meta.all.sig <- subset(meta.all, q.value<0.05)
	meta.sig.annot <- annots[meta.all.sig$MarkerName, ]

	fitnames <- names(fit.files[[cell]])
	for (i in seq(3)) {
		load(file=fit.files[[cell]][[i]])
		if (!(fitnames[i] == 'ImmVar')) {fit$q.value <- p.adjust(fit$p.value, method='fdr'); fits[[cell]][[fitnames[i]]] <- data.frame(fit)
		} else {eb.fit$q.value <- p.adjust(eb.fit$p.value, method='fdr');fits[[cell]][[fitnames[i]]] <- data.frame(eb.fit)}
	}

	sig.weight.table <- table(meta.all.sig$Weight)

	weight.data <- list()
	for (weight in names(sig.weight.table)) {
		sig.genes <- subset(meta.all.sig, Weight==weight)
		sig.annot <- annots[sig.genes$MarkerName, ]
		
		dat <- NULL
		for (gene in sig.genes$MarkerName) {
		x <- data.frame(row.names=gene)
		for (i in seq(3)) {
			fit.name <- names(fits[[cell]])[i]
			fit <- fits[[cell]][[i]] 
			if (gene%in%rownames(fit)){
			x[i] <- fit[gene, 'q.value']
			} else { x[i] <- NA }
		}
		dat <- rbind(dat, x)
		}
	weight.data[[weight]] <- dat
	}
}
