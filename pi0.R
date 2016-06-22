library(qvalue)
library(limma)

load('/group/stranger-lab/immvar_rep/mappings/annot.Robj')
load('/group/stranger-lab/immvar_data/fit.joint.CD14.Robj')
immvar.cd14 <- eb.fit
immvar.cd14$q.value <- p.adjust(immvar.cd14$p.value, method='fdr')
rm(eb.fit)

load('/group/stranger-lab/immvar_data/fit.joint.CD4.Robj')
immvar.cd4 <- eb.fit
immvar.cd4$q.value <- p.adjust(immvar.cd4$p.value, method='fdr')
rm(eb.fit)

rep.names <- c('Fairfax', 'GenCord', 'MESAT', 'MESAM')
rep.fit.objects <- c('emtab2232/fairfax_fit.Robj', 'GenCord/gencord_fit.Robj', 'GSE56580/mesa_tcells_fit.Robj', 'GSE56045/gencord_fit.Robj')
rep.dir <- '/group/stranger-lab/immvar_rep/'

rep.fits <- list()
rep.qvals <- list()
for (i in seq(length(rep.fit.objects))) {
	load(file=paste(rep.dir, rep.fit.objects[i], sep=''))	
	rep.fits[[i]] <- fit
	rep.qvals[[i]] <- qvalue(fit$p.value)
}

immvar.cd14.sig <- rownames(subset(data.frame(immvar.cd14), q.value<0.05))
immvar.cd4.sig <- rownames(subset(data.frame(immvar.cd4), q.value<0.05))

for (i in seq(length(rep.fits))) {
	immvar.cd14.sig.in.rep <- qvalue(rep.fits[[i]][rownames(rep.fits[[i]])%in%immvar.cd14.sig,]$p.value)
	immvar.cd4.sig.in.rep <- qvalue(rep.fits[[i]][rownames(rep.fits[[i]])%in%immvar.cd4.sig,]$p.value)
	print(paste('ImmVar CD14+ vs', rep.names[i]))
	print(immvar.cd14.sig.in.rep$pi0)
	print(paste('ImmVar CD4+ vs', rep.names[i]))
	print(immvar.cd4.sig.in.rep$pi0)
}

pdf(file='/group/stranger-lab/czysz/ImmVar/plots/rep.qvalue.pdf')
lapply(rep.qvals, plot)
dev.off()
