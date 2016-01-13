cell.types <- c('CD14', 'CD4')
load('/group/stranger-lab/immvar_rep/mappings/annot.Robj')

meta.files <- list(CD14='/group/stranger-lab/immvar/meta/cd14_meta_all1.txt', CD4='/group/stranger-lab/immvar/meta/cd4_meta_all1.txt')
fit.files <- list(CD14=list(MESA='GSE56045/gencord_fit.Robj', ImmVar='/group/stranger-lab/immvar_data/fit.joint.CD14.Robj', Fairfax='emtab2232/fairfax_fit.Robj'), CD4=list())

for (cell in cell.types) {
	cell <- 'CD14'
	meta.all <- read.table(file=meta.files[[cell]], header=T)
	meta.annots <- annots[meta.all$MarkerName,]


	meta.all$q.value <- p.adjust(meta.all$P.value, method='fdr')
	meta.all.sig <- subset(meta.all, q.value<0.05)
	meta.sig.annot <- annots[meta.all.sig$MarkerName, ]
}
