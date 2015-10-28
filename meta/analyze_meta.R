library(limma)
library(stats)

setwd('/group/stranger-lab/immvar/meta')
cd4.meta <- read.table(file='cd4_meta1.txt', header=T)
cd14.meta <- read.table(file='cd14_meta1.txt', header=T)

bf.cd4 <- 0.05/nrow(cd4.meta)
bf.cd14 <- 0.05/nrow(cd14.meta)

qval.cd4 <- p.adjust(cd4.meta$P.value, method="BH")
cd4.meta <- cbind(cd4.meta, Q.value = qval.cd4)
qval.cd14 <- p.adjust(cd14.meta$P.value, method="BH")
cd14.meta <- cbind(cd14.meta, Q.value = qval.cd14)

meta_sig_cd4_bon <- cd4.meta[cd4.meta$P.value<=bf.cd4, ]
meta_sig_cd14_bon <- cd14.meta[cd14.meta$P.value<=bf.cd14, ]

meta_sig_cd4_bh <- cd4.meta[cd4.meta$Q.value<=0.05, ]
meta_sig_cd14_bh <- cd14.meta[cd14.meta$Q.value<=0.05, ]

setwd('/group/stranger-lab/immvar_data')

cd14.fits<-list()
cd4.fits<-list()
for ( pop in c('Caucasian', 'African-American', 'Asian')) {
for ( cell in c('CD14', 'CD4') ) {

	if (cell=='CD14') { load(paste('fit', pop, cell, 'Robj', sep='.'));cd14.fits[[pop]]<-cd14.eb.fit }
	else { load(paste('fit', pop, cell, 'Robj', sep='.'));cd4.fits[[pop]]<-cd4.eb.fit }
}
}

sig_cd14_bh<-list()
sig_cd14_bon<-list()

sig_cd4_bh<-list()
sig_cd4_bon<-list()

for (i in seq(3)) {
	sig_cd14_bh[[i]]<-topTable(cd14.fits[[i]], number=Inf, p.value=0.05)
	sig_cd14_bon[[i]]<-topTable(cd14.fits[[i]], number=Inf, adjust.method="bonferroni", p.value=0.05)
	
	sig_cd4_bh[[i]]<-topTable(cd4.fits[[i]], number=Inf, p.value=0.05)
	sig_cd4_bon[[i]]<-topTable(cd4.fits[[i]], number=Inf, adjust.method="bonferroni", p.value=0.05)
}
