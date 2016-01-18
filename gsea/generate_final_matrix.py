import os

files = os.listdir()

filetypes = ["matrix.txt", "filter.nompval.0.05.matrix.txt","filter.fdr.0.25.matrix.txt"]
gene_set_types = ["c1.all","c2.cgp","c2.cp.biocarta","c2.cp.kegg","c2.cp.reactome","c3.mir","c3.tft","c4.cgn","c4.cm","c5.bp","c5.cc","c5.mf","c6.all","c7.all", "hg.all"]

for gene_set in gene_set_types:
	for ft in filetypes:
		with open('%s.top20_gene_sets.na_pos.%s' % (gene_set, ft)) as f:
			lines = f.readlines()	
			header = line[0]
			for line in lines[1:]:
				values=line.split.strip()
