load('/group/stranger-lab/moliva/ImmVar/data/mRNAexp/EXP_files/stats.Robj')

stats.all[stats.all$ensembl_gene == '-', "ensembl_gene"] <- NA
stats.all[stats.all$symbol_gene == '-', "symbol_gene"] <- NA

ens_comp <- numeric(0)
for (id in rownames(stats.all)) {
	old_ens <- as.character(stats.all[id, "annot2$ens"])
	new_ens <- as.character(stats.all[id, "ensembl_gene"])
	new_ens <- unlist(strsplit(new_ens, split='[.]'))[1]

	# Save comparisons as -1 for different, 0 for same, and binary for NA
	if (!(is.na(old_ens) | is.na(new_ens))) {
		if (old_ens == new_ens) { ens_comp <- c(ens_comp, 0) } 
		else { ens_comp <- c(ens_comp, -1) }
	} else {
		if (is.na(old_ens) & is.na(new_ens)) { ens_comp <- c(ens_comp, 3) } else {
		ens_comp <- c(ens_comp, sum( c(is.na(old_ens), is.na(new_ens)) * 2^c(1,0))) }
	}
}

stats.all <- cbind(stats.all, ens_comp=ens_comp)

sym_comp <- numeric(0)
for (id in rownames(stats.all)) {
	old_sym <- as.character(stats.all[id, "annot2$sym"])
	new_sym <- as.character(stats.all[id, "symbol_gene"])

	# Save comparisons as -1 for different, 0 for same, and binary for NA
	if (!(is.na(old_sym) | is.na(new_sym))) {
		if (old_sym == new_sym) { sym_comp <- c(sym_comp, 0) } 
		else { sym_comp <- c(sym_comp, -1) }
	} else {
		if (is.na(old_sym) & is.na(new_sym)) { sym_comp <- c(sym_comp, 3) } else {
		sym_comp <- c(sym_comp, sum( c(is.na(old_sym), is.na(new_sym)) * 2^c(1,0))) }
	}
}
stats.all <- cbind(stats.all, sym_comp=sym_comp)

save(stats.all, file='/group/stranger-lab/immvar_data/stats.Robj')
