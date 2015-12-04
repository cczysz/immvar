#THIS SCRIPT EXECUTES AN ANALYSIS OF EIGHT STUDIES
#THE RESULTS FOR EACH STUDY ARE STORED IN FILES Inputfile1.txt THROUGH Inputfile8.txt

#LOAD THE FIRST EIGHT INPUT FILES

# UNCOMMENT THE NEXT LINE TO ENABLE GenomicControl CORRECTION
# GENOMICCONTROL ON

# === DESCRIBE AND PROCESS THE FIRST INPUT FILE ===
MARKER GENES
ALLELE REF ALT
STRANDLABEL STRAND
EFFECT EFFECT
PVALUE PVALUE 
WEIGHT N

PROCESS /home/t.cri.cczysz/meta_f.Caucasian.CD14.txt 
PROCESS /home/t.cri.cczysz/meta_f.African-American.CD14.txt 
PROCESS /home/t.cri.cczysz/meta_f.Asian.CD14.txt 

OUTFILE immvar_cd14_meta .txt

# Minweight gives minimum total samples needed to output meta results
MINWEIGHT 10
ANALYZE 
CLEAR

MARKER GENES
ALLELE REF ALT
STRANDLABEL STRAND
EFFECT EFFECT
PVALUE PVALUE 
WEIGHT N

PROCESS /home/t.cri.cczysz/meta_f.Caucasian.CD4.txt 
PROCESS /home/t.cri.cczysz/meta_f.African-American.CD4.txt 
PROCESS /home/t.cri.cczysz/meta_f.Asian.CD4.txt 

OUTFILE immvar_cd4_meta .txt

MINWEIGHT 10
ANALYZE 
QUIT
