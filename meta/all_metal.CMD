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

PROCESS /home/t.cri.cczysz/meta_f.Gencord.txt
PROCESS /home/t.cri.cczysz/meta_f.MesaT.txt
PROCESS /home/t.cri.cczysz/meta_f.joint.CD4.txt

OUTFILE cd4_meta_weight .txt

# Minweight gives minimum total samples needed to output meta results
MINWEIGHT 719
ANALYZE 
CLEAR

MARKER GENES
ALLELE REF ALT
STRANDLABEL STRAND
EFFECT EFFECT
PVALUE PVALUE 
WEIGHT N

PROCESS /home/t.cri.cczysz/meta_f.MesaM.txt
PROCESS /home/t.cri.cczysz/meta_f.Fairfax.txt
PROCESS /home/t.cri.cczysz/meta_f.joint.CD14.txt

OUTFILE cd14_meta_weight .txt

MINWEIGHT 2097
ANALYZE 
QUIT
