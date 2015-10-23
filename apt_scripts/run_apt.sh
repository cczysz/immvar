#PBS -N czysz_apt
#PBS -S /bin/bash
#PBS -l walltime=5:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=15gb

#PBS -o $HOME/apt.out
#PBS -e $HOME/apt.err

APT_DIR=/group/stranger-lab/apt-1.18.0-x86_64-intel-linux/bin
ANNOT_DIR=/group/stranger-lab/moliva/ImmVar/probes_mapping/Hugene_info
#ANNOT_DIR=/group/stranger-lab/nicolel/HuGene-1_0-st-v1-3
HDR=HuGene-1_0-st-v1.r4

if false; then
$APT_DIR/apt-probeset-summarize \
-a rma \
-p $ANNOT_DIR/$HDR.pgf \
-c $ANNOT_DIR/$HDR.clf \
-b $ANNOT_DIR/$HDR.bgp \
-m $ANNOT_DIR/$HDR.mps \
--qc-probesets $ANNOT_DIR/$HDR.qcc \
-o /group/stranger-lab/czysz/ImmVar/all_rma_nokill \
--cel-files /group/stranger-lab/immvar/apt_scripts/cel_file_lists/cd4_cau.txt;
fi

if true; then
$APT_DIR/apt-probeset-summarize \
-a rma \
-p $ANNOT_DIR/$HDR.pgf \
-c $ANNOT_DIR/$HDR.clf \
-b $ANNOT_DIR/$HDR.bgp \
-m $ANNOT_DIR/$HDR.mps \
--qc-probesets $ANNOT_DIR/$HDR.qcc \
--kill-list /group/stranger-lab/immvar/apt_scripts/out.txt \
--cel-files /group/stranger-lab/immvar/apt_scripts/cel_file_lists/cd4_cau.txt \
-o /group/stranger-lab/czysz/ImmVar/all_rma_kill;
fi

if false; then
$APT_DIR/apt-probeset-summarize \
-p $ANNOT_DIR/$HDR.pgf \
-c $ANNOT_DIR/$HDR.clf \
-b $ANNOT_DIR/$HDR.bgp \
--qc-probesets $ANNOT_DIR/$HDR.qcc \
-m $ANNOT_DIR/$HDR.mps \
-o /group/stranger-lab/czysz/ImmVar/all_rma_kill \
--cel-files /home/t.cri.cczysz/celfiles.txt \
--kill-list /group/stranger-lab/immvar/kill_list.txt \
-a rma;
fi
