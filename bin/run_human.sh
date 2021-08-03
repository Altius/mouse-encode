#!/usr/bin/env bash

METADATA=/home/jlazar/proj/mouse_encode/data_04_18/mouse-encode.human-samples.full.txt
OUTPUT_DIR=/home/jlazar/proj/mouse_encode/data_04_18/human_data

mkdir -p $OUTPUT_DIR

HUMAN_MASTER=/home/meuleman/work/projects/ENCODE3/WM20180312_construct_masterlist_665samples/\
masterlist_DHSs_WM20180313_nonovl_any_indexIDs.txt




./global_motif_densities.py $METADATA $OUTPUT_DIR/motifs hg38
./global_motif_densities.py $METADATA $HUMAN_MASTER $OUTPUT_DIR/motifs_master hg38
####

#HUMAN_MASTER=/home/erynes/topics/ENCODEpaper2017/Final827/hotspotDirs/hotspotsFDR0.0010/masterList.827samples.FDR0.0010.hg38.bed3

#
# can make the gini here
#
OUTPUT_DIR=/home/jlazar/proj/mouse_encode/human_gini/cis_merge
mkdir -p $OUTPUT_DIR
curdir=$(pwd)
filename=$curdir/cis_accessibility_human.py
#sbatch $filename $METADATA 5000 $OUTPUT_DIR


cells=(Heart Muscle-B Muscle-A Bcell2 Bcell Monocyte Endoderm TH2 TH1 Adipocyte Dendritic NK-Cell Islet TC CD34-erythroid \
CD34-erythroid-hq NPC)
cells=(CD34-erythroid-hq)
for cell in "${cells[@]}"
do
    echo $cell
    metadata=/home/jlazar/proj/mouse_encode/human_gini/$cell.diff-samples.full.txt
    output_dir=/home/jlazar/proj/mouse_encode/human_gini/tc_cis/$cell
    mkdir -p $output_dir
    sbatch $filename $metadata 5000 $output_dir

done



#CD34_meta=/home/jlazar/proj/mouse_encode/human_gini/cd34-diff.full.txt
#CD34_dir=/home/jlazar/proj/mouse_encode/human_gini/cd34_cis

#mkdir -p $CD34_dir

#motif_file=$OUTPUT_DIR/merged-peaks.bed.motifs
#./cis_accessibility_human.py $CD34_meta 5000 $CD34_dir
#CD34_dir=/home/jlazar/proj/mouse_encode/human_gini/cd34_trans
#./global_motif_density_count.py $CD34_meta $CD34_dir hg38 $motif_file

#muscle_meta=/home/jlazar/proj/mouse_encode/human_gini/muscle-diff.full.txt
#muscle_dir=/home/jlazar/proj/mouse_encode/human_gini/muscle_cis

#mkdir -p $muscle_dir

#./cis_accessibility_human.py $muscle_meta 5000 $muscle_dir


#fibro_meta=/home/jlazar/proj/mouse_encode/human_gini/fibroblast-diff.full.txt
#fibro_dir=/home/jlazar/proj/mouse_encode/human_gini/fibroblast_cis

#mkdir -p $fibro_dir

#./cis_accessibility_human.py $fibro_meta 5000 $fibro_dir