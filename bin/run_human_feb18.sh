#!/usr/bin/env bash

METADATA=/home/jlazar/proj/mouse_encode/metadata/ENCODE_665-sample_list_WM20180305_ENCODE3_metadata.tsv
OUTPUT_DIR=/home/jlazar/proj/mouse_encode/data_04_18/human_data

mkdir -p $OUTPUT_DIR

./convert_metadata.py $METADATA $METADATA.full.txt &&
METADATA=$METADATA.full.txt

HUMAN_MASTER=/home/meuleman/work/projects/ENCODE3/WM20180312_construct_masterlist_665samples/\
masterlist_DHSs_WM20180313_nonovl_any_indexIDs.txt

cut -f1-3 $HUMAN_MASTER | awk '{print $0"\t"NR-1}' > $OUTPUT_DIR/human_master.bed &&

HUMAN_MASTER=$OUTPUT_DIR/human_master.bed

./map_motifs_genome.sh $HUMAN_MASTER hg38

./global_motif_densities.py $METADATA $OUTPUT_DIR/peak_motifs hg38

### wait until motifs mapped
while [ ! -f "$HUMAN_MASTER.motifs" ] ; do
    sleep 1
done
mdir=`dirname $HUMAN_MASTER`
while [ -f "$mdir/taipale.db" ] ; do
    sleep 1
done


mkdir -p $OUTPUT_DIR/motifs_master
mkdir -p $OUTPUT_DIR/motifs_master_cluster
./global_motif_densities.py $METADATA $OUTPUT_DIR/motifs_master hg38 $HUMAN_MASTER
./global_motif_densities.py $METADATA $OUTPUT_DIR/motifs_master_cluster hg38 $HUMAN_MASTER cluster
####


mkdir -p $OUTPUT_DIR/cis_master
curdir=$(pwd)
filename=$curdir/cis_accessibility_human.py
sbatch $filename $METADATA 5000 $OUTPUT_DIR/cis_master $HUMAN_MASTER


