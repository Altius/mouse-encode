#!/usr/bin/env bash

set -o pipefail

OUTPUT_DIR=/home/jlazar/proj/mouse_encode/data_05_19
METADATA=/home/jlazar/proj/mouse_encode/metadata/mouse-final-samples.04_18.txt
MASTER_LIST=$OUTPUT_DIR/master-peaks.mouse.bed4
CODE_DIR=/home/jlazar/proj/mouse_encode/bin

mkdir -p $OUTPUT_DIR

erynes_master=/home/erynes/topics/Mouse_2017/Autumn2018/TotalRecall/masterLists/masterlist_DHSs_173samples_ER20180914_nonovl_any.bed1
cut -f1-3 $erynes_master | awk '{print $0"\t"NR-1}' > $MASTER_LIST &&

CUT_DIR=$OUTPUT_DIR/cut_counts

curdir=$(realpath $(dirname $0))

./map_motifs.sh $MASTER_LIST

### should find a way to wait....

### deal with the bad sample
rep_file=$OUTPUT_DIR/midbrain.11.5.replicate_peak_file.bed
./merge_hindbrain.py $METADATA $rep_file &&

### wait until motifs mapped
while [ ! -f "$MASTER_LIST.motifs" ] ; do
    sleep 1
done
mdir=`dirname $MASTER_LIST`
while [ -f "$mdir/taipale.db" ] ; do
    sleep 1
done

### create trans
./global_motif_densities.py $METADATA $OUTPUT_DIR/peak_motifs mm10
./global_motif_densities.py $METADATA $OUTPUT_DIR/master_motifs_cluster mm10 $MASTER_LIST cluster
./global_motif_densities.py $METADATA $OUTPUT_DIR/master_motifs mm10 $MASTER_LIST

### create cis
filename=$curdir/cDHS_accessibility.py
ssh sched0 sbatch $filename $CUT_DIR None $MASTER_LIST $OUTPUT_DIR/cis-signal.5kb.in-DHS.txt 5000


#./cell_data.py $METADATA $MASTER_LIST $OUTPUT_DIR/peak_summary

#./run_all_conservation.sh $MASTER_LIST
#./run_all_conservation.sh $MASTER_LIST.motifs

#phylop_scores=/net/fileserv0/vol7/annotations/data/mm10/bed/phyloP60way/phyloP60way.starch
#awk '$5 < 1e-5 { print; }' $MASTER_LIST.motifs | bgzip > $MASTER_LIST.motifs.1e-5.gz
#zcat $MASTER_LIST.motifs.1e-5.gz | bedmap --faster --echo --skip-unmapped --bp-ovr 1 $phylop_scores - | \
#bgzip > $MASTER_LIST.motifs.phylop_scores.gz
#tabix -p bed $MASTER_LIST.motifs.1e-5.gz
#tabix -p bed $MASTER_LIST.motifs.phylop_scores.gz

#./run_human.sh

mkdir -p $OUTPUT_DIR/peak_summary/brain_anova
brain_metadata=$OUTPUT_DIR/peak_summary/brain.metadata.txt
brain_counts=$OUTPUT_DIR/peak_summary/brain.raw-counts.txt
brain_calls=$OUTPUT_DIR/peak_summary/brain.peak-ind.txt
./tc_anova.py $brain_metadata $brain_counts $brain_calls $OUTPUT_DIR/peak_summary/brain_anova

return 0

### there was some stuff on loops that doesn't look like it is going to make the final cut...