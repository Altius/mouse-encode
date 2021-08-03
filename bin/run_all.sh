#!/usr/bin/bash



MASTER_LIST=/home/jlazar/proj/mouse_encode/data/final_list/master-peaks.mouse.bed4
CDHS_FILE=/home/jlazar/proj/mouse_encode/data/final_list/corrs_distalsFirst_above0.7_celltypes_500kb.bed8
CDHS_FILE=/home/jlazar/proj/mouse_encode/data/final_list/cDHSs.from-bytestore.bed

MATRIX_FILE=/home/jlazar/proj/mouse_encode/data/final_list/mouse-peaks.mouse.matrix.txt
METADATA=/home/jlazar/proj/mouse_encode/data/mouse-final-samples.full.txt
#METADATA_STRICT=/home/jlazar/proj/mouse_encode/data/mouse-strict-samples_resequence.txt
CUT_DIR=/home/jlazar/proj/mouse_encode/data/cut_counts
CODE_DIR=/home/jlazar/proj/mouse_encode/bin

OUTPUT_DIR=/home/jlazar/proj/mouse_encode/data/nov_17

mkdir -p $OUTPUT_DIR

curdir=$(realpath $(dirname $0))

#./map_motifs.sh $MASTER_LIST

#./global_motif_densities.py $METADATA $OUTPUT_DIR/peak_motifs mm10
./global_motif_densities.py $METADATA $OUTPUT_DIR/master_motifs_cluster mm10 $MASTER_LIST cluster
#./global_motif_density_count.py $METADATA $MASTER_LIST $OUTPUT_DIR/score_motif mm10

#brain_metadata=/home/jlazar/proj/mouse_encode/data/mouse-brain-samples.full.txt
#heart_metadata=/home/jlazar/proj/mouse_encode/data/mouse-heart-samples.full.txt

#./cell_distances.py $brain_metadata $MASTER_LIST $OUTPUT_DIR/brain_mds
#./cell_distances.py $heart_metadata $MASTER_LIST $OUTPUT_DIR/heart_mds
#./cell_distances.py $developmental_data $MASTER_LIST $OUTPUT_DIR/develop_mds

#developmental_data=/home/jlazar/proj/mouse_encode/data/mouse-final-samples-extra.full.txt
#./library_test.py $developmental_data $MASTER_LIST $OUTPUT_DIR/batch_library
#library_metadata=$OUTPUT_DIR/batch_library/library-metadata.txt
#./global_motif_densities.py $library_metadata $OUTPUT_DIR/batch_library/motifs mm10 $MASTER_LIST




####
echo $curdir
filename=$curdir/cDHS_accessibility.py
echo $filename
#ssh sched0 sbatch $filename $CUT_DIR $CDHS_FILE $MASTER_LIST $OUTPUT_DIR/cis-signal.cDHSs.plus-1kb.txt 1000
#ssh sched0 sbatch $filename $CUT_DIR $CDHS_FILE $MASTER_LIST $OUTPUT_DIR/cis-signal.cDHSs.txt
#ssh sched0 sbatch $filename $CUT_DIR $CDHS_FILE $MASTER_LIST $OUTPUT_DIR/cis-signal.cDHSs.plus-5kb.txt 5000
#ssh sched0 sbatch $filename $CUT_DIR None $MASTER_LIST $OUTPUT_DIR/cis-signal.5kb.in-DHS.txt 5000

#sbatch $filename $CUT_DIR $CDHS_FILE $OUTPUT_DIR/cis-signal.cDHSs.25kb.txt 25000
#sbatch $filename $CUT_DIR $CDHS_FILE $OUTPUT_DIR/cis-signal.cDHSs.50kb.txt 50000
#sbatch $filename $CUT_DIR $CDHS_FILE $OUTPUT_DIR/cis-signal.cDHSs.75kb.txt 75000
#sbatch $filename $CUT_DIR $CDHS_FILE $OUTPUT_DIR/cis-signal.cDHSs.100kb.txt 100000
#sbatch $filename $CUT_DIR $CDHS_FILE $OUTPUT_DIR/cis-signal.cDHSs.150kb.txt 150000
#sbatch $filename $CUT_DIR $CDHS_FILE $OUTPUT_DIR/cis-signal.cDHSs.250kb.txt 250000
#ssh sched0 'sbatch $curdir/cDHS_accessibility.py $CUT_DIR $CDHS_FILE $OUTPUT_DIR/cis-signal.cDHSs.75kb.txt 75000'
#ssh sched0 'sbatch $curdir/cDHS_accessibility.py $CUT_DIR $CDHS_FILE $OUTPUT_DIR/cis-signal.cDHSs.100kb.txt 100000'
#ssh sched0 'sbatch $curdir/cDHS_accessibility.py $CUT_DIR $CDHS_FILE $OUTPUT_DIR/cis-signal.cDHSs.150kb.txt 150000'
#ssh sched0 'sbatch $curdir/cDHS_accessibility.py $CUT_DIR $CDHS_FILE $OUTPUT_DIR/cis-signal.cDHSs.250kb.txt 250000'

filename=$curdir/cis_accessibility.py
#./cis_accessibility.py $CUT_DIR 1000 $OUTPUT_DIR/cis-signal.1kb.txt
#./cis_accessibility.py $CUT_DIR 2500 $OUTPUT_DIR/cis-signal.2.5kb.txt
#ssh sched0 sbatch $filename $CUT_DIR 5000 $OUTPUT_DIR/cis-signal.5kb.txt
#./cis_accessibility.py $CUT_DIR 10000 $OUTPUT_DIR/cis-signal.10kb.txt
#./cis_accessibility.py $CUT_DIR 10000 $OUTPUT_DIR/cis-signal.10kb.txt
#./cis_accessibility.py $CUT_DIR 25000 $OUTPUT_DIR/cis-signal.25kb.txt

#filename=$curdir/tf_info.py
#$filename $MATRIX_FILE $MASTER_LIST $OUTPUT_DIR/all-tissues --ds_file $METADATA
#./motif_density.py $MATRIX_FILE $MASTER_LIST $OUTPUT_DIR/all-tissues --ds_file $METADATA
