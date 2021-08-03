#!/bin/bash

script_dir=`dirname "$(readlink -f "$0")"`
echo $script_dir

file=$1
output=$file.motifs
outdir=`dirname $file`


motif_dir=/net/seq/data/projects/motifs/fimo

module load bedops

taipale=$motif_dir/hg19.taipale.1e-4/fimo.combined.1e-4.parsed.starch
uniprobe=$motif_dir/hg19.uniprobe.1e-4/fimo.combined.1e-4.parsed.starch
xfac=$motif_dir/hg19.xfac.1e-4/fimo.combined.1e-4.parsed.starch
jaspar=$motif_dir/hg19.jaspar.1e-4/fimo.combined.1e-4.parsed.starch

echo $outdir
echo $jaspar
echo $uniprobe
echo $xfac
echo $taipale
echo $file

jid1=$(sbatch $script_dir/map_batch.sh $taipale $file $outdir/taipale.db)
jid1=${jid1##* }
jid2=$(sbatch $script_dir/map_batch.sh $uniprobe $file $outdir/uniprobe.db)
jid2=${jid2##* }
jid3=$(sbatch $script_dir/map_batch.sh $xfac $file $outdir/xfac.db)
jid3=${jid3##* }
jid4=$(sbatch $script_dir/map_batch.sh $jaspar $file $outdir/jaspar.db)
jid4=${jid4##* }

t=$outdir/taipale.db
u=$outdir/uniprobe.db
x=$outdir/xfac.db
j=$outdir/jaspar.db

echo $jid1
echo $jid2


jid5=$(sbatch -o $output --dependency=afterany:$jid1:$jid2:$jid3:$jid4  $script_dir/merge_db.sh $t $u $x $j)






