#!/bin/bash

file=$1
genome=$2
#genome=hg19
output=$file.motifs
outdir=`dirname $file`
cdir=`dirname $0`

motif_dir=/net/seq/data/projects/motifs/fimo

module load bedops

taipale=$motif_dir/$genome.taipale.1e-4/fimo.combined.1e-4.parsed.starch
uniprobe=$motif_dir/$genome.uniprobe.1e-4/fimo.combined.1e-4.parsed.starch
xfac=$motif_dir/$genome.xfac.1e-4/fimo.combined.1e-4.parsed.starch
jaspar=$motif_dir/$genome.jaspar.1e-4/fimo.combined.1e-4.parsed.starch


jid1=$(sbatch $cdir/map_batch.sh $taipale $file $outdir/taipale.db)
jid1=${jid1##* }
jid2=$(sbatch $cdir/map_batch.sh $uniprobe $file $outdir/uniprobe.db)
jid2=${jid2##* }
jid3=$(sbatch $cdir/map_batch.sh $xfac $file $outdir/xfac.db)
jid3=${jid3##* }
jid4=$(sbatch $cdir/map_batch.sh $jaspar $file $outdir/jaspar.db)
jid4=${jid4##* }

t=$outdir/taipale.db
u=$outdir/uniprobe.db
x=$outdir/xfac.db
j=$outdir/jaspar.db

echo $jid1
echo $jid2


jid5=$(sbatch -o $output --dependency=afterany:$jid1:$jid2:$jid3:$jid4  $cdir/merge_db.sh $t $u $x $j)






