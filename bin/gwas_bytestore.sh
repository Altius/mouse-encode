#!/bin/sh

outdir=../data/jan_18/gwas
MASTER_LIST=/home/jlazar/proj/mouse_encode/data/final_list/master-peaks.mouse.bed4
GWAS_HITS=/home/jlazar/proj/mouse_encode/data/jan_18/gwas/10-6.gwas.snps.bed
LIFT_FILE=/home/jlazar/liftOver/hg19ToMm10.over.chain.gz

mkdir -p $outdir

MOUSE_HITS=$outdir/10-6.gwas.snps.mm10.bed
#/home/rsandstrom/liftOver/liftOver $GWAS_HITS $LIFT_FILE $MOUSE_HITS $MOUSE_HITS.unmapped.bed


GWAS_MASTER=$outdir/master_list.gwas.bed
sort-bed $MOUSE_HITS | bedmap --skip-unmapped --echo --min --delim '\t' --sci $MASTER_LIST - > $GWAS_MASTER 

GWAS_CDHS=$outdir/mouse.gwas-hits.cdhss.bed


/net/seq/data/projects/bytestore/github/byte-store/share/query-bytestore.sh --mm10-198sample-dnaseI-pearsonr-092017 --within-range 0.7:1.0 $GWAS_MASTER > $GWAS_CDHS

bedmap --echo --count --exact --delim '\t' $GWAS_CDHS | cut -f 1-3,8 | uniq > $GWAS_CDHS.count


MATCH=$outdir/master_list.gwas.control.bed
MATCH_CDHS=$outdir/mouse.gwas-hits.control-cdhss.bed


/net/seq/data/projects/bytestore/github/byte-store/share/query-bytestore.sh --mm10-198sample-dnaseI-pearsonr-092017 --within-range 0.7:1.0 $MATCH > $MATCH_CDHS

bedmap --echo --count --exact --delim '\t' $MATCH_CDHS | cut -f 1-3,8 | uniq > $MATCH_CDHS.count



#rm extra.bed


#/net/seq/data/projects/bytestore/github/byte-store/share/query-bytestore.sh --mm10-198sample-dnaseI-pearsonr-092017 --within-range 0.7:1.0 $MASTER_LIST > $outdir/mouse.cdhss.bed
