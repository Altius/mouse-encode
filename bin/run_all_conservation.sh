#!/usr/bin/env bash

dirn=`dirname $0`
echo $dirn

#matched_human_networks
#test_dir=/home/sjn/proj/encode3/breeze/networks4jlazar/hg38/results.regulators.all/buffer.5000/LN2346
#test_out=/home/jlazar/proj/mouse_encode/data/conservation/LN2346_3

human_phylop=/home/sjn/proj/encode3/vierstra/footprints/ftd/get.conservation/data/phylop100way.final/phyloP100way.starch

mouse_phylop=/net/fileserv0/vol7/annotations/data/mm10/bed/phyloP60way/phyloP60way.starch
#master_file=/home/jlazar/proj/mouse_encode/data/final_list/master-peaks.mouse.bed4
master_file=$1

#sbatch $dirn/ffl_footprints.py $test_dir $test_out $human_phylop
sbatch $dirn/ffl_footprints.py $master_file $mouse_phylop $master_file.conservation
