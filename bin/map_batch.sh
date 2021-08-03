#!/bin/bash
#

#SBATCH --mem=1G                 # memory pool for all cores
#SBATCH -e slurm.%N.%j.err        # STDERR

motif_file=$1
bed_file=$2
output=$3

bedops -e -1 $motif_file $bed_file > $output