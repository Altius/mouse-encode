#!/bin/bash
#
#SBATCH -N 1                      # number of nodes
#SBATCH -n 1                      # number of cores
#SBATCH --mem=1G                 # memory pool for all cores
#SBATCH -e slurm.%N.%j.err        # STDERR

bedops -u "$@"
rm "$@"


