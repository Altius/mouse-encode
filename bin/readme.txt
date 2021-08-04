Code for generating data matrices for mouse encode manuscript

run_all.sh - main driver script to contain data matrices
run_all_conservation.sh - script called by run_all.sh to run conservation analyses
run_all_human.sh - script called by run_all.sh to run human-mouse comparisons
map_motifs.sh - script to map TF motifs to DHS sequences
map_motifs_genome.sh - script to map TF motifs to DHS sequences for arbitrary genome

global_motif_densities.py - script to identify TF motif densities in DHS peaks for each mouse DNaseI sample
cDHS_accessibility.py - script to identify DNaseI cleavage density in DHSs around gene TSS for each mouse DNaseI sample
cell_data.py - script to identify DNaseI cleavage density int DHSs for each mouse DNaseI sample
merge_hindbrain.py - script to identify replicate peaks in technical replicates of hindbrain samples
tc_anova.py - script to identify variable DHSs in mouse fetal brain samples

kmeans - scripts for running kmeans clusterin on DNaseI index data
lims - scripts for extracting DnaseI experiment metadata from the Altius Institute's Lab Information Management System API
