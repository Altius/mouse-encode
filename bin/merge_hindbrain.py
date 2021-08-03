#!/usr/bin/env python
import sys,pandas as pd,subprocess
metadata_file = sys.argv[1]
output_file = sys.argv[2]

replicate_peak_file = '/net/seq/data/aggregations/LN27886/aggregation-2245/peaks/\
LN27886.mm10-encode3-male.uniques.sorted.peaks.starch'

metadata = pd.read_table(metadata_file,index_col = 0,header = 0)
gx = [x for x in metadata.index.values
      if metadata.loc[x,'Tissue']=='midbrain' and metadata.loc[x,'Time']=='11.5']
print gx
if len(gx)!=1: raise IndexError
fl = metadata.loc[gx[0],'Peak-File']

if fl!=output_file:

    with open(output_file,'w') as f1:
        p1 = subprocess.Popen(['bedops','-e','25%',fl,replicate_peak_file],stdout =f1)
        p1.wait()

    metadata.loc[gx[0],'Peak-File'] = output_file
    metadata.to_csv(metadata_file,sep = '\t')





