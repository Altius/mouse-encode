#!/usr/bin/env python
import pandas as pd,numpy as np,os
import map_dnase
import scipy.stats,sys

def main(metadata_file,master_list,output_dir):
    if not os.path.exists(output_dir): os.mkdir(output_dir)
    metadata = pd.read_table(metadata_file,index_col = 0,header = 0)

    raw_matrix = get_counts(master_list,metadata,output_dir)
    raw_file = os.path.join(output_dir,'raw-counts.txt')
    raw_matrix.to_csv(raw_file,sep = '\t')

    reps = sample_replicates(metadata)

    output = []
    for x in reps:
        if not reps[x]: continue
        corvs = []
        for x2 in reps[x]:
            corvs.append(scipy.stats.pearsonr(np.log(raw_matrix[x].values+1),
                                            np.log(raw_matrix[x2].values + 1))[0])
        output.append('%s\t%f' %(x,np.mean(corvs)))
    output_file = os.path.join(output_dir,'replicate-correlations.txt')
    with open(output_file,'w') as f1: f1.write('\n'.join(output)+'\n')



def sample_replicates(metadata):
    reps = {}
    rep_func = lambda x1,x2: (metadata.loc[x1,'Tissue']==metadata.loc[x2,'Tissue'] and
                              metadata.loc[x1,'Time']==metadata.loc[x2,'Time'])
    for x1 in metadata.index.values:
        reps[x1] = []
        for x2 in metadata.index.values:
            if x1 != x2 and rep_func(x1,x2):
                reps[x1].append(x2)
    return reps

def get_counts(master_list,metadata,output_dir):
    count_matrix = None
    tfiles = {}
    for ds in metadata.index.values:
        tfile = os.path.join(output_dir,'%s.cc.txt' %ds)
        map_dnase.submit2cluster(metadata.loc[ds,'Cut-Count-File'],master_list,tfile)
        tfiles[ds] = tfile
    while any([not os.path.exists(tfiles[r]) for r in tfiles]):
        continue
    for ds in tfiles:
        vls = np.loadtxt(tfiles[ds],usecols=[1])
        if count_matrix is None:
            count_matrix = pd.DataFrame(index = np.arange(len(vls)),columns = metadata.index.values)
        count_matrix[ds] = vls
        os.remove(tfiles[ds])
    return count_matrix

if __name__ == "__main__":
    main(*sys.argv[1:])