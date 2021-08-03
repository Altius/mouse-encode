#!/usr/bin/env python
import numpy as np, pandas as pd, subprocess,os,sys,tempfile
import hierarchical_clustering, deseq_norm,map_dnase

def main(metadata_file,master_list,output_dir):
    if not os.path.exists(output_dir): os.mkdir(output_dir)
    master_list = os.path.realpath(master_list)

    output_dir = os.path.realpath(output_dir)
    metadata = pd.read_table(metadata_file,index_col = 0,header = 0)

    raw_matrix = get_counts(master_list,metadata,output_dir)
    raw_file = os.path.join(output_dir,'raw-counts.txt')
    raw_matrix.to_csv(raw_file,sep = '\t')

    peak_matrix = get_peaks(master_list,metadata)
    raw_file = os.path.join(output_dir,'peak-ind.txt')
    peak_matrix.to_csv(raw_file,sep = '\t')

    subset = np.where(np.sum(peak_matrix.values,axis = 1)>0)

    norm_data = normalize(raw_matrix,peak_matrix)
    norm_file = os.path.join(output_dir,'norm-counts.txt')
    norm_data.to_csv(norm_file,sep = '\t')

    ### add in the "brain" here..
    brain_ids = [x for x in metadata.index.values
                    if metadata.loc[x,'Tissue'] in ['midbrain','hindbrain','forebrain']]

    metadata.loc[brain_ids].to_csv(os.path.join(output_dir,'brain.metadata.txt'),sep = '\t')
    raw_matrix[brain_ids].to_csv(os.path.join(output_dir,'brain.raw-counts.txt'),sep = '\t')
    peak_matrix[brain_ids].to_csv(os.path.join(output_dir, 'brain.peak-ind.txt'),sep = '\t')

    #hierarchical_clustering.clustering(norm_file,output_dir,subset = subset,sep = '\t')

def get_called_subset(master_list,peak_files):
    p1 = subprocess.Popen(['bedops','-e','25%',master_list]+list(peak_files),stdout = subprocess.PIPE)
    subset = []
    for line in iter(p1.stdout.readline,''):
        subset.append(int(line.strip().split()[3]))
    p1.stdout.close()
    return np.array(subset)

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
        print ds,len(vls)
        if count_matrix is None:
            count_matrix = pd.DataFrame(index = np.arange(len(vls)),columns = metadata.index.values)
        count_matrix[ds] = vls
        os.remove(tfiles[ds])
    return count_matrix


def qnorm(dt):

    vls = np.sort(dt[dt.columns.values[0]].values)
    norm_data = pd.DataFrame(np.zeros(dt.shape),columns = dt.columns.values,
                             index = dt.index.values)

    for i,cl in enumerate(dt.columns.values):
        order = np.argsort(dt[cl].values)
        zs = np.zeros(len(order))
        zs[order] = vls
        norm_data[cl] = zs
    return norm_data

def normalize(raw_matrix,peak_matrix):
    density_vectors = {c:raw_matrix[c].values for c in raw_matrix}
    called_vectors = {c:peak_matrix[c].values.astype(bool) for c in peak_matrix}

    gm_means = np.zeros(len(density_vectors[density_vectors.keys()[0]]))
    counts = np.zeros(len(gm_means))
    for d in density_vectors:
        gm_means += np.log(density_vectors[d])
        counts[called_vectors[d]] += 1
    gm_means /= float(len(density_vectors))

    norm_vectors = {}
    good = (counts == np.max(counts)) & (gm_means>0)
    for d in density_vectors:
        s = np.median(gm_means[good]-np.log(density_vectors[d][good]))
        norm_vectors[d] = density_vectors[d] * np.exp(s)

    norm_matrix = pd.DataFrame(index = raw_matrix.index.values,columns = raw_matrix.columns.values)
    for c in norm_vectors: norm_matrix[c] = norm_vectors[c]
    return norm_matrix


def get_peaks(master_list,metadata):
    peak_matrix = None
    for ds in metadata.index.values:
        vls = indicator(master_list,metadata.loc[ds,'Peak-File'])
        if peak_matrix is None:
            peak_matrix = pd.DataFrame(index = np.arange(len(vls)),columns = metadata.index.values)
        peak_matrix[ds] = vls
    return peak_matrix

def indicator(master_list,peak_file):
    run = subprocess.Popen(['bedmap', '--indicator','--fraction-either','0.25', master_list, peak_file],
                           stdout=subprocess.PIPE)
    values = []
    for line in iter(run.stdout.readline, ''):
        value = line.strip()
        if value.upper() == 'NAN' or not value:
            values.append(0)
        else:
            values.append(int(float(value)))
    run.stdout.close()
    return np.array(values)

def count(master_list,cc_file):
    run = subprocess.Popen(['bedmap', '--echo', '--sum', master_list, cc_file],
                           stdout=subprocess.PIPE)
    values = []
    for line in iter(run.stdout.readline, ''):
        ref, value = line.strip().split('|')
        if value.upper() == 'NAN' or not value:
            values.append(0)
        else:
            values.append(int(float(value)))
    run.stdout.close()
    return np.array(values)

if __name__ == "__main__":
    print sys.argv
    main(*sys.argv[1:])