#!/usr/bin/env python
import pandas as pd, numpy as np,os
from scipy.spatial.distance import pdist,squareform
import scipy.cluster.hierarchy
from sklearn.manifold import MDS
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt


def clustering(count_file,output_dir,subset = None,sep = ','):
    count_data = pd.read_table(count_file,header = 0,index_col = 0,sep = sep)
    good_columns = [c for c in count_data.columns.values if c not in 
            ['DS51473','DS51479','DS51492']]
    count_data = count_data[good_columns]
    if subset is not None:
        count_data = count_data.loc[subset]

    distance_base = os.path.join(output_dir,'cell-distances.txt')
    distances = get_log_distances(count_data)
    distances.to_csv(distance_base,sep = '\t')

    mds_output = os.path.join(output_dir,'cell-mds.txt')
    mds_points = run_mds(distances.values)
    with open(mds_output,'w') as f1:
        lines = []
        for i in range(mds_points.shape[0]):
            lines.append('%s\t%f\t%f' %(distances.columns.values[i],mds_points[i,0],mds_points[i,1]))
        f1.write('\n'.join(lines) + '\n')

    mds_output = os.path.join(output_dir,'cell-pca.txt')
    mds_points = run_mds(distances.values)
    with open(mds_output,'w') as f1:
        lines = []
        for i in range(mds_points.shape[0]):
            lines.append('%s\t%f\t%f' %(distances.columns.values[i],mds_points[i,0],mds_points[i,1]))
        f1.write('\n'.join(lines) + '\n')

    clustering_out = os.path.join(output_dir,'cell-clustering.pdf')
    run_hierarchical(distances,clustering_out)

def run_mds(distances):
    mds_pre = MDS(n_components = 2, dissimilarity = 'precomputed')
    mds_points =  mds_pre.fit_transform(distances)
    return mds_points

def run_pca(distances):
    mds_pre = PCA()
    mds_points =  mds_pre.fit_transform(distances)
    return mds_points

def run_hierarchical(distances,output_file):
    clst = scipy.cluster.hierarchy.ward(distances.values)
    f = plt.figure()
    
    ax = plt.subplot(111)
    c = scipy.cluster.hierarchy.dendrogram(clst,no_labels=True, 
                               link_color_func=lambda k:'k',
                               color_threshold=np.inf)
    ticks = set([])
    for k,icor in enumerate(c['icoord']):
        for i,xlink in enumerate(icor):
            if c['dcoord'][k][i]==0:
                ticks.add(xlink)
    xticks = sorted(list(ticks))
    labels = distances.columns.values
    xt = plt.xticks(xticks, labels, rotation='vertical')
    r = f.axes[0].tick_params(length=0,labelsize='small')
    plt.tight_layout()
    f.savefig(output_file)



def get_log_distances(count_matrix):

    data = np.log(count_matrix.values.astype(float)+1)
    print data.shape
    #data = data.divide(np.std(data[good_columns].values.astype(float),axis = 1),axis = 'index')

    #transformed_distances = squareform(pdist(np.transpose(np.log2(data[good_columns].values[subset] + 1)),
    #                        metric = 'euclidean'))
    transformed_distances = squareform(pdist(np.transpose(data),
                            metric = 'euclidean'))
    print transformed_distances.shape
    distance_frame = pd.DataFrame(transformed_distances,columns = count_matrix.columns.values,
        index = count_matrix.columns.values)
    return distance_frame

