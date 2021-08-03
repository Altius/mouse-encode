#!/usr/bin/env python
import argparse
import os
import subprocess

import matrix_norm_memmap
import numpy as np,pandas as pd

'''
The main driver script for clustering the mouse master list matrix
will take in the file, log transform

Then for a range of K values will estimate how good it is using MiniBatch kmeans + gap statistic
Final clustering at a given K will be done w/o mini-batching

'''

RANDOM_STATE = 100
SIMULATION_NUMBER = 100
MAX_CPU = 20
SEED = np.random.RandomState(RANDOM_STATE)
BASE_DIR = os.path.split(os.path.realpath(__file__))[0]
MINIBATCH = os.path.join(BASE_DIR,'random_inertia.py')
INERTIA = os.path.join(BASE_DIR,'run_random_inertia.py')
RUN_PCA = os.path.join(BASE_DIR,'pca_matrix.py')
KMEANS = os.path.join(BASE_DIR,'run_kmeans.py')

def main(matrix_file,k_range,parse_matrix = False,output_dir = '',pca_file = None,
         subset = None):
    if not output_dir:
        output_dir = os.path.split(matrix_file)[0]
    if not os.path.exists(output_dir): os.mkdir(output_dir)
    name = os.path.split(matrix_file)[1].rstrip('.txt')

    if parse_matrix:
        tmp_file = os.path.join(output_dir,'%s.npy' %name)
        data = matrix_norm_memmap.main(matrix_file, tmp_file, subset_number=subset)
    else:
        tmp_file = matrix_file
        if not tmp_file.endswith('npy'):
            print 'warning file is not numpy array'

    if pca_file is None:
        pca_file = tmp_file.rstrip('.npy') + '.pca.pickle'
        p1 = subprocess.Popen([RUN_PCA, tmp_file, pca_file])
        p1.communicate()
    elif not os.path.exists(pca_file):
        print 'PCA file does not exist'
        raise
    while not os.path.exists(pca_file): continue

    #random_storage = os.path.join(output_dir,'randomized')
    #if not os.path.exists(random_storage): os.mkdir(random_storage)

    #### right now we are repeating a bunch of work... it would probably be better to s
    ## split out the PCA/data randomization(?) into a separate script
    ### PCA would be easy, randomization might be an bit iffy
    ssh = lambda cmd: ' '.join(['ssh','sched0','sbatch']+cmd)
    for k in k_range:
        true_output = os.path.join(output_dir,'%d-kmeans' %k)
        random_output = os.path.join(true_output,'randomized')
        if not os.path.exists(true_output): os.mkdir(true_output)
        if not os.path.exists(random_output): os.mkdir(random_output)

        #p1 = subprocess.Popen([INERTIA,matrix_file,str(k),pca_file,str(SIMULATION_NUMBER),
        #                  random_output])
        #p1.communicate()
        kmeans_cmd = [KMEANS,tmp_file,str(k),os.path.join(true_output,'%s.kmeans' %k),
                      '--save','--seed','%d' %RANDOM_STATE]
        p1 = subprocess.Popen(ssh(kmeans_cmd),shell = True)
        p1.communicate()

### memory fix: do it by columns: update qnorm,avg,
def parse_matrix_file(matrix_file):
    first_line = open(matrix_file).readline()
    shape = len(first_line.strip().split())
    print shape

    data = np.loadtxt(matrix_file,usecols = np.arange(3,shape))
    
    data = data[SEED.choice(data.shape[0],size = 300000,replace = False)].astype(float)

    return data



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run Kmeans Clustering + Gap Statistic')
    parser.add_argument('-m','--matrix_file',dest = 'matrix_file', type=str,default = None,
                        help='matrix file to kmeans')
    parser.add_argument('-o','--output_directory',dest = 'output_dir', type=str,default = '',
                        help='matrix file to kmeans')
    parser.add_argument('-p','--parse_matrix',dest = 'parse', action='store_true',default = True,
                        help='should you parse the matrix (may include sub-sampling)')
    parser.add_argument('--pca',dest = 'pca',
                        help='pca object')

    ### later might want to make this a bit more robust
    parser.add_argument('-k','--max-k',type = int,default = 30,dest = 'k')
    parser.add_argument('--min-k', type=int, default=1, dest='min_k')
    parser.add_argument('-s', '--subset-number', type=int, default=None, dest='subset_number')

    args = parser.parse_args()
    if args.matrix_file is None:
        matrix_file = os.path.join(args.output_dir,'tmp-matrix-file.txt')
    else:
        matrix_file = args.matrix_file
    print args.subset_number
    krange = np.arange(5,70,5)
    #krange = np.arange(5,105,10)
    main(matrix_file,krange,parse_matrix= args.parse,pca_file = args.pca,subset = args.subset_number,
         output_dir=args.output_dir)



##### don't really use this anymore
def parse_matrix_file_lm(matrix_file):  ## low_memory
    first_line = open(matrix_file).readline()
    shape = len(first_line.strip().split())
    print shape

    values = []
    for i in range(shape - 3):
        print i
        cdata = np.loadtxt(matrix_file, usecols=[i+3])
        if i == 0:
            subset = SEED.choice(len(cdata),size = 200000,replace = False)
        values.append(cdata[subset])
    data = np.column_stack(values)

    return data


def normalize_matrix(data):
    data = qnorm(data)
    data = np.log2(data+1)  ## normalize data by row?  ## other option would be using DESeq to do variance stabilization
    return data

def qnorm(data):
    ref_dist = get_average_distribution(data) ## make things the "typical"  distribution for cleavage
    ref_data = np.sort(ref_dist)
    for i in range(data.shape[1]):
        data[np.argsort(data[:,i]),i] = ref_data
    return data

def get_average_distribution(data):
    ref_dist = np.zeros(data.shape[0])
    for i in range(data.shape[1]):
        ref_dist += data[np.argsort(data[:,i]),i]
    ref_dist /= data.shape[1]
    return ref_dist