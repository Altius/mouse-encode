#!/usr/bin/env python
import subprocess
import pandas as pd,numpy as np,sys,os
from collections import defaultdict
import motifs
import scipy.stats

TSS_FILE = '/home/jlazar/annotations/gencode.v19.all_tss.ensembl_id.bed'
GENE_NAME_FILE = '/home/jlazar/annotations/ensg2gene_names.txt'
DISTANCE = 12500

def main(input_dir,background_file,output_base,use_tfs = True):
    c2g = motifs.cluster2gene_mapping(motifs.CLUSTER_MAPPINGS_FILE)
    file_dict = parse_dir(input_dir)

    all_tfs = set([])
    for m in c2g: all_tfs.update(c2g[m])
    all_tfs = list(all_tfs)

    shape = (len(all_tfs),len(file_dict.keys()))


    enr_data = pd.DataFrame(np.zeros(shape), index=all_tfs,
                            columns=file_dict.keys())
    pvalue_data = pd.DataFrame(np.zeros(shape), index=all_tfs,
                            columns=file_dict.keys())

    background_counts = get_peak_count(background_file)
    background_size = get_line_count(background_file)
    for i,name in enumerate(file_dict.keys()):
        print name
        peak_counts = get_peak_count(file_dict[name])
        sample_size = get_line_count(file_dict[name])
        for tf in all_tfs:
            pseudo = (background_counts[tf]+1)/float(background_size+1)
            enr_data.loc[tf,name] = np.log2((((peak_counts[tf]+1)/float(sample_size+1))+pseudo)/
                                            (((background_counts[tf]) / float(background_size + 1))+pseudo))
            pv = 1-scipy.stats.hypergeom.cdf(peak_counts[tf]-.1,background_size,
                                               background_counts[tf],sample_size)
            pvalue_data.loc[tf,name] = -1*np.log10(pv)

    enr_data.to_csv(output_base + '.enr.txt',sep = '\t')
    pvalue_data.to_csv(output_base + '.log10-pvalue.txt', sep='\t')


def parse_dir(dir):
    data = {}
    for file in os.listdir(dir):
        name = file.split('forFIMO')[0].rstrip('_')
        data[name] = os.path.join(dir,file)
    return data


def get_peak_count(bed_file,ensg = False):
    gene_counts = defaultdict(int)
    gene_names = ''
    if not ensg:
        ensg2gene = pd.read_table(GENE_NAME_FILE,header = None,index_col = 0)
    run = subprocess.Popen(['bedmap','--range',str(DISTANCE),'--echo','--count',TSS_FILE,bed_file],
                           stdout = subprocess.PIPE)
    for line in iter(run.stdout.readline,''):
        ref,count = line.strip().split('|')
        name = ref.split('\t')[3]
        if ensg:
            gene_counts[name] += int(count)
        else:
            if name not in ensg2gene.index:
                print name,line
            else:
                gene_counts[ensg2gene.loc[name,1]] += int(count)

    run.stdout.close()
    return gene_counts

def get_line_count(file):
    count = 0
    with open(file) as f1:
        for line in f1: count += 1
    return count


if __name__ == "__main__":
    input_dir = sys.argv[1]
    bed_file = sys.argv[2]
    output = sys.argv[3]
    main(input_dir,bed_file,output)