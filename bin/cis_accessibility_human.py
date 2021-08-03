#!/usr/bin/env python
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G                 # memory pool for all cores
#SBATCH -e slurm.%N.%j.err        # STDERR

from collections import defaultdict
import pandas as pd,subprocess,os,numpy as np,sys
import multiprocessing as mp,tempfile


TSS_FILE = '/home/sjn/proj/encode3/dhs-paper/annotate-master-list/data/gencode.v26-basic-tss.bed'

def main(metadata_file,distance,output_dir,master_file = None):
    metadata = pd.read_table(metadata_file,header = 0)
    cut_counts_files = {}
    peak_files = {}
    for ind in metadata.index.values:
        file = metadata.loc[ind,'Cut-Count-File']
        if pd.isnull(file) or not file.endswith('.starch'): continue
        name = metadata.loc[ind,'DS-ID']
        cut_counts_files[name] = file
        peak_files[name] = metadata.loc[ind,'Peak-File']
        if not os.path.exists(cut_counts_files[name]): print cut_counts_files[name], 'shouldnt happen'


    if master_file is None:
        master_file = os.path.join(output_dir,'merged.bed')
        merge_peak_files([peak_files[p] for p in peak_files],master_file)
    tss_file = os.path.join(output_dir,'tss-peaks.bed')
    create_combined_cdhs(master_file, tss_file, distance)

    name_order = cut_counts_files.keys()

    obase = os.path.join(output_dir, '%s.tss-dhs-data.txt')
    pool = mp.Pool(processes = max(32,mp.cpu_count()-8))
    results = [pool.apply_async(get_gene_values,args=(cut_counts_files[name],tss_file,obase %name))
                for name in name_order]
    pool.close()
    pool.join()
    output = [s.get() for s in results]

    output_file = os.path.join(output_dir,'summary-data.txt')
    out_data = None
    for i,nm in enumerate(cut_counts_files.keys()):

        gene_counts = output[i]
        numbers = np.sum([gene_counts[g] for g in gene_counts])

        shape = (len(gene_counts),len(cut_counts_files))
        if out_data is None: out_data = pd.DataFrame(np.zeros(shape),index = gene_counts.keys(),columns=
                                        cut_counts_files.keys())

        for g in gene_counts: out_data.loc[g,nm] = gene_counts[g] / float(numbers) * 1000000
    out_data.to_csv(output_file,sep = '\t')


def parse_number_file(file):
    return int(float(open(file).readline().split()[3]))

def merge_peak_files(peak_files,output_file,max_number = 500):
    if len(peak_files) <= max_number:
        with open(output_file,'w') as f1:
            p1 = subprocess.Popen(['bedops','-m']+[x for x in peak_files],stdout = f1)
            p1.communicate()
    else:
        tmp_files = []
        for i in range(0,len(peak_files),max_number):
            x,tmpf = tempfile.mkstemp()
            merge_peak_files(peak_files[i:i+max_number],tmpf)
            tmp_files.append(tmpf)
        merge_peak_files(tmp_files,output_file)
        for tmpf in tmp_files:
            os.remove(tmpf)

def create_combined_cdhs(master_list_file,output_file,distance):
    with open(output_file,'w') as f1:
        p1 = subprocess.Popen(['bedmap','--range','%d' %distance,'--echo','--echo-map-id-uniq',master_list_file,TSS_FILE],
                              stdout=subprocess.PIPE)
        p2 = subprocess.Popen(['sort-bed','-'],stdin = subprocess.PIPE,stdout = subprocess.PIPE)
        p3 = subprocess.Popen(['uniq'],stdin = p2.stdout,stdout = f1)
        for line in iter(p1.stdout.readline,''):
            ln,genes = line.strip().split('|')
            cdata = ln.split()[0:3]
            for r in genes.split(';'):
                if not r: continue
                p2.stdin.write('\t'.join(cdata+['%s' %r,'.'])+'\n')
        f2 = open(output_file)
        for line in f2:
            dt = line.strip().split()[0:4] + ['.']
            p2.stdin.write('\t'.join(dt)+'\n')
        f2.close()
        p1.stdout.close()
        p2.stdin.close()
        p2.stdout.close()
        p3.wait()


def get_gene_values(cut_count_file,cdhs_file,ensg = True):
    gene_counts = defaultdict(int)

    run = subprocess.Popen(['bedmap','--echo','--sum',cdhs_file,cut_count_file],
                           stdout = subprocess.PIPE)
    for line in iter(run.stdout.readline,''):
        ref,value = line.strip().split('|')
        if value.upper() == 'NAN' or not value: continue
        name = ref.split('\t')[3]
        gene_counts[name.upper()] += int(float(value))

    run.stdout.close()
    return gene_counts

if __name__ == "__main__":
    main(sys.argv[1],int(sys.argv[2]),sys.argv[3])