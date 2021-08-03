#!/usr/bin/env python
#
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G                 # memory pool for all cores
#SBATCH -e slurm.%N.%j.err        # STDERR


from collections import defaultdict
import pandas as pd,subprocess,os,numpy as np,sys
import multiprocessing as mp


TSS_FILE = '/home/jlazar/annotations/gencode.vM9.mouse.all_tss.bed'
GENE_NAME_FILE = '/home/jlazar/annotations/ensg2gene_names.mouse.txt'
DISTANCE = 10000

def main(cut_directory,cdhs_file,master_list_file,output_file,distance = None):
    cut_counts_files = {}
    number_files = {}
    out_data = None
    for file in os.listdir(cut_directory):
        if not file.endswith('.starch'): continue
        try:
            name = file.split('.')[1]
        except IndexError:
            print file, 'huh?'
            continue
        cut_counts_files[name] = os.path.join(cut_directory,file)
        number_files[name] = os.path.join(cut_directory,'%s-counts.txt' %name)
        if not os.path.exists(number_files[name]): print number_files[name]

    print len(cut_counts_files)

    if distance is not None:
        distance = int(distance)
        nm = os.path.split(cdhs_file)[1]
        out_dir = os.path.split(output_file)[0]
        #stripped_cdhs_file = os.path.join(out_dir,nm.rstrip('.bed8')+'within-%dkb.bed' %(distance/1000))
        #create_stripped_cdhs(cdhs_file,stripped_cdhs_file,distance)
        additional_cdhs_file = os.path.join(out_dir,nm.rstrip('.bed8')+'.plus-%dkb.bed' %(distance/1000))
        create_combined_cdhs(cdhs_file,master_list_file,additional_cdhs_file,distance)
        cdhs_file = additional_cdhs_file

    name_order = cut_counts_files.keys()
    pool = mp.Pool(processes = max(32,mp.cpu_count()-8))
    results = [pool.apply_async(get_gene_values,args=(cut_counts_files[name],cdhs_file))
                for name in name_order]
    pool.close()
    pool.join()
    output =  [s.get() for s in results]


    for i,nm in enumerate(cut_counts_files.keys()):

        gene_counts = output[i]
        #numbers = np.sum([gene_counts[g] for g in gene_counts])
        numbers = int(open(number_files[nm]).readline().split()[3])

        shape = (len(gene_counts),len(number_files))
        if out_data is None: out_data = pd.DataFrame(np.zeros(shape),index = gene_counts.keys(),columns=
                                        number_files.keys())

        for g in gene_counts: out_data.loc[g,nm] = gene_counts[g] / float(numbers) * 1000000
    out_data.to_csv(output_file,sep = '\t')

def create_combined_cdhs(cdhs_file,master_list_file,output_file,distance):
    print cdhs_file
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
        if cdhs_file != 'None':
            f2 = open(cdhs_file)
            for line in f2:
                dt = line.strip().split()[0:4] + ['.']
                p2.stdin.write('\t'.join(dt)+'\n')
            f2.close()
        p1.stdout.close()
        p2.stdin.close()
        p2.stdout.close()
        p3.wait()

def create_stripped_cdhs(cdhs_file,output_file,distance):
    p1 = subprocess.Popen(['bedmap','--range','%d' %distance,'--echo','--echo-map-id-uniq',cdhs_file,TSS_FILE],
                          stdout=subprocess.PIPE)
    f1 = open(output_file,'w')
    for line in iter(p1.stdout.readline,''):
        ln,genes = line.strip().split('|')
        nm = line.split()[3]
        if any([nm==r for r in genes.split(';')]):
            f1.write(ln+'\n')
    f1.close()

def combine_tss(cdhs_file,output_file,distance):
    with open(output_file,'w') as f1:
        p1 = subprocess.Popen(['bedmap','--range','%d' %distance,'--echo','--echo-map-id-uniq',cdhs_file,TSS_FILE],
                              stdout=subprocess.PIPE)
        p1 = subprocess.Popen(['bedops','-u','-',cdhs_file,TSS_FILE],
                              stdout=subprocess.PIPE)
        for line in iter(p1.stdout.readline,''):
            ln,genes = line.strip().split('|')
            nm = line.split()[3]
            if not any([nm==r for r in genes.split(';')]):
                pass

    pass

def parse_number_file(file):
    return int(float(open(file).readline().split()[3]))

def get_gene_values(cut_count_file,cdhs_file,ensg = True):
    gene_counts = defaultdict(int)
    gene_names = ''
    if not ensg:
        ensg2gene = pd.read_table(GENE_NAME_FILE,header = None,index_col = 0)

    run = subprocess.Popen(['bedmap','--echo','--sum',cdhs_file,cut_count_file],
                           stdout = subprocess.PIPE)
    for line in iter(run.stdout.readline,''):
        ref,value = line.strip().split('|')
        if value.upper() == 'NAN' or not value: continue
        name = ref.split('\t')[3]
        gene_counts[name.upper()] += int(float(value))

        #if ensg:
        #    prev = gene_counts[name.upper()]
        #    gene_counts[name.upper()] = max(prev,int(float(value)))
        #else:
        #    if name not in ensg2gene.index:
        #        print name,line
        #    else:
        #        prev = gene_counts[ensg2gene.loc[name,1].upper()]
        #        gene_counts[ensg2gene.loc[name,1].upper()] = max(prev,int(float(value)))

    run.stdout.close()
    return gene_counts

if __name__ == "__main__":
    #main(sys.argv[1],sys.argv[2],sys.argv[3])
    main(*sys.argv[1:])