#!/usr/bin/env python

#SBATCH --mem=8G                 # memory pool for all cores
#SBATCH -e slurm.%N.%j.err        # STDERR

import subprocess,tempfile,os,shutil,sys

MAX_FILES = 500


def main(network_dir,output_base,phylop_scores):
    if not os.path.exists(output_base): os.mkdir(output_base)
    fp_base = os.path.join(output_base,'%d-node.ffl.footprints.bed')
    cons_base = os.path.join(output_base, '%d-node.ffl.footprints.avg-conservation.txt')

    all_fp  = os.path.join(output_base,'all.footprints.bed')
    fp_file = os.path.join(network_dir,'motifs.in.universe')
    merge_files([fp_file],all_fp)
    all_conservation = os.path.join(output_base,'all.footprints.avg-conservation.txt')

    conservation_subset = os.path.join(output_base,'conservation.in-fp.bed')
    with open(conservation_subset,'w') as f1:
        p1 = subprocess.Popen(['bedmap','--faster','--echo','--skip-unmapped','--bp-ovr','1',phylop_scores,all_fp],
                              stdout = f1)
        p1.communicate()

    get_avg_conservation(all_fp,conservation_subset,all_conservation)

    for loop in [1,2,3]:
        output_file = fp_base %loop
        get_loop_fp(network_dir,loop,output_file)
        get_avg_conservation(output_file,conservation_subset,cons_base %loop)
    #os.remove(conservation_subset)

def get_loop_fp(network_dir,loop_number,output_file):
    tmp_outdir = tempfile.mkdtemp()

    network_file = os.path.join(network_dir,'cleaned.network')
    network = parse_network_file(network_file)

    footprint_file = os.path.join(network_dir,'motifs.in.universe.mapped.2.genes')
    old_tss_file = os.path.join(os.path.split(network_dir)[0],'all.tss.buff')
    tss_file = tempfile.mktemp(dir = tmp_outdir)
    shutil.copy(old_tss_file,tss_file)


    if loop_number == 1:
        one_nodes = get_1node(network)
        tf_dict = {tf:[tf] for tf in one_nodes}
    elif loop_number == 2:
        two_nodes = get_2node(network)
        tf_dict = {}
        for e1,e2 in two_nodes:
            if e2 not in tf_dict: tf_dict[e2] = set([])
            tf_dict[e2].add(e1)
    elif loop_number == 3:
        three_nodes = get_3node(network)
        tf_dict = {}
        for e1,e2 in three_nodes:
            if e2 not in tf_dict: tf_dict[e2] = set([])
            tf_dict[e2].add(e1)
    else:
        raise

    tf_output = os.path.join(tmp_outdir,'%s.temp-merged.fp.bed')
    get_footprints( tf_dict, tss_file, footprint_file, output_file,name = '%d' %loop_number)
    os.remove(tss_file)


def parse_network_file(network_file):
    network = {}
    with open(network_file) as f1:
        for line in f1:
            data = line.strip().split()
            data[0] = data[0].upper()
            data[1] = data[1].upper()
            if data[1] not in network: network[data[1]] = set([])
            if data[0] not in network: network[data[0]] = set([])
            network[data[1]].add(data[0])
    return network

def get_1node(network):
    return set([x for x in network if x in network[x]])

def get_2node(network):
    return set([(x,y) for x in network for y in network[x] if x in network[y] and x!=y])

def get_3node(network):
    edges = set([])
    for tf1 in network:
        for tf2 in network[tf1]:
            if tf2 == tf1: continue
            for tf3 in network[tf2]:
                if tf3 == tf2 or tf3 == tf1: continue
                if tf1 in network[tf3]:
                    edges.add((tf1,tf2))
                    edges.add((tf2, tf3))
                    edges.add((tf3, tf1))
    return edges

def merge_files(files,output_file):
    with open(output_file,'w') as f1:
        p1 = subprocess.Popen(['bedops','-m']+files,stdout = f1)
        p1.communicate()

def get_footprints(tf_dict,tss_file,footprint_file,output_file,name = ''):
    upstream_genes = {tf: set([x for x in tf_dict[tf]]) for tf in tf_dict}
    get_fps = ['bedmap','--echo','--echo-map','--bp-ovr','1',tss_file,footprint_file]
    sort_elements = ['sort-bed','-']
    merge_fps = ['bedops', '-m', '-']
    f1 = open(output_file,'w')

    p1 = subprocess.Popen(get_fps,stdout = subprocess.PIPE)
    p2 = subprocess.Popen(sort_elements,stdin = subprocess.PIPE,stdout = subprocess.PIPE)
    p3 = subprocess.Popen(merge_fps,stdin = p2.stdout,stdout = subprocess.PIPE)
    for line in iter(p1.stdout.readline,''):
        element,mapped = line.strip().split('|')
        nm = element.split()[3]
        if nm not in tf_dict: continue

        mapped_elements = mapped.split(';')
        for m in mapped_elements:
            if not m: continue
            dt = m.strip().split()
            if dt[4] in tf_dict[nm]:
                p2.stdin.write('\t'.join(dt[0:3])+'\n')
    p1.stdout.close()
    p2.stdin.close()
    for line in iter(p3.stdout.readline,''):
        data = line.strip().split()
        f1.write('\t'.join(data+[name])+'\n')
    p2.stdout.close()
    p3.stdout.close()

    f1.close()


def get_avg_conservation(bed_file,phylop_file,output_file):
    avg = ['bedmap','--faster','--mean','--bp-ovr','1',bed_file,phylop_file]
    with open(output_file,'w') as f1:
        p1 = subprocess.Popen(avg,stdout = f1)
        p1.communicate()

if __name__ == "__main__":
    get_avg_conservation(*sys.argv[1:])
    #network_dir = sys.argv[1]
    #output_dir = sys.argv[2]
    #scores = sys.argv[3]
    #main(network_dir,output_dir,scores)