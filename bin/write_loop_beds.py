#!/usr/bin/env python
import sys,motifs,numpy as np,os,subprocess,pandas as pd
from collections import defaultdict
DISTANCE = 5000

def main(master_file,tss_file,pk_file,node_number,output_dir):

    cDHSs = get_cis_dhs(master_file,tss_file)
    fdict = {'1': get_1_node_peaks_byTF,'2':get_2_node_peaks_byTF,'3':get_3_node_peaks_byTF}
    pks = pd.read_table(pk_file)
    pks = pks[pks.columns.values[0]].values.astype(bool)

    mlocs = motifs.motif_locations(master_file,master_file + '.motifs',
                                    motif_clusters = False)

    if node_number == 'all':
        tf2pks = {tf:reduce(np.union1d,[cDHSs[x][pks[cDHSs[x]] & mlocs[tf][cDHSs[x]]] for x in mlocs if x in cDHSs])
                  for tf in mlocs if tf in cDHSs}
    else:
        fnc = fdict[node_number]
        tf2pks = fnc(pks,cDHSs,mlocs)
    for tf in tf2pks:
        if len(tf2pks[tf])>0:
            write_subset(tf2pks[tf],master_file,
                         os.path.join(output_dir,'%s' %tf))


def write_subset(subset, bed_file, output_file):
    count = 0
    subset = set([x for x in subset])
    fout = open(output_file, 'w')
    with open(bed_file) as f1:
        for line in f1:
            if count in subset: fout.write(line)
            count += 1
    fout.close()


def get_1_node_peaks_byTF(pks,cDHSs,motifs):
    return {tf:cDHSs[tf][pks[cDHSs[tf]] & motifs[tf][cDHSs[tf]]] for tf in motifs if tf in cDHSs}

def get_2_node_peaks_byTF(pks,cDHSs,motifs):
    reg_tfs = {tf:set([x for x in motifs if np.sum(pks[cDHSs[tf]] & motifs[x][cDHSs[tf]])>0 and x!=tf
                       and x in cDHSs])
               for tf in motifs if tf in cDHSs}

    regs = {tf:reduce(np.union1d,[cDHSs[tf2][pks[cDHSs[tf2]] & motifs[tf][cDHSs[tf2]]]
                                  for tf2 in reg_tfs[tf] if tf2 in cDHSs and tf2 !=tf])
                    for tf in reg_tfs if len(reg_tfs[tf])>0}
    return regs

def get_3_node_peaks_byTF(pks,cDHSs,motifs):
    regs = {tf:{tf2:cDHSs[tf2][pks[cDHSs[tf2]] & motifs[tf][cDHSs[tf2]]]
                for tf2 in motifs if tf2 in cDHSs} for tf in motifs if tf in cDHSs}
    loop_pks = defaultdict(set)
    for tf1 in regs:
        for tf2 in regs[tf1]:
            if tf1 == tf2: continue
            for tf3 in regs[tf2]:
                if tf3==tf2 or tf3==tf1: continue
                if tf1 in regs[tf3] and len(regs[tf3][tf1])>0:
                    loop_pks[tf3].update(list(regs[tf3][tf1]))
                    loop_pks[tf2].update(list(regs[tf2][tf3]))
                    loop_pks[tf1].update(list(regs[tf1][tf2]))
    return loop_pks

def get_cis_dhs(bed_file,tss_file,distance = DISTANCE):
    gene_counts = defaultdict(set)

    run = subprocess.Popen(['bedmap', '--range','%d' %distance,'--echo','--echo-map-id', bed_file, tss_file],
                           stdout=subprocess.PIPE)
    for line in iter(run.stdout.readline, ''):
        ref, value = line.strip().split('|')
        if not value: continue
        id = int(ref.split('\t')[3])
        for r in value.split(';'):
            gene_counts[r.upper()].add(id)
    run.stdout.close()
    return {g: np.array(list(gene_counts[g])) for g in gene_counts}

if __name__ == "__main__":
    args = sys.argv[1:]
    print args
    main(*args)

