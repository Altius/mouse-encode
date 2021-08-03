#!/usr/bin/env python
import human2mouse_chain
import dev_clusters
import pandas as pd,numpy as np,os

def create_summary(outdir):
    os.chdir('/home/jlazar/proj/mouse_encode/data')
    match_file = 'mouse2human_cor_match.txt'
    sorted_mouse = 'human2mouse_mapping/lifted.over.human.to.mouse.txt.sort.bed'
    sorted_human = 'human2mouse_mapping/lifted.over.mouse.to.human.txt.sort.bed'
    master_file = 'final_list/master-peaks.mouse.bed4'

    conserve_data = np.loadtxt(master_file + '.conservation')
    goodc = np.where(~np.isnan(conserve_data))[0]

    call_file = 'nov_17/batch_library/peak-ind.txt'
    mdata = pd.read_table('mouse-final-samples.full.txt',index_col = 1,header = 0)
    human_metadata = pd.read_table('human_gini/mouse-encode.human-samples.full.txt',
        index_col = 6,header = 0,delimiter = '\t')

    call_data = pd.read_table(call_file,index_col = 0,header = 0)
    rep_calls,cm_data = combined_replicate_data(call_data,mdata,rep_limit = 2)

    if not os.path.exists(outdir): 
        os.mkdir(outdir)
        os.mkdir(os.path.join(outdir,'mouse_bed'))
        os.mkdir(os.path.join(outdir,'human_bed'))
    #motif_locs = motifs.motif_locations(master_file,master_file + '.motifs',
    #                                by_motif = False,motif_clusters = False)

    tss_file = '/home/jlazar/annotations/gencode.vM9.mouse.all_tss.bed'
    cDHSs = dev_clusters.get_cis_dhs(master_file,tss_file,distance = 5000)
    
    matches = [tuple(line.strip().split('\t')) for line in open(match_file).readlines()]
    print matches
    summary_dt = pd.DataFrame(index = np.arange(len(matches)))
    count = 0
    for mt,ht in matches:
        if ',' in ht: ht = ht.split(',')[0]
        summary_dt.loc[count,'Mouse-Sample'] = mt
        summary_dt.loc[count,'Human-Sample'] = ht

        if mt not in rep_calls.columns.values: 
            print '%s odes not exist' %mt
            continue

        mouse_pk_file = os.path.join(outdir,'mouse_bed','%s.rep-peaks.bed' %mt) 
        all_pks = np.where(rep_calls[mt].values)[0]
        
        ## have the peaks for the networks
        write_subset(all_pks,master_file,mouse_pk_file)
        try:
            human_pk_file = [human_metadata.loc[r,'Peak-File'] for r in human_metadata.index.values
                             if human_metadata.loc[r,'Name'] == ht.split(';')[0] 
                             and (human_metadata.loc[r,'Age']==ht.split(';')[1] or ht.split(';')[1]=='nan')][0]
        except IndexError:
            print 'Tried %s-%s' %(mt,ht)
            continue
        try:
            os.symlink(human_pk_file,os.path.join(outdir,'human_bed',os.path.split(human_pk_file)[1]))
        except OSError:
            pass

        ## parse all
        recip_file = os.path.join(outdir,'%s.reciprocal.bed' %(mt))
        human2mouse_chain.get_reciprocal_peak(sorted_mouse,sorted_human,recip_file,mouse_file = mouse_pk_file,
                        human_file = human_pk_file)

        all_pks = rep_calls[mt].values.astype(bool)
        one_node = get_1node_peaks(all_pks,cDHSs,motif_locs)
        two_node = get_2node_peaks(all_pks,cDHSs,motif_locs)
        #wt2a = np.zeros(len(all_pks))
        #for pk in wt2: wt2a[pk] = wt2[pk]
        three_node = get_3node_peaks(all_pks,cDHSs,motif_locs)
        cpeaks = get_all_peaks(all_pks,cDHSs,motif_locs)

        names = ['One-Node','Two-Node','Three-Node','All-TF','All']
        for i,pkset in enumerate([one_node,two_node,three_node,cpeaks,np.where(all_pks)[0]]):
            subset_peaks = os.path.join(outdir,'%s.%s.mouse-peaks.bed' %(mt,names[i]))
            write_subset(pkset,master_file,subset_peaks)

            recip_file = os.path.join(outdir,'%s.%s-reciprocal.bed' %(mt,names[i]))
            human2mouse_chain.get_reciprocal_peak(sorted_mouse,sorted_human,recip_file,mouse_file = subset_peaks,
                            human_file = human_pk_file)
            #if i == 1:
            gp = np.intersect1d(goodc,pkset)
            #summary_dt.loc[count,'%s-conservation' %names[i]] = np.sum(conserve_data[gp]*wt2a[gp])/np.sum(wt2a[gp])
            #else:
            summary_dt.loc[count,'%s-conservation' %names[i]] = np.mean(conserve_data[np.intersect1d(goodc,pkset)])
            summary_dt.loc[count,'%s-reciprocal-peaks' %names[i]] = len(open(recip_file).readlines())/float(len(pkset))


        count += 1
    ofile = os.path.join(outdir,'convseration-summary.txt')
    summary_dt.to_csv(ofile,sep = '\t',index = False)

import subprocess
def overlap_fraction(subset_file,other_file):
    p1 = subprocess.Popen(['bedmap','--indicator',subset_file,other_file],
                          stdout = subprocess.PIPE)
    vls = []
    for line in iter(p1.stdout.readline,''):
        vls.append(int(line.strip()))
    p1.stdout.close()
    return np.mean(vls)

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


def main(human_network_file,mouse_network_file,output_file):
    human_network = parse_network_file(human_network_file)
    mouse_network = parse_network_file(mouse_network_file)

    lines = ['Loop\tHuman\tMouse\tShared']
    h_total,m_total,shared_total = get_jaccard_index(human_network,mouse_network)
    lines.append('Total\t%d\t%d\t%d' %(h_total,m_total,shared_total))
    for n in range(1,4,1):
        h_num,m_num,shared_num  = get_shared_autoreg(human_network,mouse_network,n)
        lines.append(['%d\t%d\t%d\t%d'] %(n,h_num,m_num,shared_num))

def get_shared_autoreg(network1,network2,node):
    if node == 1:
        edges1 = get_1node(network1)
        edges2 = get_1node(network2)
    elif node == 2:
        edges1 = get_2node(network1)
        edges2 = get_2node(network2)
    elif node == 3:
        edges1 = get_3node(network1)
        edges2 = get_3node(network2)
    else:
        raise ValueError
    return len(edges1),len(edges2),len(edges1 & edges2)

import numpy as np
def get_1node(network):
    return set([x for x in network if x in network[x]])

def get_1_node_peaks_byTF(pks,cDHSs,motifs):
    return {tf:cDHSs[tf][pks[cDHSs[tf]] & motifs[tf][cDHSs[tf]]] for tf in motifs if tf in cDHSs}

def get_2_node_peaks_byTF(pks,cDHSs,motifs):
    regs = {tf:reduce(np.union1d,[cDHSs[tf2][pks[cDHSs[tf2]] & motifs[tf][cDHSs[tf2]]]
                                  for tf2 in motifs if tf2 in cDHSs and tf2 !=tf])
                    for tf in motifs if tf in cDHSs}
    return regs

def get_3_node_peaks_byTF(pks,cDHSs,motifs):
    regs = {tf:[(tf2,cDHSs[tf2][pks[cDHSs[tf2]] & motifs[tf][cDHSs[tf2]]])
                for tf2 in motifs if tf2 in cDHSs] for tf in motifs if tf in cDHSs}
    loop_pks = defaultdict(set)
    for tf1 in regs:
        for tf2 in regs[tf1]:
            if tf1 == tf2: continue
            for tf3 in regs[tf2]:
                if tf3==tf2 or tf3==tf1: continue
                if tf1 in regs[tf3] and len(regs[tf3][tf1][1])>0:
                    loop_pks[tf1].update(list(regs[tf1][tf2][1]))
    return loop_pks



def get_1node_peaks(pks,cDHSs,motifs):
    return reduce(np.union1d,[cDHSs[tf][pks[cDHSs[tf]] & motifs[tf][cDHSs[tf]]]
                            for tf in motifs if tf in cDHSs])

def get_all_peaks(pks,cDHSs,motifs):
    return reduce(np.union1d,[cDHSs[tf][pks[cDHSs[tf]]]
                            for tf in motifs if tf in cDHSs])


def get_2node(network):
    return set([(x,y) for x in network for y in network[x] if x in network[y] and x!=y])

from collections import defaultdict
def get_2node_peaks(pks,cDHSs,motifs):
    pk_weights = defaultdict(int)
    regs = {tf:[tf2 for tf2 in motifs if tf2 in cDHSs and np.sum(pks[cDHSs[tf2]] & motifs[tf][cDHSs[tf2]])>0]
                    for tf in motifs if tf in cDHSs}
    for tf in regs:
        for tf2 in regs[tf]:
            if tf2 == tf: continue
            for pk in cDHSs[tf][pks[cDHSs[tf]] & motifs[tf2][cDHSs[tf]]]: pk_weights[pk] += 1
    #for tf in motifs:
    #    regs[tf] = [tf2 for tf2 in motifs if tf2!=tf and np.sum(motifs[tf][cDHSs[tf2]] & 
    lps = reduce(np.union1d,[cDHSs[tf][pks[cDHSs[tf]] & motifs[tf2][cDHSs[tf]]] for tf in regs
                            for tf2 in regs[tf] if tf!=tf2])#,pk_weights
    return lps

def get_3node_peaks(pks,cDHSs,motifs):
    regs = {tf:set([tf2 for tf2 in motifs if tf2 in cDHSs and np.sum(pks[cDHSs[tf2]] & motifs[tf][cDHSs[tf2]])>0])
                    for tf in motifs if tf in cDHSs}
    loop_pks = set([])
    for tf1 in regs:
        for tf2 in regs[tf1]:
            if tf1 == tf2: continue
            for tf3 in regs[tf2]:
                if tf3==tf2 or tf3==tf1: continue
                if tf1 in regs[tf3]:
                    pv = list(cDHSs[tf1][pks[cDHSs[tf1]] & motifs[tf3][cDHSs[tf1]]])
                    if len(pv)<1: raise
                    loop_pks.update(pv)
    return np.array([r for r in loop_pks])


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
def get_jaccard_index(network1,network2):
    n1_total = 0
    n2_total = 0
    shared = 0
    for tf in network1:
        n1_total += len(network1[tf])
        if tf in network2: shared += len(network1[tf] & network2[tf])
    for tf in network2:
        n2_total += len(network2[tf])
    return n1_total,n2_total,shared

def write_subset(subset,bed_file,output_file):
    count = 0
    subset = set([x for x in subset])
    fout = open(output_file,'w')
    with open(bed_file) as f1:
        for line in f1:
            if count in subset: fout.write(line)
            count += 1
    fout.close()


def combined_replicate_data(dataframe,metadata,columns = ['Tissue','Time'],rep_limit = None):
    replicate_ids = {}
    
    cm_data = {}
    for ds in metadata.index.values:
        name = ';'.join([str(metadata.loc[ds,x]) for x in columns])
        if name not in replicate_ids: replicate_ids[name] = []
        correct_ds = [x for x in dataframe.columns.values if ds==x or ds[0:7]==x[0:7]]
        #if len(correct_ds)!=1: 
        #    print ds,correct_ds
        #    continue
        replicate_ids[name].append(sorted(correct_ds)[-1])
        cm_data[name] = [metadata.loc[ds,'Tissue'],metadata.loc[ds,'System']]
    outd = pd.DataFrame(columns = replicate_ids,index = dataframe.index.values)
    for nm in replicate_ids:
        if 'forebrain' in nm: print replicate_ids[nm],nm
        if rep_limit is None:
            outd[nm] = np.mean(dataframe[replicate_ids[nm]].values,axis = 1)
        else:
            outd[nm] = np.sum(dataframe[replicate_ids[nm]].values,axis = 1)>=min(len(replicate_ids[nm]),
                                                                                 rep_limit)
        
    return outd,cm_data

import tempfile,motifs
def cell_conservation(pks,master_file,motif_file,motif_conservation_file,motif_locs,cDHSs):
    one_loop_pks = get_1_node_peaks_byTF(pks,cDHSs,motif_locs)
    two_loop_pks = get_1_node_peaks_byTF(pks, cDHSs, motif_locs)
    three_loop_pks = get_1_node_peaks_byTF(pks, cDHSs, motif_locs)

    one_conservation = combine_loops(one_loop_pks,master_file,motif_file,motif_conservation_file)
    return one_conservation

import multiprocessing as mp
def combine_loops(tf2pks,master,motif_file,motif_conservation_file):
    num = 0
    denom = 0
    pool = mp.Pool(processes=max(32, mp.cpu_count() - 8))
    results = [pool.apply_async(get_tf_conservation, args=(tf2pks[tf],tf, master,motif_file,motif_conservation_file))
               for tf in tf2pks]
    pool.close()
    pool.join()
    output = [s.get() for s in results]

    for avg,number in output:
        if number == 0: continue
        num += avg * number
        denom += number
    return num / float(denom)


def get_tf_conservation(pks,tf,master_file,motif_file,motif_conservation_file):
    if not len(pks): return np.nan,0

    motif2tf = motifs.get_tfs_from_file(motifs.MOTIF_MAPPINGS_FILE,False)
    tfile1 = tempfile.mktemp()
    tfile2 = tempfile.mktemp()

    write_subset(pks,master_file,tfile1)
    f1 = open(tfile2,'w')
    pm = subprocess.Popen(['bedops','-m','-'],stdin = subprocess.PIPE,stdout = f1)
    p1 = subprocess.Popen(['tabix','-R',tfile1,motif_file],stdout = subprocess.PIPE)
    for line in iter(p1.stdout.readline,''):
        mt = line.strip().split()[3].upper()
        if mt in motif2tf and tf in motif2tf[mt]:
            pm.stdin.write(line)
    pm.stdin.close()
    p1.stdout.close()
    pm.communicate()
    f1.close()

    vls = 0
    count = 0
    p1 = subprocess.Popen(['tabix','-R',tfile2,motif_conservation_file],stdout = subprocess.PIPE)
    for line in iter(p1.stdout.readline,''):
        mt = line.strip().split()
        try:
            cons = float(mt[4])
            vls += cons
        except:
            continue
        count += 1
    p1.stdout.close()

    os.remove(tfile1)
    os.remove(tfile2)
    return vls/float(count),count
