#!/usr/bin/env python
import motifs
import pandas as pd,scipy.stats,numpy as np,os,subprocess,sys
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
import multiprocessing as mp


def obs_exp(network_file,set1,set2,motif_locs = None,subset = None,seed = 1000):
    dt = pd.read_table(network_file,header = None,names = ['TFA','TFB','Peaks'])

    out_edges = np.array([x in set1 for x in dt['TFA'].values],dtype = bool)
    in_edges = np.array([x in set2 for x in dt['TFB'].values],dtype = bool)

    in_number = len(flatten_pks(dt.ix[in_edges,'Peaks'].values))
    '''
    if motif_locs is None:
        out_number = len(flatten_pks(dt.ix[out_edges,'Peaks'].values))
        total_number = len(flatten_pks(dt['Peaks'].values))
        fraction = out_number/float(total_number)
    else:
        if subset is None: subset = np.ones(len(motif_locs[motif_locs.keys()[0]]),dtype = bool)
        no_motif = np.mean(np.sum(np.column_stack([motif_locs[r][subset] for r in
                                            set1]),axis = 1)==0)
        fraction = np.mean(np.sum(np.column_stack([motif_locs[r][subset] for r in
                                            set1]),axis = 1)>=1)
        print fraction,no_motif
        fraction -= no_motif
    '''
    #all_peak_scores = score_peaks(motif_locs,set1)
    all_peak_scores = score_peaks(dt, set1)
    all_cis_peaks,cis_scores = get_cis_score(dt)

    observed_pks = np.array([int(r) for r in flatten_pks(dt.ix[in_edges & out_edges,'Peaks'])])

    print np.max(observed_pks),len(all_peak_scores),len(cis_scores)
    permutation_values = permutation_expectation(dt,observed_pks,all_cis_peaks,cis_scores,all_peak_scores,nsim = 500,
                                                 subset = observed_pks,seed = seed)

    #observed = len(flatten_pks(dt.ix[in_edges & out_edges,'Peaks']))
    #expected = in_number * fraction
    expected = np.mean(permutation_values)
    observed = np.sum(all_peak_scores[observed_pks])

    return observed,expected,(observed-expected+1)/float(expected + 1)

def score_peaks(dt,tfs):
    all_pks = [int(x) for x in flatten_pks(dt['Peaks'].values)]
    scores = np.zeros(np.max(all_pks)+1)
    for x in dt.index.values:
        pks = np.array([int(r) for r in dt.loc[x,'Peaks'].split(';')],dtype = int)
        scores[pks] += 1
    return scores

    #return np.sum(np.column_stack([motif_locs[r] for r in tfs]),axis = 1)

def get_cis_score(dt):
    pk_dict = {}
    for x in dt.index.values:
        for pk in dt.loc[x,'Peaks'].split(';'):
            if int(pk) not in pk_dict: pk_dict[int(pk)] = set([])
            pk_dict[int(pk)].add(dt.loc[x,'TFB'])
    pk_ids = np.array(sorted(pk_dict.keys()))
    cis_scores = np.zeros(np.max(pk_ids)+1)
    cis_scores[pk_ids] = np.array([len(pk_dict[x]) for x in pk_ids])
    return pk_ids,cis_scores

def permutation_expectation(dt,score_peaks,cis_peaks,cis_scores,peak_scores,nsim = 1000,subset = None,
                            seed = 1000):

    perm_values = []
    if subset is None: subset = peak_scores>0
    r = np.random.RandomState(seed)

    for x in xrange(nsim):
        rpk = np.zeros(len(peak_scores))
        rpk[subset] = r.permutation(peak_scores[subset])
        perm_values.append(np.sum(rpk[score_peaks]*cis_scores[score_peaks]))
    return perm_values

def flatten_pks(pks_lines):
    vls = []
    for x in pks_lines:
        vls.extend(x.split(';'))
    return vls

def main(metadata_file,master_file,motif_density,peak_signal,peak_calls,tss_file,
         mouse2human_file,output_dir,fdr_limit = 0.05,distance = 5000):
    if not os.path.exists(output_dir): os.mkdir(output_dir)
    metadata = pd.read_table(metadata_file,index_col = 1,header = 0)

    motif_locs = motifs.motif_locations(master_file,master_file + '.motifs',
                                        motif_clusters=False)
    peak_signal = pd.read_table(peak_signal,index_col = 0,header = 0)
    peak_signal = peak_signal.rename(columns = {c:'DS'+c.split('DS')[-1]
                                                for c in peak_signal.columns.values})
    peak_signal = combined_replicate_data(peak_signal,metadata)
    motif_density = pd.read_table(motif_density,index_col= 0,header = 0)
    motif_density = motif_density.rename(columns = {c:'DS'+c.split('DS')[-1]
                                                for c in motif_density.columns.values})
    motif_density = combined_replicate_data(motif_density,metadata)

    output_base = os.path.join(output_dir,'%s.correlated_dhs.bed')
    motif_locs = adjust_motifs(motif_locs,motif_density,peak_signal,master_file,output_base,fdr_limit = fdr_limit)

    cDHSs = get_tf_cdhs(master_file,tss_file,mouse2human_file,distance = distance)

    full_network = create_network(motif_locs,cDHSs)
    with open(os.path.join(output_dir,'full_network.txt'),'w') as f1:
        f1.write('\n'.join(full_network)+'\n')


    peak_calls = pd.read_table(peak_calls,index_col = 0,header = 0)
    peak_calls = peak_calls.rename(columns = {c:'DS'+c.split('DS')[-1]
                                                for c in peak_calls.columns.values})
    peak_calls = combined_replicate_data(peak_calls,metadata,func='replicated')
    for cl in peak_calls.columns.values:
        subset = peak_calls[cl].values.astype(bool)
        ntw = create_network(motif_locs,cDHSs,subset = subset)
        with open(os.path.join(output_dir,'%s_network.txt' %cl),'w') as f1:
            f1.write('\n'.join(ntw)+'\n')


def get_tf_cdhs(master_file,tss_file,mouse2human,distance = 5000):
    output = {}
    mouse2human = pd.read_table(mouse2human,index_col = 1,header = 0)
    mouse2human = mouse2human[~mouse2human.index.duplicated(keep='first')]
    p1 = subprocess.Popen(['bedmap','--range','%d' %distance,'--echo','--echo-map-id-uniq',tss_file,master_file],
                          stdout=subprocess.PIPE)
    for line in iter(p1.stdout.readline,''):
        ln,pks = line.strip().split('|')
        pks = [int(x) for x in pks.split(';') if x]
        if not pks: continue
        #pk = int(ln.split()[3])
        gn = ln.split()[3]
        if gn not in mouse2human.index.values:
            print gn
            continue
        gn_name = mouse2human.loc[gn,'Human']
        if pd.isnull(gn_name): gn_name = gn.upper()
        if gn_name not in output: output[gn_name] = set([])
        output[gn_name].update(pks)

    return {x:np.array(list(output[x]),dtype = int) for x in output}

def create_network(motif_locs,cDHSs,subset = None):
    lines = []
    if subset is None: subset = np.ones(len(motif_locs[motif_locs.keys()[0]])).astype(bool)
    for tf1 in motif_locs:
        for tf2 in motif_locs:
            if tf2 not in cDHSs: continue
            pks = cDHSs[tf2][(motif_locs[tf1][cDHSs[tf2]] & subset[cDHSs[tf2]])]
            if len(pks): lines.append('\t'.join([tf1,tf2,';'.join(['%d' %x for x in pks])]))
    return lines

import multiprocessing as mp
def adjust_motifs(motif_locs,motif_density,signal,master_list,output_base,fdr_limit = 0.05,
                  cpu_number = 16):
    signal = np.log(signal+1)[motif_density.columns.values]

    for tf in motif_locs:
        good = np.where(motif_locs[tf])[0]
        if tf not in motif_density.index.values:
            motif_locs[tf] = np.zeros(len(motif_locs[tf]),dtype = bool)
            continue

        splits = np.array_split(good,cpu_number*3)
        pool = mp.Pool(processes=min(cpu_number, mp.cpu_count() - 8))
        results = [pool.apply_async(get_correlations, args=(motif_density.loc[tf].values,signal.loc[s].values))
                   for s in splits]
        pool.close()
        pool.join()
        output = [s.get() for s in results]
        pvs = np.hstack(output)
        print len(pvs),len(good),pvs.shape,good.shape


        #pvs = []
        #for pk in good:
        #    pvs.append(scipy.stats.pearsonr(motif_density.loc[tf].values,
        #                                    signal.loc[pk,motif_density.columns.values])[1])

        pss,fdrs = fdrcorrection0(pvs)
        print len(fdrs),fdrs.shape
        motif_locs[tf] = np.zeros(len(motif_locs[tf]))
        motif_locs[tf][good[fdrs<fdr_limit]] = 1
        motif_locs[tf] = motif_locs[tf].astype(bool)
        write_subset(good[fdrs<fdr_limit],master_list,output_base %tf)


    return motif_locs

def get_correlations(vector,matrix):
    pvs = []
    for i in xrange(matrix.shape[0]):
        cor,pv = scipy.stats.pearsonr(vector,matrix[i])
        pvs.append(pv)
    return np.array(pvs)



def write_subset(subset,bed_file,output_file):
    count = 0
    subset = set([x for x in subset])
    fout = open(output_file,'w')
    with open(bed_file) as f1:
        for line in f1:
            if count in subset: fout.write(line)
            count += 1
    fout.close()

def combined_replicate_data(dataframe,metadata,columns = ['Tissue','Time'],func = np.mean):
    replicate_ids = {}
    for ds in metadata.index.values:
        name = ';'.join([str(metadata.loc[ds,x]) for x in columns])
        if name not in replicate_ids: replicate_ids[name] = []
        correct_ds = [x for x in dataframe.columns.values if ds==x or ds[0:7]==x[0:7]]
        if len(correct_ds)!=1: print ds,correct_ds
        replicate_ids[name].append(sorted(correct_ds)[-1])
    outd = pd.DataFrame(columns = replicate_ids,index = dataframe.index.values)
    for nm in replicate_ids:
        if func == 'replicated':
            vls = np.sum(dataframe[replicate_ids[nm]].values,axis = 1)
            outd[nm] = vls >= max(1,min(2,len(replicate_ids[nm])))
        else:
            outd[nm] = func(dataframe[replicate_ids[nm]].values,axis = 1)
    return outd


if __name__ == "__main__":
    main(*sys.argv[1:])