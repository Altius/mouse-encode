#!/usr/bin/env python
import numpy as np,pandas as pd,sys
import statsmodels.api as sm
from statsmodels.stats.api import anova_lm
from statsmodels.formula.api import ols
import multiprocessing as mp,subprocess
import motifs



def average_data(metadata,count_data,collapse_on = ['Time']):
    dtypes = defaultdict(set) #(list(metadata[collapse_on].values))
    for x in metadata.index.values:
        dtypes[';'.join(metadata.loc[x,collapse_on].values.astype(str))].add(x)
    nmeta = pd.DataFrame(index = dtypes.keys(),columns = collapse_on)
    ncounts = pd.DataFrame(index = count_data.index.values,columns = dtypes.keys())
    for k in dtypes:
        ncounts[k] = np.mean(count_data[list(dtypes[k])].values,axis = 1)
        nmeta.loc[k] = metadata.loc[list(dtypes[k])[0],collapse_on]
    return ncounts,nmeta

def main(metadata_file,norm_count_file,bed_file,tss_file,output_file):
    metadata = pd.read_table(metadata_file,index_col = 0,header = 0)
    counts = pd.read_table(norm_count_file,index_col = 0,header = 0)
    log_counts = np.log(counts + 1)
    print log_counts.shape
    log_counts,metadata = average_data(metadata,log_counts)
    #time_design,restricted_design = build_design(metadata, normal_values=['Time'])
    print log_counts.shape
    cis_peaks = get_cis_dhs(bed_file,tss_file)
    cis_peaks = {g: np.intersect1d(cis_peaks[g],log_counts.index.values) for g in cis_peaks}
    print cis_peaks.keys()[0:10]
    motif_locations = motifs.motif_locations(bed_file,bed_file +'.motifs',motif_clusters = False)


    tfs = motif_locations.keys()


    #out_data = pd.DataFrame(index = tfs,
    #                        columns = ['Motif-number','Motif-Pvalue','Cis-number','Cis-Pvalue'])

    #pool = mp.Pool(processes=min(32, mp.cpu_count() - 8))
    #out_data['Motif-number'] = [np.sum(motif_locations[tf]) for tf in out_data.index.values]
    #results = [pool.apply_async(run_anova, args=(time_design, restricted_design,
    #                                    log_counts.loc[np.where(motif_locations[tf])[0]].values))
    #               for tf in out_data.index.values if np.sum(motif_locations[tf])]
    #pool.close()
    #pool.join()
    #pvalues = [s.get() for s in results]
    pvalues = []
    rep_number = log_counts.shape[1]
    #for tf in tfs:
    for tf in ['NFIB']:
        if tf != 'NFIB': continue
        pks = cis_peaks[tf]
        pks = np.intersect1d(np.where(motif_locations[tf])[0],log_counts.index.values)
        print np.column_stack([pks]*rep_number).shape
        print np.column_stack([[metadata.loc[x, 'Time']] * len(pks) for x in log_counts.columns.values]).shape
        dt = {'DHS':log_counts.loc[pks].values.flatten(),'Peak': np.column_stack([pks]*rep_number).flatten(),
              'Time':np.column_stack([[metadata.loc[x,'Time']]*len(pks) for x in log_counts.columns.values]).flatten()}
        for d in dt: print d,dt[d].shape
        pd.DataFrame.from_dict(dt).to_csv(output_file,sep = '\t',index = False)
        return
    #return
    #out_data['Motif-Pvalue'] = pvalues
    cis_data = pd.DataFrame(index = tfs,columns = metadata.index.values)
    for tf in tfs:
        if tf not in cis_peaks or len(cis_peaks[tf]) == 0:
            continue

            cis_data.loc[tf] = 0
        else:
            cis_data.loc[tf]= np.mean(log_counts.loc[cis_peaks[tf],cis_data.columns.values],axis = 0)
    cis_data.to_csv(output_file,sep = '\t')
    #pool = mp.Pool(processes=min(32, mp.cpu_count() - 8))
    #out_data['Cis-number'] = [len(cis_peaks[tf]) for tf in out_data.index.values]
    #results = [pool.apply_async(run_anova, args=(time_design, restricted_design,
    #                                             log_counts.loc[cis_peaks[tf]].values))
    #               for tf in out_data.index.values if len(cis_peaks[tf])]
    #pool.close()
    #pool.join()
    #pvalues = [s.get() for s in results]
    #out_data['Cis-Pvalue'] = pvalues

    #out_data.to_csv(output_file, sep='\t')



from collections import defaultdict
def get_cis_dhs(bed_file,tss_file,distance = 5000):
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


def run_anova(full_design_matrix,restricted_design_matrix,yvalues):
    yvalues = (yvalues - np.mean(yvalues,axis = 1)).flatten()
    print full_design_matrix.shape,restricted_design_matrix.shape,yvalues.shape
    f1 = sm.OLS(yvalues,full_design_matrix).fit()
    fr = sm.OLS(yvalues, restricted_design_matrix).fit()
    return f1.compare_f_test(fr)[1]



def build_design(metadata,normal_values = ['Time','Tissue'],interaction_values = []):
    base = []
    for n in normal_values:
        options = set([x for x in metadata[n].values])
        base.append(np.column_stack([metadata[n].values==r for r in options]))
    for i,n1 in enumerate(interaction_values):
        for j,n2 in enumerate(interaction_values):
            if i>=j: continue
            options1 = set([x for x in metadata[n1].values])
            options2 = set([x for x in metadata[n2].values])
            base.append(np.column_stack([(metadata[n1].values==r1) & (metadata[n2].values==r2)
                                         for r1 in options1 for r2 in options2]))
    base = np.column_stack(base).astype(int)
    ones = np.ones(base.shape[0])
    return np.column_stack([base,ones]),np.array([ones])
    wtf_value = 'Batch'
    wtf_options = []
    weights = []
    for r in metadata[wtf_value].values:
        vls = r.split(';')
        wtf_options.extend(vls)
        weights.append(len(vls))
    wtf_options = set(wtf_options)

    weights = np.array(weights).astype(float)
    batch_matrix = []
    for n in wtf_options: batch_matrix.append(np.array([int(n in rv)/weights[i] for i,rv
                                               in enumerate(metadata[wtf_value].values)]))
    batch_matrix = np.column_stack(batch_matrix)
    return np.column_stack([base,batch_matrix]),batch_matrix

if __name__ == "__main__":
    main(*sys.argv[1:])