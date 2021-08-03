#!/usr/bin/env python
import numpy as np,pandas as pd,sys,os
import statsmodels.api as sm
from statsmodels.stats.api import anova_lm
from statsmodels.formula.api import ols
import multiprocessing as mp
import cell_data


def main(metadata_file,count_file,output_file):
    metadata = pd.read_table(metadata_file,index_col = 0,header = 0)
    counts = pd.read_table(count_file,index_col = 0,header = 0)

    log_counts = np.log(counts + 1)
    full_design,restricted_design = build_design(metadata,normal_values = ['Time','Tissue'],
                                                 interaction_values =['Time','Tissue'])
    time_design,restricted_design = build_design(metadata, normal_values=['Time'])
    tissue_design, restricted_design = build_design(metadata, normal_values=['Tissue'])

    out_data = pd.DataFrame(index = log_counts.index.values,
                            columns = ['Average','Full-Pvalue','Tissue-Pvalue','Time-Pvalue'])
    avg_values = np.mean(log_counts.values,axis = 1)
    out_data['Average'] = avg_values

    colnames = ['Full-Pvalue','Tissue-Pvalue','Time-Pvalue']
    for jj,bk_design in enumerate([restricted_design,time_design,tissue_design]):
        pool = mp.Pool(processes = max(32,mp.cpu_count()-8))
        results = [pool.apply_async(run_anova,args=(full_design,bk_design,log_counts.loc[i].values))
                    for i in xrange(out_data.shape[0])]
        pool.close()
        pool.join()
        pvalues =  [s.get() for s in results]
        out_data[colnames[jj]] = pvalues
        out_data.to_csv(output_file,sep = '\t')
    #for i in xrange(log_counts.shape[0]):
    #    pvalue = run_anova(full_design,restricted_design,log_counts.iloc[i].values)
    #    pvalues.append(pvalue)
    #    avg_values.append(np.mean(log_counts.iloc[i].values))
    #return np.array(pvalues),avg_values
    #return np.array(pvalues),np.array(avg_values)

def run_anova(full_design_matrix,restricted_design_matrix,yvalues):
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
    ## should do the normalize here
    metadata,count_data,call_data,output_dir = sys.argv[1:]
    raw_matrix = pd.read_table(count_data,index_col = 0,header = 0)
    peak_matrix = pd.read_table(call_data,index_col = 0,header = 0)

    tmm_output = os.path.join(output_dir,'tmm-norm.txt')
    tmm_norm = cell_data.normalize(raw_matrix,peak_matrix)
    tmm_norm.to_csv(tmm_output,sep = '\t')
    main(metadata,tmm_output,os.path.join(output_dir,'peak-anova.tmm.txt'))

    q_output = os.path.join(output_dir,'quantile-norm.txt')
    qnorm = cell_data.qnorm(raw_matrix)
    qnorm.to_csv(q_output,sep = '\t')
    main(metadata,q_output,os.path.join(output_dir,'peak-anova.qnorm.txt'))

