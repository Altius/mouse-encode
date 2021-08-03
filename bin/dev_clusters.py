#!/usr/bin/env python
import motifs,os,pandas as pd,numpy as np,subprocess
from collections import defaultdict,OrderedDict
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats
## generate the motif data for the developmental clusters
##

### quick and dirty stats: hypergeometric for motif values
### 
MASTER_LIST = '/home/jlazar/proj/mouse_encode/data/final_list/master-peaks.mouse.bed4'
TSS_FILE = '/home/jlazar/annotations/gencode.vM9.mouse.all_tss.bed'
def main(tissues,metadata_file,call_file,count_file,output_dir,master_file = MASTER_LIST,
				fdr_peaks = None,rep_factors = ['Tissue','Time'],tissue_name = None,motif_locs = None,
				tss_file = TSS_FILE):
	if not os.path.exists(output_dir): os.mkdir(output_dir)
	#metadata = pd.read_table(metadata_file,index_col = 0,header = 0)
	#count_data = pd.read_table(count_file,index_col = 0,header = 0)
	#call_data = pd.read_table(call_file,index_col = 0,header = 0)
	metadata = metadata_file
	count_data = count_file
	call_data = call_file

	rep_peak_data = average_data(count_data,metadata,rep_factors)

	if tissue_name is None: tissue_name = tissues[0]

	tissue_ids = [x for x in metadata.index.values if 
	any([t in metadata.loc[x,'Tissue'].lower() for t in tissues]) and 'Glia' not in metadata.loc[x,'Tissue']]
	tissue_ids = sorted(tissue_ids,key = lambda k: float(metadata.loc[k,'Time']))

	differential_peaks,all_peaks = get_peaks(call_data,metadata,tissue_ids,rep_factors)
	dev_clusters,peak_order = make_clusters(tissues,rep_peak_data.loc[differential_peaks],metadata)

	output_base = os.path.join(output_dir,'%s.cluster-peaks.bed')
	for cl in dev_clusters:
		write_subset(dev_clusters[cl],master_file,output_base %('-'.join([str(x) for x in cl])))

	f1 = make_heatmap(rep_peak_data.loc[peak_order],dev_clusters)
	f1.savefig(os.path.join(output_dir,'%s.cluster-heatmap.pdf' %tissue_name))


	#for cl in dev_clusters:
	trans_enr_total,trans_pv_total = get_trans_enr(dev_clusters,master_file,motif_locs = motif_locs,
		background = None)
	trans_enr_total.to_csv(os.path.join(output_dir,'%s.cluster-enr.to-all.txt' %tissue_name),sep = '\t')
	trans_pv_total.to_csv(os.path.join(output_dir,'%s.cluster-pv.to-all.txt' %tissue_name),sep = '\t')

	trans_enr,trans_pv = get_trans_enr(dev_clusters,master_file,motif_locs = motif_locs,
		background = 'var')
	trans_enr.to_csv(os.path.join(output_dir,'%s.cluster-enr.to-variable.txt' %tissue_name),sep = '\t')
	trans_pv.to_csv(os.path.join(output_dir,'%s.cluster-pv.to-variable.txt' %tissue_name),sep = '\t')
		

	pk_data = rep_peak_data
	cis_enr_total,cis_pv_total = get_cis_enr(pk_data,master_file,tss_file,subset = trans_pv.index.values)
	cis_pv_total.to_csv(os.path.join(output_dir,'%s.cDHS-pv.total.txt' %tissue_name),sep = '\t')

	cl_order = []
	for cl in dev_clusters.keys():
		cl_order.append([x for x in rep_peak_data.columns.values if x[0]==cl[0] and x[1]==cl[1]][0])

	pk_data = rep_peak_data[cl_order]
	cis_enr,cis_pv = get_cis_enr(pk_data,master_file,tss_file,subset = trans_pv.index.values)
	cis_enr.to_csv(os.path.join(output_dir,'%s.cDHS-score.txt' %tissue_name),sep = '\t')
	cis_pv.to_csv(os.path.join(output_dir,'%s.cDHS-pv.txt' %tissue_name),sep = '\t')
	#cis_pv = pd.read_table(os.path.join(output_dir,'%s.cDHS-pv.txt' %tissue_name),header = 0,index_col = 0)
	#cis_enr = pd.read_table(os.path.join(output_dir,'%s.cDHS-score.txt' %tissue_name),header = 0,index_col = 0)

	### create TF list

	cluster_tfs = get_sig_tfs(cis_pv,cis_pv_total,trans_pv,trans_pv_total)
	with open(os.path.join(output_dir,'%s.cluster-tfs.txt' %tissue_name),'w') as f1:
		for cl in cluster_tfs: f1.write('%s\t%s\n' %(cl,';'.join(cluster_tfs[cl])))

	f1 = tf_heatmap(cluster_tfs,cis_enr,trans_enr)
	f1.savefig(os.path.join(output_dir,'%s.cluster-tfs.heatmap.pdf' %tissue_name))


def tf_heatmap(cluster_tfs,tf_cis,tf_trans,per_cl = 3):
	from matplotlib.colors import ListedColormap
	my_cmap = sns.blend_palette(["#333333","#F1CC81"],n_colors = 20,as_cmap = True)


	tf_order = []
	cluster_order = sorted(cluster_tfs.keys(),key = lambda k: float(k.strip('()').split(',')[1]))
	for cl in cluster_order:
		clt =  cluster_tfs[cl][0:per_cl]

		tf_order += sorted(clt,key = lambda k: (np.sum(tf_cis.loc[k].values*np.arange(tf_cis.shape[1]))/
								np.sum(tf_cis.loc[k].values)))

	f1 = plt.figure(figsize = (3.5,8))
	for i,tf in enumerate(tf_order):
	    ax = plt.subplot(len(tf_order),1,i+1)
	    hd = np.vstack([tf_trans.loc[tf].values,tf_cis.loc[tf].values]).astype(float)
	    print hd.shape
	    hd = (hd-np.mean(hd,axis = 1)[:,np.newaxis])/np.std(hd,axis = 1)[:,np.newaxis]
	    sns.heatmap(hd,ax = ax,cbar = False,yticklabels = ['Motif','TSS'],xticklabels= False,
	                cmap = my_cmap,linewidth = .5)
	    ax.yaxis.tick_right()
	    ax.set_ylabel(tf)
	plt.tight_layout()
	return f1




from statsmodels.sandbox.stats.multicomp import fdrcorrection0   
def get_sig_tfs(cis_data,cis_data_total,trans_data,trans_data_total,fdr_limit = 0.05):
	sig_tfs = {}
	trans_data.fillna(1,inplace = True)
	trans_data[np.isnan(trans_data.values)] =1 
	trans_data = trans_data.rename(columns = {x:str(x) for x in trans_data.columns.values})
	trans_data_total = trans_data_total.rename(columns = {x:str(x) for x in trans_data_total.columns.values})

	cis_data.fillna(1,inplace = True)
	cis_data[np.isnan(cis_data.values)] = 1
	cis_data = cis_data.rename(columns = {x:str(x) for x in cis_data.columns.values})
	good_tf = np.array([x for x in trans_data.index.values if x in cis_data.index])
	cis_data_total = cis_data_total.rename(columns = {x:str(x) for x in cis_data_total.columns.values})



	for i,cl in enumerate(trans_data.columns.values):
		pss,fdrs = fdrcorrection0(trans_data.loc[good_tf,cl].values)
		pss,total_fdrs = fdrcorrection0(trans_data_total.loc[good_tf,cl].values)
		cfdr = np.maximum(fdrs,total_fdrs)
		trans_enr = good_tf[cfdr<fdr_limit]
		otx = good_tf=='OTX2'
		print fdrs[otx],total_fdrs[otx],cl

		pss,cis_fdrs = fdrcorrection0(cis_data.loc[trans_enr,cl].values)
		pss,cis_fdrs_total = fdrcorrection0(cis_data_total.loc[trans_enr,cl].values)
		cis_fdrs = np.maximum(cis_fdrs,cis_fdrs_total)
		cl_sig = trans_enr[cis_fdrs<fdr_limit]
		if 'OTX2' in list(trans_enr):
			otx = list(trans_enr).index('OTX2')
			print cis_fdrs_total[otx],cis_fdrs[otx],cl

		scores = (1-cis_fdrs[cis_fdrs<fdr_limit])*(1-cfdr[cfdr<fdr_limit][cis_fdrs<fdr_limit])
		sig_tfs[cl] = cl_sig[np.argsort(-1*scores)]
	return sig_tfs



    ### link DHSs to genes for GO terms ...

def get_cis_enr(count_data,master_file,tss_file,subset = None):
	cdhs = get_cis_dhs(master_file,tss_file)
	if subset is None: 
		subset = cdhs.keys()
	else:
		subset = np.array([x for x in subset if x in cdhs])
	
	cis_enr = pd.DataFrame(index = subset,columns = count_data.columns.values)
	cis_pv = pd.DataFrame(index = subset,columns = count_data.columns.values)

	for gn in subset:
		cis_enr.loc[gn] = np.mean(np.log(count_data.loc[cdhs[gn]].values+1),axis = 0)
		for i,cl in enumerate(count_data.columns.values):
			ncl = [x for x in count_data.columns.values if x != cl]
			cell_data = np.log(count_data.loc[cdhs[gn]].values+1)
			ttest = scipy.stats.ttest_rel(cell_data[:,i],
						np.mean(np.log(count_data.loc[cdhs[gn],ncl].values+1),axis = 1))
			ttest = max(ttest[1]/2., int(ttest[0]<0))
			cis_pv.loc[gn,cl] = ttest
	return cis_enr,cis_pv



def get_trans_enr(clusters,master_file,motif_locs = None,background = None):
	if motif_locs is None:
		motif_locs = motifs.motif_locations(master_file,master_file + '.motifs',
					motif_clusters = True,by_motif = True)
	if background is None:
		all_peaks = np.arange(len(motif_locs[motif_locs.keys()[0]]))
	elif type(background) is str:
		all_peaks = reduce(np.union1d,[clusters[c] for c in clusters])
	else:
		all_peaks = background

	print len(all_peaks)
	cluster_names = motifs.cluster_names(motifs.CLUSTER_MAPPINGS_FILE)
	cluster_names = {r:r for r in motif_locs.keys()}

	enr_data = pd.DataFrame(index = [cluster_names[r] for r in motif_locs.keys()],columns = clusters.keys())
	pv_data = pd.DataFrame(index = [cluster_names[r] for r in motif_locs.keys()],columns = clusters.keys())


	for m in motif_locs:
		#mn = cluster_names[m]
		mn = m
		for cl in clusters:
			enr,pv = hypergeometric(motif_locs[m][clusters[cl]],motif_locs[m][all_peaks])
			enr_data.loc[mn,cl] = enr
			pv_data.loc[mn,cl] = pv

	return enr_data,pv_data

def hypergeometric(set_motif,background_motif):
    background_size = len(background_motif)
    sample_size = len(set_motif)

    pseudo = (np.sum(background_motif)+1)/float(background_size+1)
    enr = np.log2((((np.sum(set_motif)+1)/float(sample_size+1))+pseudo)/
                                    (((np.sum(background_motif)+1) / float(background_size + 1))+pseudo))
    pv = 1-scipy.stats.hypergeom.cdf(np.sum(set_motif)-.1,background_size,
                                       np.sum(background_motif),sample_size)
    return enr,pv


def average_data(count_data,metadata,rep_factors):
	reps = defaultdict(list)
	for r in count_data.columns.values:
	    reps[tuple(metadata.loc[r,rep_factors])].append(r)
	rep_peak_data = pd.DataFrame(index = count_data.index.values,columns = reps.keys())
	for r in reps:
	    rep_peak_data[r] = np.mean(count_data[reps[r]].values,axis = 1)
	return rep_peak_data


def get_peaks(call_data,metadata,tissue_ids,rep_factors,rep_limit = 2):
	reps = defaultdict(list)
	for r in tissue_ids:
	    reps[tuple(metadata.loc[r,rep_factors])].append(r)
	replicate_data = pd.DataFrame(index = call_data.index.values,columns = reps.keys())
	any_data = pd.DataFrame(index = call_data.index.values,columns = reps.keys())
	for r in reps:
	    replicate_data[r] = np.sum(call_data[reps[r]].values,axis = 1) >= rep_limit
	    any_data[r] = np.sum(call_data[reps[r]].values,axis = 1) >= 1

	good = np.sum(replicate_data.values,axis = 1) > 0
	sep = np.sum(any_data.values,axis = 1) < any_data.shape[1]

	return call_data.index.values[good & sep],call_data.index.values[good]


def make_clusters(tissues,count_data,metadata):

	tissue_order = [x for x in count_data.columns.values if x[0] in tissues]
	tissue_order = sorted(sorted(tissue_order,key = lambda k: tissues.index(k[0])),key = lambda k: k[1])

	time_vector = sorted([r[1] for r in tissue_order])
	time_vector = np.array([time_vector.index(r[1]) for r in tissue_order])
	tissue_vectors = np.array([tissues.index(r[0])*.01 for r in tissue_order])

	weights = count_data[tissue_order].values
	wavg = np.sum(weights*time_vector,axis = 1) / np.sum(weights,axis = 1)


	worder = np.argsort(wavg)
	max_point = np.argmax(count_data[tissue_order].values[worder],axis = 1)
	clstrs = OrderedDict()
	for tp in sorted(set(list(max_point))):
	    name = count_data[tissue_order].columns.values[tp]
	    clstrs[name] = count_data.index.values[worder][max_point==tp]

	    
	pk_scores = time_vector[max_point]+tissue_vectors[max_point]
	morder = np.argsort(pk_scores,kind = 'mergesort')
	order = count_data.index.values[worder][morder]	
	return clstrs,order


def make_heatmap(cdata,clstrs,fdr_peaks = None, figsize = (5,6),xoffset = -0.5,set_rot = 'vertical',vmin = -0.5,vmax = 2):

	zscores = cdata[clstrs.keys()]
	std_values = np.std(zscores.values,axis = 1)
	mean_values = np.mean(zscores.values,axis = 1)
	zscores = zscores.subtract(mean_values,axis = 'rows')
	zscores = zscores.divide(std_values,axis = 'rows')

	zscores = zscores.rename(columns = {c:'%s-%s' %c for c in zscores.columns.values})

	altius = sns.color_palette(["#333333", "#F1CC81"])
	from matplotlib.colors import ListedColormap
	my_cmap = ListedColormap(altius.as_hex())
	my_cmap = sns.blend_palette(["#333333","#F1CC81"],n_colors = 20,as_cmap = True)

	with sns.axes_style('white'):
	    f1 = plt.figure(figsize = figsize)
	    ax = plt.subplot(111)
	    r = sns.heatmap(zscores,cmap = my_cmap,linewidths = 0,rasterized = True,vmin = -.5,
	                vmax = 2,yticklabels = False,xticklabels = True)

	    #ax.set_ylabel('%s DHSs' %len(good_peaks),color = 'orange')
	    ax.set_ylabel('')
	    ax.set_xlabel('')

	    start = 0
	    top,bottom = ax.get_ylim()
	    ax.set_xlim([xoffset,ax.get_xlim()[1]])
	    for vls in clstrs:
	        nmb = len(clstrs[vls])
	        length = nmb
	        ax.plot([-.1,-.1],[start+0.1*length,start+.9*length],lw = 3,color = 'orange')
	        ax.text(-.15,start + .5*length,'%d DHSs' %nmb,rotation = set_rot,color = 'orange',
	                verticalalignment='center',horizontalalignment = 'right',fontsize = 8)
	        start += nmb
	    ax.text(xoffset,.5*top,'%s DHSs' %(zscores.shape[0]),rotation = 'vertical',color = 'orange',
	                verticalalignment='center',horizontalalignment = 'right',fontsize = 12)
	    plt.tight_layout()
	return f1

MASTER_LIST='final_list/master-peaks.mouse.bed4'
def write_subset(subset,bed_file,output_file):
    count = 0
    subset = set([x for x in subset])
    fout = open(output_file,'w')
    with open(bed_file) as f1:
        for line in f1:
            if count in subset: fout.write(line)
            count += 1
    fout.close()

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