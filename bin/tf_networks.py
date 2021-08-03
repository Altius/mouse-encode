#!/usr/bin/env python
import dev_clusters
from collections import defaultdict

def main(tissues,tf_set1,tf_set2,call_data,metadata,master_file,tss_file,output_file,
	motif_locs = None,rep_factors = ['Tissue','Time']):
	set_ids = {}
	for tf in tf_set1: set_ids[tf] = 1
	for tf in tf_set2: set_ids[tf] = 2

	cDHS = dev_clusters.get_cis_dhs(master_file,tss_file)
	tissue_ids = [x for x in metadata.index.values if 
	any([t in metadata.loc[x,'Tissue'].lower() for t in tissues]) and 'Glia' not in metadata.loc[x,'Tissue']]
	tissue_ids = sorted(tissue_ids,key = lambda k: float(metadata.loc[k,'Time']))
	print metadata.loc[tissue_ids]
	differential_peaks,all_peaks = dev_clusters.get_peaks(call_data,metadata,tissue_ids,rep_factors)
	print len(differential_peaks),len(all_peaks)

	if motif_locs is None:
		motif_locs = motifs.motif_locations(master_file,master_file + '.motifs',
					motif_clusters = False,by_motif = False)

	header = ['TFA','TFB','Connection-Number','Set']
	lines = ['\t'.join(header)]
	edges = get_edges(motif_locs,cDHS,differential_peaks,tf_set1 + tf_set2)
	#edges = get_edges(motif_locs,cDHS,all_peaks,tf_set1 + tf_set2)
	for e1,e2 in edges:
		c = edges[(e1,e2)]
		edge_set = '%d-%d' %(set_ids[e1],set_ids[e2])

		lines.append('\t'.join([e1,e2,'%d' %c,edge_set]))
	with open(output_file,'w') as f1:
		f1.write('\n'.join(lines))

def get_expected(tissues,tf_sets,call_data,metadata,master_file,tss_file,output_file,
	motif_locs = None,rep_factors = ['Tissue','Time']):


	cDHS = dev_clusters.get_cis_dhs(master_file,tss_file)
	tissue_ids = [x for x in metadata.index.values if 
	any([t in metadata.loc[x,'Tissue'].lower() for t in tissues]) and 'Glia' not in metadata.loc[x,'Tissue']]
	tissue_ids = sorted(tissue_ids,key = lambda k: float(metadata.loc[k,'Time']))
	differential_peaks,all_peaks = dev_clusters.get_peaks(call_data,metadata,tissue_ids,rep_factors)
	#all_cdhs = reduce(np.union1d,[cDHS[r] for r in cDHS])
	#differential_peaks = np.intersect1d(differential_peaks,all_cdhs)

	if motif_locs is None:
		motif_locs = motifs.motif_locations(master_file,master_file + '.motifs',
					motif_clusters = False,by_motif = False)

	header = ['In-Set','Out-Set','Observed','Expected','Peaks']
	lines = ['\t'.join(header)]
	for tfs1 in tf_sets:
		for tfs2 in tf_sets:
			obs,exp,pks = get_numbers(tf_sets[tfs1],tf_sets[tfs2],cDHS,all_peaks,motif_locs)
			lines.append('\t'.join([tfs1,tfs2,'%d' %obs,'%f' %exp,'%d' %pks]))
	with open(output_file,'w') as f1: f1.write('\n'.join(lines) + '\n')


def get_numbers(tf_seta,tf_setb,cDHSs,differential_peaks,motif_locs):
	dif_peaks = set([x for x in differential_peaks])

	total_pks = 0
	pk_number = 0
	for tf in tf_setb:
		pks = cDHSs[tf]
		for pk in pks:
			if pk not in dif_peaks: continue
			total_pks += 1
			#if any([motif_locs[tf2][pk] for tf2 in tf_seta]): pk_number += 1
			pk_number += np.sum([motif_locs[tf2][pk] for tf2 in tf_seta])
	exp = total_pks * np.mean(np.sum(np.column_stack([motif_locs[r][differential_peaks] for r in tf_seta]),axis = 1))
	return pk_number,exp,total_pks


import numpy as np
def get_edges(motif_locs,cDHSs,differential_peaks,tf_sets):
	dif_peaks = set([x for x in differential_peaks])
	edges = defaultdict(int)
	for tf in tf_sets:
		pks = cDHSs[tf]

		for pk in pks:
			if pk not in dif_peaks: continue
			for tf2 in tf_sets:
				if motif_locs[tf2][pk]: edges[(tf2,tf)] += 1
	return edges