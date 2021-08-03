#!/usr/bin/env python
import pandas as pd,subprocess,tempfile,os,time
## goal go from a spreadsheet of DS #s/peak files to all

MAP_MOTIFS = os.path.join(os.path.split(__file__)[0],'map_motifs.sh')
MOTIF_DENSITY_CMD = os.path.join(os.path.split(__file__)[0],'peak_motif_density.py')

def main(experiment_data,output_dir):
    if not os.path.exists(output_dir): os.mkdir(output_dir)
    data = pd.read_table(experiment_data,index_col = 0,header = 0)
    all_peak_files = data['Peak-File'].values

    merge_file = os.path.join(output_dir,'merged-peaks.bed')
    merge_peak_files(all_peak_files,merge_file)

    p1 = subprocess.Popen([MAP_MOTIFS,merge_file])
    p1.communicate()

    motif_file = merge_file + '.motifs'
    while not os.path.exists(motif_file):
        time.sleep(600)

    motif_dir = os.path.join(output_dir,'sample-motif-densities')
    if not os.path.exists(motif_dir): os.mkdir(motif_dir)
    for i,pk in enumerate(all_peak_files):
        id = data['DS-ID'].values[i]+'-%s' %(str(data['Aggregation-ID'].values[i]))
        submit_motif_density(pk,motif_file,'motif_dir/%s.motif-density.txt' %id)

    summary_file = os.path.join(output_dir,'motif-density-summary.txt')
    combined_data = summarize(motif_dir)
    combined_data.to_csv(summary_file,sep = '\t')

def summarize(data_dir):
    data = None
    for fl in os.listdir(data_dir):
        if 'motif-density' in fl:
            nm = fl.split('.')[0]
            fl_data = pd.read_table(os.path.join(data_dir,fl),index_col = 0,header = None)
            if data is None:
                data = pd.DataFrame(fl_data.values[:,1],columns = [nm],index = fl_data.index.values)
            else:
                data[nm] = fl_data.loc[data.index.values]
    return data

def submit_motif_density(peak_file,motif_file,output_file):

    cmd = ['sbatch', '--job-name=' + os.path.split(output_file)[1], MOTIF_DENSITY_CMD, peak_file, motif_file,
               output_file]
    ssh = ['ssh', 'sched0', '\'' + ' '.join(cmd) + '\'']
    p1 = subprocess.Popen(' '.join(ssh), shell=True)
    p1.communicate()



def merge_peak_files(peak_files,output_file,max_number = 500):
    if len(peak_files) <= max_number:
        with open(output_file,'w') as f1:
            p1 = subprocess(['bedops','-m']+[x for x in peak_files],stdout = f1)
            p1.communicate()

    else:
        tmp_files = []
        for i in range(0,len(peak_files),max_number):
            tmpf = tempfile.mkstemp()
            merge_peak_files(peak_files[i:i+max_number],tmpf)
            tmp_files.append(tmpf)
        merge_peak_files(tmp_files,output_file)
        for tmpf in tmp_files:
            os.remove(tmpf)


def parsing_human_ds_numbers():
    all_encode_file = '/home/jlazar/proj/mouse_encode/bin/lims/human.encode.datasets.txt'
    curated_fetal_time_points = '/home/jlazar/proj/mouse_encode/bin/lims/fetal_human_samples.txt'
    combined_file = '/home/jlazar/proj/mouse_encode/bin/lims/human_motif_density_datasets.txt'
    encode2dir = '/net/seq/data/data-release/all_DNase_peaks'

    peak_loc_data = pd.read_table(all_encode_file,index_col = 0,header = 0)
    curated_data = pd.read_table(curated_fetal_time_points,index_col = 0,header = 0)

    for fl in os.listdir(encode2dir):
        try:
            nm,ds = fl.split('.')[0].split('-')
        except ValueError:
            print fl
            continue
        if not any([ds in x for x in curated_data.index.values]):
            print ds,nm
            curated_data.loc[ds] = [nm,'nan','nan']
    for ds in curated_data.index.values:
        if ds not in peak_loc_data.index.values:
            good = [x for x in peak_loc_data.index.values if x[0:7]==ds]
            if not good:
                print ds
                continue
            else:
                ds2 = sorted(good)[-1]
        else:
            ds2 = ds
        curated_data.loc[ds,'Peak-File'] = peak_loc_data.loc[ds2,'Peak-File']
        curated_data.loc[ds, 'Cut-Count-File'] = peak_loc_data.loc[ds2, 'Cut-Count-File']
        curated_data.loc[ds,'Aggregation-ID'] = peak_loc_data.loc[ds2,'Aggregation-ID']
    curated_data.to_csv(combined_file,sep = '\t')


