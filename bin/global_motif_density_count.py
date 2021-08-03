#!/usr/bin/env python
import pandas as pd, subprocess, tempfile, os, time, numpy as np, sys
from collections import defaultdict

## goal go from a spreadsheet of DS #s/peak files to all

MOTIF_DENSITY_CMD = os.path.join(os.path.split(os.path.realpath(__file__))[0], 'peak_motif_density_score.py')


def main(metadata_file,master_file, output_dir, do_motifs=False, motif_cluster=False, genome='hg38'):
    output_dir = os.path.realpath(output_dir)

    if not os.path.exists(output_dir): os.mkdir(output_dir)

    data = pd.read_table(metadata_file, index_col=0, header=0)



    motif_file = master_file + '.motifs'

    if not os.path.exists(motif_file): make_motif_mapping_file(master_file, genome=genome)
    while not os.path.exists(motif_file):
        time.sleep(600)

    motif_dir = os.path.join(output_dir, 'sample-motif-densities')
    if not os.path.exists(motif_dir): os.mkdir(motif_dir)

    motif_output_files = []

    for ds in data.index.values:
        if not type(data.loc[ds, 'Cut-Count-File']) is str: continue
        try:
            id = ds
        except ValueError:
            id = ds
        if os.path.exists(data.loc[ds, 'Cut-Count-File']):
            output_file = os.path.join(motif_dir, '%s.motif-density.txt' % id)
            if not os.path.exists(output_file):
                submit_motif_density(master_file, data.loc[ds, 'Cut-Count-File'], motif_file, output_file,
                                     motif_cluster=motif_cluster)
                #submit_motif_density(merge_file,data.loc[ds, 'Cut-Count-File'], motif_file, output_file, motif_cluster=motif_cluster)
                motif_output_files.append(output_file)

    while any([not os.path.exists(x) for x in motif_output_files]):
        time.sleep(600)

    if True: #motif_output_files:
        summary_file = os.path.join(output_dir, 'motif-density-summary.txt')
        combined_data = summarize(motif_dir)
        combined_data.to_csv(summary_file, sep='\t')


def make_motif_mapping_file(peakset_filename, genome='hg19'):
    dirname = os.path.split(os.path.realpath(__file__))[0]
    MAP_FILE = os.path.join(dirname, 'map_motifs_genome.sh')
    cmd = [MAP_FILE, peakset_filename, genome]
    ssh = ['ssh', 'sched0', '\'' + ' '.join(cmd) + '\'']
    p1 = subprocess.Popen(' '.join(ssh), shell=True)
    p1.communicate()


def jarch2starch(jarch_file, out_file):
    with open(out_file, 'w') as f1:
        p1 = subprocess.Popen([UNJARCH, jarch_file], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(['starch', '-'], stdin=p1.stdout, stdout=f1)
        p1.stdout.close()
        p2.wait()


def parse_bam_output(bam_out_file, tss_file, output_file, ensg=True):
    names = [x.split()[4] for x in open(tss_file).readlines()]
    gene_counts = defaultdict(int)
    with open(bam_out_file) as f1:
        for line in f1:
            count, value = line.strip().split('|')
            name = names[count]
            prev = gene_counts[name.upper()]
            gene_counts[name.upper()] = max(prev, int(float(value)))

    lines = ['%s\t%d' % (g, gene_counts[g]) for g in gene_counts]
    with open(output_file, 'w') as f1: f1.write('\n'.join(lines) + '\n')


def summarize(data_dir, file_stem='motif-density'):
    data = None
    for fl in os.listdir(data_dir):
        if file_stem in fl:
            nm = fl.split('.motif')[0]
            try:
                fl_data = pd.read_table(os.path.join(data_dir, fl), index_col=0, header=None)
            except:
                continue
            if data is None:
                data = pd.DataFrame(fl_data.values[:, 0], columns=[nm], index=fl_data.index.values)
            else:
                data[nm] = fl_data.loc[data.index.values,1]
    return data


def submit_motif_density(peak_file,cut_count_file,motif_file, output_file, motif_cluster=False):
    print motif_file,os.path.exists(motif_file)
    if motif_cluster:
        cmd = ['sbatch', '--job-name=' + os.path.split(output_file)[1], MOTIF_DENSITY_CMD, peak_file, cut_count_file,
               motif_file, output_file, 'cluster']
    else:
        cmd = ['sbatch', '--job-name=' + os.path.split(output_file)[1], MOTIF_DENSITY_CMD, peak_file, cut_count_file,
               motif_file, output_file]

    ssh = ['ssh', 'sched0', '\'' + ' '.join(cmd) + '\'']

    # t_cmd = [MOTIF_DENSITY_CMD, peak_file, motif_file,
    # output_file]
    p1 = subprocess.Popen(' '.join(ssh), shell=True)
    # p1 = subprocess.Popen(t_cmd)
    p1.communicate()


def merge_peak_files(peak_files, output_file, max_number=500):
    if len(peak_files) <= max_number:
        with open(output_file, 'w') as f1:
            count = 0
            p1 = subprocess.Popen(['bedops', '-m'] + [x for x in peak_files], stdout=subprocess.PIPE)
            for line in iter(p1.stdout.readline,''):
                f1.write(line.strip() + '\t%d\n' %count)
                count += 1
            p1.stdout.close()
    else:
        tmp_files = []
        for i in range(0, len(peak_files), max_number):
            x, tmpf = tempfile.mkstemp()
            merge_peak_files(peak_files[i:i + max_number], tmpf)
            tmp_files.append(tmpf)
        merge_peak_files(tmp_files, output_file)
        for tmpf in tmp_files:
            os.remove(tmpf)


def parsing_human_ds_numbers():
    all_encode_file = '/home/jlazar/proj/mouse_encode/bin/lims/human.encode.datasets.txt'
    curated_fetal_time_points = '/home/jlazar/proj/mouse_encode/bin/lims/fetal_human_samples.txt'
    combined_file = '/home/jlazar/proj/mouse_encode/bin/lims/human_motif_density_datasets.txt'
    encode2dir = '/net/seq/data/data-release/all_DNase_peaks'

    peak_loc_data = pd.read_table(all_encode_file, index_col=0, header=0)
    curated_data = pd.read_table(curated_fetal_time_points, index_col=0, header=0)

    for fl in os.listdir(encode2dir):
        try:
            nm, ds = fl.split('.')[0].split('-')
        except ValueError:
            print fl
            continue
        if not any([ds in x for x in curated_data.index.values]):
            print ds, nm
            curated_data.loc[ds] = [nm, 'nan', 'nan']
    for ds in curated_data.index.values:
        if ds not in peak_loc_data.index.values:
            good = [x for x in peak_loc_data.index.values if x[0:7] == ds]
            if not good:
                print ds
                continue
            else:
                ds2 = sorted(good)[-1]
        else:
            ds2 = ds
        curated_data.loc[ds, 'Peak-File'] = peak_loc_data.loc[ds2, 'Peak-File']
        curated_data.loc[ds, 'Cut-Count-File'] = peak_loc_data.loc[ds2, 'Cut-Count-File']
        curated_data.loc[ds, 'Aggregation-ID'] = peak_loc_data.loc[ds2, 'Aggregation-ID']
    curated_data.to_csv(combined_file, sep='\t')


if __name__ == "__main__":
    # metadata_file = '/home/jlazar/proj/mouse_encode/bin/lims/human_motif_density_datasets.txt'
    # output_dir = '/home/jlazar/proj/mouse_encode/data/human_comp'

    # metadata_file = '/home/jlazar/lebowski-mirror/LSD1/results-12-16/lineage_MDS/all-encode.4motif.txt'
    # output_dir = '/home/jlazar/lebowski-mirror/LSD1/results-12-16/lineage_MDS/motifs4waddington'
    metadata_file = sys.argv[1]
    master_file = sys.argv[2]
    output_dir = sys.argv[3]
    genome = sys.argv[4]


    main(metadata_file,master_file,output_dir, genome=genome)
