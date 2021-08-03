#!/usr/bin/env python
import sys,numpy as np,tempfile,os,subprocess,os
sys.path.append('/home/jlazar/sandbox/cancer')
import motifs

def main(peak_file,master_file,motif_file,output_file,motif_clusters = False):
    if motif_clusters != False: motif_clusters = True
    autosomal_peaks = get_overlaps(peak_file,master_file)
    motif_locs = motifs.motif_locations(autosomal_peaks,motif_file,motif_clusters = motif_clusters,
                                        by_motif = motif_clusters)
    lines = ['%s\t%f' %(tf,np.mean(motif_locs[tf])) for tf in motif_locs]
    with open(output_file,'w') as f1: f1.write('\n'.join(lines)+'\n')
    os.remove(autosomal_peaks)

def line_number(peaks):
    if peaks.endswith('.starch'):
        unzip = ['unstarch',peaks]
    else:
        unzip = ['zcat','-f', peaks]
    p1 = subprocess.Popen(unzip,stdout = subprocess.PIPE)
    count = 0
    for line in iter(p1.stdout.readline,''):
        count += 1
    p1.stdout.close()
    return count

def get_overlaps(peak_file,master_file):
    x,output_file = tempfile.mkstemp()
    if peak_file.endswith('.starch'):
        unzip = ['unstarch',peak_file]
    else:
        unzip = ['zcat','-f', peak_file]
    overlap = ['bedops','-e','25%',master_file,'-']
    tmp = tempfile.NamedTemporaryFile(delete=False)
    p1 = subprocess.Popen(unzip,stdout = subprocess.PIPE)
    p2 = subprocess.Popen(overlap,stdin = p1.stdout,stdout = subprocess.PIPE)
    for line in iter(p2.stdout.readline,''):
        length = line.strip().split()
        if len(length) < 4:
            length.append('1')
        tmp.write('\t'.join(length)+'\n')

    p1.stdout.close()
    p2.stdout.close()
    tmp.close()
    return tmp.name


if __name__ == "__main__":
    main(*sys.argv[1:])

