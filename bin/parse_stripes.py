#!/usr/bin/env python
import subprocess

thresh = 0.7
strip_file = '/home/areynolds/proj/bstripes/results/promoter_distal_distance_pairs_cis_trans/stripes.cis.mm10-198sample-dnaseI-pearsonr-092417.starch'
promoter_file= '/home/areynolds/proj/bstripes/data/mm10PromoterDHSs.092417.bed'
TSS_BED= '/home/areynolds/proj/bstripes/data/gencode.vM9.mouse.all_tss.bed'

promoter_dhss = set([tuple(x.strip().split()[0:3]) for x in iter(open(promoter_file).readline,'')])

cf = subprocess.Popen(['closest-features','--shortest',promoter_file,TSS_BED],stdout = subprocess.PIPE)
closest_promoter = {tuple(x.strip().split('|')[0].strip().split()):x.strip().split('|')[1].strip().split()[3] for x in iter(cf.stdout.readline,'')}
cf.stdout.close()

unstarch = subprocess.Popen(['unstarch',strip_file],stdout = subprocess.PIPE)
sort_bed = subprocess.Popen(['sort-bed','-'],stdin = subprocess.PIPE)
for line in iter(unstarch.stdout.readline,''):
    data = line.strip().split()
    if float(data[-1])<thresh: continue
    dtup = tuple(data[0:3])
    if tuple(data[0:3]) in closest_promoter:
        sort_bed.stdin.write('\t'.join(data[3:6] + [closest_promoter[dtup]] + data[0:3] + data[-1:])+'\n')
    else:
        dtup = tuple(data[3:6])
        sort_bed.stdin.write('\t'.join(data[0:3] + [closest_promoter[dtup]] + data[3:6] + data[-1:])+'\n')
unstarch.stdout.close()
sort_bed.stdin.close()
sort_bed.wait()
