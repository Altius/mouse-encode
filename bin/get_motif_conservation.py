#!/usr/bin/env python
import motifs,tempfile,subprocess,os,sys

def main(pk_file,tf,motif_file,conservation_file):

    if not os.path.getsize(pk_file): return 0,0

    motif2tf = motifs.get_tfs_from_file(motifs.MOTIF_MAPPINGS_FILE,False)
    tfile2 = tempfile.mktemp()

    f1 = open(tfile2,'w')
    pm = subprocess.Popen(['bedops','-m','-'],stdin = subprocess.PIPE,stdout = f1)
    p1 = subprocess.Popen(['tabix','-R',pk_file,motif_file],stdout = subprocess.PIPE)
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
    p1 = subprocess.Popen(['tabix','-R',tfile2,conservation_file],stdout = subprocess.PIPE)
    for line in iter(p1.stdout.readline,''):
        mt = line.strip().split()
        try:
            cons = float(mt[4])
            vls += cons
        except:
            continue
        count += 1
    p1.stdout.close()

    os.remove(tfile2)
    if count == 0:
        return 0,count
    else:
        return vls/float(count),count

if __name__ == "__main__":
    c,n = main(*sys.argv[1:])
    sys.stdout.write('%f\t%d\n' %(c,n))