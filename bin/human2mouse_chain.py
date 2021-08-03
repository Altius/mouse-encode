#!/usr/bin/env python
import subprocess,os

mouse_file = '/home/cbreeze/bnmapper_analysis_new_masterlist/lifted.over.human.to.mouse.txt'
human_file = '/home/cbreeze/bnmapper_analysis_new_masterlist/mouse.to.human.alignment/lifted.over.mouse.to.human.txt'
outdir = '/home/jlazar/proj/mouse_encode/data/human2mouse_mapping'

def main(human_file,mouse_file,output_dir):
	if not os.path.exists(output_dir): os.mkdir(output_dir)

	human_file_name = os.path.split(human_file)[1]
	mouse_file_name = os.path.split(mouse_file)[1]
	sort_human = os.path.join(output_dir,
		human_file_name + '.sort.bed')
	sort_mouse = os.path.join(output_dir,
		mouse_file_name + '.sort.bed')

	with open(sort_human,'w') as f1:
		p1 = subprocess.Popen(['sort-bed',human_file],stdout = f1)
		p1.wait()
	with open(sort_mouse,'w') as f1:
		p1 = subprocess.Popen(['sort-bed',mouse_file],stdout = f1)
		p1.wait()


	human_output = os.path.join(output_dir,'human_reciprocal_dhs2.bed')
	mouse_output = os.path.join(output_dir,'mouse_reciprocal_dhs2.bed')

	find_reciprocal(sort_human,sort_mouse,human_output)
	find_reciprocal(sort_mouse,sort_human,mouse_output)

import tempfile
def get_reciprocal_peak(sort_mouse,sort_human,output_file,mouse_file = None,
				human_file = None):
	if mouse_file is not None:
		mtemp = tempfile.mktemp()
		subset_peaks(mouse_file,sort_mouse,mtemp)
	else:
		mtemp = sort_mouse
	if human_file is not None:
		htemp = tempfile.mktemp()
		subset_peaks(human_file,sort_human,htemp)
	else:
		htemp = sort_human

	find_reciprocal(mtemp,htemp,output_file)
	if mouse_file is not None: os.remove(mtemp)
	if human_file is not None: os.remove(htemp)	

def subset_peaks(pk_file,sort_file,output_file):
	with open(output_file,'w') as f1:
		p1 = subprocess.Popen(['bedmap','--echo','--skip-unmapped','--fraction-either','0.25',sort_file,
						pk_file],stdout = f1)
		p1.wait()



def find_reciprocal(speciesa,speciesb,output_file):
	
	#f1 = open(output_file + '.atmp','w')
	mapped_region = subprocess.Popen('''cut -f1-4 %s''' %speciesa,
			stdout = subprocess.PIPE,shell = True)
	cmd = '''awk -F '["\t":-]' '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3}' '''
	parsed_region = subprocess.Popen([cmd],
			stdout = subprocess.PIPE,stdin = mapped_region.stdout,shell = True)	
	sort_r = subprocess.Popen(['sort-bed','-'],stdin = parsed_region.stdout,
				stdout = subprocess.PIPE)

	reciprocal = subprocess.Popen(['bedmap','--echo','--indicator','--echo-map','-',speciesb],
		stdin = sort_r.stdout,stdout = subprocess.PIPE)


	inds = []
	for line in iter(reciprocal.stdout.readline,''):
		og,ind,maps = line.strip().split('|')
		if int(ind) == 0:
			inds.append(0)
		else:
			inds.append(overlap(og.split()[3:],maps))

	mapped_region.stdout.close()
	parsed_region.stdout.close()
	sort_r.stdout.close()
	reciprocal.stdout.close()

	fold = open(speciesa)
	with open(output_file,'w') as f1:
		count = 0
		for line in fold:
			if inds[count]: f1.write(line)
			count += 1
	fold.close()

def overlap(original_values,maps):
	anyv = 0
	for dt in maps.split(';'):
		mapped_loc = dt.split()[3]
		mchr,mbp = mapped_loc.split(':')
		start,end = mbp.split('-')

		m_start = int(start)
		m_end = int(end)

		o_start = int(original_values[1])
		o_end = int(original_values[2])

		if mchr != original_values[0]: continue
		if o_end < m_start or o_start > m_end: continue

		anyv += 1
	return anyv > 0


def find_reciprocal_old(speciesa,speciesb,output_file):
	
	mapped_region = subprocess.Popen('''cut -f4 %s''' %speciesa,
			stdout = subprocess.PIPE,shell = True)
	cmd = '''awk -F '[:-]' '{print $1"\t"$2"\t"$3}' '''
	parsed_region = subprocess.Popen([cmd],
			stdout = subprocess.PIPE,stdin = mapped_region.stdout,shell = True)

	
	sort_r = subprocess.Popen(['sort-bed','-'],stdin = parsed_region.stdout,
				stdout = subprocess.PIPE)
	reciprocal = subprocess.Popen(['bedmap','--indicator','-',speciesb],
		stdin = sort_r.stdout,stdout = subprocess.PIPE)

	inds = []
	for line in iter(reciprocal.stdout.readline,''):
		inds.append(int(line.strip()))
	mapped_region.stdout.close()
	parsed_region.stdout.close()
	sort_r.stdout.close()
	reciprocal.stdout.close()

	fold = open(speciesa)
	with open(output_file,'w') as f1:
		count = 0
		for line in fold:
			if inds[count]: f1.write(line)
			count += 1
	fold.close()



if __name__ == "__main__":
	main(human_file,mouse_file,outdir)