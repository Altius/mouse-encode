#!/usr/bin/env python
import HTSeq,subprocess

priority = ['promoter','exon','intron','utr','intergenic']

def main(bed_file,gencode_file):
	
	## create promoters
	## create exons
	## create utrs
	## create introns

	## map all those files

	map_name(bed_file,gencode_file + '.bed')


	pass


def create_tss(gencode_file,output_file):
	with open(output_file,'w') as f1:
		sb = subprocess.Popen(['sort-bed','-'],stdin = subprocess.PIPE,
			stdout = f1)
		for rg in iter(HTSeq.GFF_Reader(gencode_file)):
			if rg.type != 'transcript': continue

			name = rg.attr['gene_id']
			if 'gene_name' in rg.attr: name = rg.attr['gene_name']

			sb.stdin.write('%s\t%d\t%d\t%s\t%s\n' %
				(rg.iv.chrom,rg.iv.start_d,rg.iv.start_d+1,rg.iv.strand,
					name))
		sb.stdin.close()
		sb.communicate()


def create_promoter(gencode_file,output_file,size = 1000):
	with open(output_file,'w') as f1:
		sb = subprocess.Popen(['sort-bed','-'],stdin = subprocess.PIPE,
			stdout = f1)
		for rg in iter(HTSeq.GFF_Reader(gencode_file)):
			if rg.type != 'transcript': continue

			name = rg.attr['gene_id']
			if 'gene_name' in rg.attr: name = rg.attr['gene_name']
			if rg.iv.strand == '+':
				end = rg.iv.start
				start = end - size
			if rg.iv.strand == '-':
				start = rg.iv.start_d
				end = start + size

			sb.stdin.write('%s\t%d\t%d\t%s\t%s\n' %
				(rg.iv.chrom,start,end,rg.iv.strand,
					name))
		sb.stdin.close()
		sb.communicate()

def create_promoter(gencode_file,output_file,size = 1000):
	f2 = open('test.err','w')

	with open(output_file,'w') as f1:
		sb = subprocess.Popen(['sort-bed','-'],
			stdin = subprocess.PIPE,stdout = f1,stderr = f2)
		for rg in iter(HTSeq.GFF_Reader(gencode_file)):
			if rg.type != 'transcript': continue

			name = rg.attr['gene_id']
			if 'gene_name' in rg.attr: name = rg.attr['gene_name']
			if rg.iv.strand == '+':
				end = rg.iv.start
				start = max(0,end - size)
			else:
				start = rg.iv.start_d
				end = start + size
			if start < 0 or type(end) is not int:
				print start,end
			end = max(end,start + 1)
			sb.stdin.write('%s\t%d\t%d\t%s\t%s\n' %
				(rg.iv.chrom,start,end,rg.iv.strand,
					name))
		sb.stdin.close()
		sb.communicate()

def create_generic(gencode_file,output_file,atype = 'exon'):
	with open(output_file,'w') as f1:
		sb = subprocess.Popen(['sort-bed','-'],stdin = subprocess.PIPE,
			stdout = f1)
		for rg in iter(HTSeq.GFF_Reader(gencode_file)):
			if rg.type != atype: continue

			name = rg.attr['gene_id']
			if 'gene_name' in rg.attr: name = rg.attr['gene_name']

			sb.stdin.write('%s\t%d\t%d\t%s\t%s\n' %
				(rg.iv.chrom,rg.iv.start,rg.iv.end,rg.iv.strand,
					name))
		sb.stdin.close()
		sb.communicate()




