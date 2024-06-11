"""
Anne de Jong
2023 Dec

python3 /data/bagel5/scripts/sORF.py -s /data/www/users/vicky/bagel5/classI -basename NZ_CP024795.1__AOI


"""

import sys
import pandas as pd
import argparse 
import subprocess
import os
import glob
import re
import numpy as np

parser = argparse.ArgumentParser(description='BAGEL5 small ORF calling script')
parser.add_argument('-s', dest='sessiondir', help='Session Folder', nargs='?', default='.')
parser.add_argument('-basename', dest='basename', help='basename without path, for .fna and .gff. And for result .sORF.gff and .sORF.fnn')
parser.add_argument('-min', dest='sORF_min', help='Minimum ORF size in bases', nargs='?', default=60)
parser.add_argument('-max', dest='sORF_max', help='Maximum ORF size in bases', nargs='?', default=300)
parser.add_argument('-overlap', dest='overlap', help='Overlap of ORFs in bases', nargs='?', default=2)
args, unknown = parser.parse_known_args()

complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
startcodons = 'ATG|GTG|TTG'
stopcodons  = 'TGA|TAG|TAA'
codons = 'ATG|GTG|TTG|GCC|AGT|TGT|CGA|ATC|AAC|AGC|TAC|TCG|ACA|CTG|CCG|GCA|AAG|GTT|CAC|AGA|ACC|CCA|TGG|CGC|CTC|CAG|ACG|AAA|GTA|CTT|GGA|GTC|TGC|TCA|ATT|TAT|AAT|ACT|CAA|GAC|GGT|TCC|TTT|AGG|CGT|CGG|CAT|ATA|CCC|GGG|GAG|TTA|CTA|GAT|TCT|TTC|GCG|GGC|GAA|GCT|CCT'
ORFpattern = re.compile(f"(.)({startcodons})((?:{codons}){{{args.sORF_min},{args.sORF_max}}})({stopcodons})")



gff_header  = ["chrom","db", "type", "start", "end", "name", "strand", "score", "description"]




def replace_with_X(dna, start, end):
    return dna[:start] + 'X' * (end - start) + dna[end:]

def mask(DNA, df): # Function to replace characters between start and end with 'X' and return the modified DNA
	masked_DNA = DNA
	for _, row in df.iterrows(): masked_DNA = replace_with_X(masked_DNA, row['start'], row['end'])
	return masked_DNA

def read_GFF(filename):
	GFF = pd.read_csv(filename, header=None,  comment='#',sep='\t', names=gff_header)
	convert_dict = { 'start': int, 'end': int }
	GFF = GFF.astype(convert_dict)  # be sure that start and end are integers
	GFF  = GFF.loc[GFF['type'] == 'gene'] # use sORF instead of gene
	return GFF


def read_fasta(filename):
	results = {}
	lines = open(filename, "r").read().split('\n')
	key=''
	seq=''
	for line in lines:
		if re.search("^>", line):
			if key != '': results[key] = seq
			my_list = re.match("^>(.*)", line)
			key = my_list.group(1)
			seq=''
		else:
		   seq += line
	if key != '': results[key] = seq  # add the last record
	return results

def inverse_complement(DNA):
	# add an N for unkown chars
	complemented_sequence = ''.join(complement_dict.get(base, 'N') for base in reversed(DNA))
	return complemented_sequence


def call_ORFs(DNA, strand, key):
	if strand == '-': DNA = inverse_complement(DNA) 
	index = 0
	count = 0
	df = pd.DataFrame(columns=gff_header)
	fnn = '' # fasta file for genes
	while match := ORFpattern.search(DNA[index:]):
		start = match.group(2)
		codons = match.group(3)
		stop = match.group(4)
		ORF = start+codons
		startPos = index + match.start()
		if strand == '-': startPos = len(DNA) - startPos - len(ORF) +1
		index += match.start() + args.overlap + len(ORF)  # allow some overlap between ORFs
		count += 1
		match = re.search(r'.*?(AOI.*)', key) # remove the chrom name from the locus_tag name
		locus_tag = match.group(1) + '_sORF_'+str(count).zfill(4)
		fnn+='>'+locus_tag+'\n'+ORF+'\n'
		df = df.append({"chrom":key,"db":'sORF', "type":'sORF', "start": startPos, "end": startPos+len(ORF), 
			"name": '.', "strand": strand, "score":'.', "description": 'locus_tag='+locus_tag+';seq='+ORF}, ignore_index=True)
	with open(args.sessiondir+'/'+args.basename+'.sORF.fnn', 'a') as file: file.write(fnn)
	return df

GFF = read_GFF(args.sessiondir+'/'+args.basename+'.gff')
FNA = read_fasta(args.sessiondir+'/'+args.basename+'.fna')
if os.path.exists(args.sessiondir+'/'+args.basename+'.sORF.fnn'):
	os.remove(args.sessiondir+'/'+args.basename+'.sORF.fnn') # remove old file

sORF_GFF = pd.DataFrame(columns=gff_header)
for key in FNA:
	DNA = mask(FNA[key], GFF[GFF['chrom'] == key])
	sORF_GFF = sORF_GFF.append(call_ORFs(DNA, '+', key), ignore_index=True)
	sORF_GFF = sORF_GFF.append(call_ORFs(DNA, '-', key), ignore_index=True)

#print(sORF_GFF.sort_values(by=['chrom','start']))
sORF_GFF.sort_values(by=['chrom','start']).to_csv(args.sessiondir+'/'+args.basename+'.sORF.gff', sep='\t', header=False, index=False)

	