# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 15:13:47 2023

@author: Anne

python3 /data/bagel5/scripts/extract_hmm_from_pfamdb.py -pfamdb /data/bagel5/db_resources/Pfam-A.hmm -outdir /data/bagel5/db_hmm -acc PF10417 

"""

import sys
import os
import argparse 
import re

parser = argparse.ArgumentParser(description='Routine to extract PFAM hmm from Pfam-A.hmm')
parser.add_argument('-pfamdb', dest='pfamdb', help='Pfam-A.hmm database', nargs='?', default='Pfam-A.hmm')
parser.add_argument('-acc', dest='acc', help='PFAM ACC e.g., PF00005')
parser.add_argument('-outdir', dest='outdir', help='output folder', nargs='?', default='.')

args, unknown = parser.parse_known_args()


#args.pfamdb='/data/bagel5/db_resources/Pfam-A.hmm'
#args.acc='PF00005'
#args.outdir='/data/bagel5/db_hmm'

with open(args.pfamdb, 'r') as pfamdb:
	export=False
	result='HMMER3/f [3.3 | Nov 2019]\n'
	for line in pfamdb:
		if re.search("ACC   "+args.acc, line):
			result+=prev_line
			export=True
		elif (export and re.search("//", line)):
			result+=line
			#print(result)
			with open(args.outdir+'/'+args.acc+'.hmm', 'w') as f: f.write(result)
			break
		if export: result+=line
		prev_line = line
