# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 15:13:47 2023

@author: Anne

"""

import sys
import os
import argparse 
import pandas as pd
import urllib.request

parser = argparse.ArgumentParser(description='Routine to update bacteriocin uniprot descriptions as xml')
parser.add_argument('-programdir', dest='programdir', help='Bagel5 program root', nargs='?', default='/data/bagel5')
args, unknown = parser.parse_known_args()

df = pd.read_csv(args.programdir+'/tables/bagel5_bacteriocins_annotation.table',sep="\t", index_col=False, engine='python')

def downloadXML(row):
	xml = str(row['UniProt ID'])+'.xml'
	print('request '+xml)
	url = 'http://www.uniprot.org/uniprot/'+xml
	destination_file = args.programdir+'/db_uniprot/'+xml
	try:
		urllib.request.urlretrieve(url, destination_file)
	except:
		print(xml+' not found')

df.apply(downloadXML, axis=1)
