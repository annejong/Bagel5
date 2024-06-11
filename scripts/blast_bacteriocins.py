# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 15:13:47 2023

@author: Anne

python3 /data/bagel5/scripts/blast_bacteriocins.py -query /data/bagel5/example/two

python3 /data/bagel5/scripts/blast_bacteriocins.py -query /data/bagel5/example/nisin
python3 /data/bagel5/scripts/blast_bacteriocins.py -query /data/bagel5/example/Thioalbamide
python3 /data/bagel5/scripts/blast_bacteriocins.py -query /data/www/users/vicky/bagel5_annotation/classI/Paenibacillus/NZ_CP024795.1/NZ_CP024795.1__AOI.all

blastp -outfmt 6 -max_target_seqs 1 -num_threads 6 -db /data/bagel5/db_proteins/bagel5_class12_bacteriocins_db.faa -query /data/www/users/vicky/bagel5_annotation/classI/Paenibacillus/NZ_CP024795.1/NZ_CP024795.1__AOI.all.faa -evalue 0.0001 -out /data/www/users/vicky/bagel5_annotation/classI/Paenibacillus/NZ_CP024795.1/NZ_CP024795.1__AOI.blastp


"""

import sys
import os
import argparse 
import re
import pandas as pd
import subprocess
import xml.etree.ElementTree as ET

parser = argparse.ArgumentParser(description='Blast bacteriocins and add Uniprot feature data')
parser.add_argument('-query', dest='query', help='Full path to fasta without .faa extension', nargs='?', default='.')
parser.add_argument('-db', dest='db', help='peptide database in folder db_proteins/ [bagel5_class12_bacteriocins_db.faa]', nargs='?', default='bagel5_class12_bacteriocins_db.faa')
parser.add_argument('-progdir', dest='progdir', help='Bagel5 program root [/data/bagel5]', nargs='?', default='/data/bagel5')
parser.add_argument('-cpu', type=str, dest='cpu', help='Number of threads [6]', nargs='?', default=6)
parser.add_argument('-evalue', type=str, dest='evalue', help='blastp evalue cutoff [1E-5]', nargs='?', default=0.00001)
args, unknown = parser.parse_known_args()

blastp_headers = [ 'query_id', 'subject_id', 'identity_percent', 'alignment_length','mismatches',
	'gap_opens', 'q_start', 'q_end', 's_start', 's_end','e_value', 'bit_score' ]

full_blastp_headers = ['query_id', 'subject_num', 'subject_id', 'subject_accession', 'subject_len',
	'bit_score', 'score', 'e_value', 'q_start', 'q_end', 's_start','s_end','mismatches',
	'identity_percent','positive','gap_opens','alignment_length','qseq','hseq','midline' ]

uniprot_features_headers = ['featureType', 'featureDescription', 'featurePos', 'featureStart', 'featureEnd']

LantiFeatures_headers = ['query_id','subject_id','Name','StructuralGene','Class','SubClass','SubSubClass'] +	uniprot_features_headers 



######################################################################################################
##                                     Functions                                                    ##
######################################################################################################


def run_cmd(cmd):
	result = ''
	result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
	#print("\tcommand stdout:", result.stdout)

def get_uniprot_features(uniprotID):
	# Read XML file
	with open(args.progdir+'/db_uniprot/'+str(uniprotID)+'.xml', 'r') as file: xml = file.read().replace('\n', '')
	index = 0
	features = []
	key = 0
	# Match feature elements in the XML
	while re.search(r'<feature (.*?)<\/feature>', xml[index:]):
		match = re.search(r'<feature (.*?)<\/feature>', xml[index:])
		feature = match.group(1)
		index += match.start() + 1
		key += 1
		# Extract feature attributes
		type_ = re.search(r'type=\"(.*?)\"', feature)
		type_ = type_.group(1) if type_ else ''
		description = re.search(r'description=\"(.*?)\"', feature)
		description = description.group(1) if description else ''
		beginpos = re.search(r'begin position=\"(.*?)\"', feature)
		beginpos = beginpos.group(1) if beginpos else ''
		endpos = re.search(r'end position=\"(.*?)\"', feature)
		endpos = endpos.group(1) if endpos else ''
		position = re.search(r'position position=\"(.*?)\"', feature)
		position = position.group(1) if position else ''
		# Append feature attributes to features list
		features.append((type_, description, position, beginpos, endpos))
		# Create DataFrame from features list
	return pd.DataFrame(features, columns=uniprot_features_headers)


def blast_results_to_df(xml_file):
    with open(xml_file, 'r') as file:
        xml_data = file.read()
    tree = ET.fromstring(xml_data)  # Parse the XML data
    iterations = tree.findall('.//Iteration')
    hits_data = []
    for iteration in iterations:
        query_id = iteration.find('Iteration_query-def').text
        query_len = iteration.find('Iteration_query-len').text
        hits = iteration.find('.//Iteration_hits')
        for hit in hits.findall('Hit'):
            subject_num = hit.find('Hit_num').text
            subject_id = hit.find('Hit_id').text
            subject_accession = hit.find('Hit_accession').text
            subject_len = hit.find('Hit_len').text
            hsps = hit.findall('.//Hsp')
            for hsp in hsps:
                bit_score = hsp.find('Hsp_bit-score').text
                score = hsp.find('Hsp_score').text
                e_value = hsp.find('Hsp_evalue').text
                q_start = hsp.find('Hsp_query-from').text
                q_end = hsp.find('Hsp_query-to').text
                s_start = hsp.find('Hsp_hit-from').text
                s_end = hsp.find('Hsp_hit-to').text
                identity = hsp.find('Hsp_identity').text
                positive = hsp.find('Hsp_positive').text
                gap_opens = hsp.find('Hsp_gaps').text
                alignment_length = hsp.find('Hsp_align-len').text
                qseq = hsp.find('Hsp_qseq').text
                hseq = hsp.find('Hsp_hseq').text
                midline = hsp.find('Hsp_midline').text
                mismatches = int(identity) - int(positive)
                identity_percent = int(100 * int(identity) / int(positive))

                hits_data.append([query_id, subject_num, subject_id, subject_accession, subject_len,
                                  bit_score, score, e_value, q_start, q_end, s_start, s_end, mismatches,
                                  identity_percent, positive, gap_opens, alignment_length, qseq, hseq, midline])

    return pd.DataFrame(hits_data, columns=full_blastp_headers)
	
	
def blast(basename):
	blastp_df=pd.DataFrame()
	blastp_df['db']=''
	blast_outputfile = basename+'.lanti.xml'
	cmd  = 'blastp -outfmt 5 -max_target_seqs 1 -num_threads '+str(args.cpu)
	cmd += ' -db '+args.progdir+'/db_proteins/'+args.db
	cmd += ' -query '+ basename+'.all.faa'
	cmd += ' -evalue '+str(args.evalue)
	cmd += ' -out '+ blast_outputfile
	run_cmd(cmd)
	return blast_results_to_df(blast_outputfile)	

def combine_Blast_BacAnn_UniProt(row):
    subject_id = row['subject_id']
    #print('subject_id===========>'+subject_id+'<<<<<<')  # Printing for debugging
    if subject_id in BacAnn_df['ID'].values:  # Check if subject_id exists in BacAnn_df['ID']
        BacAnn_row = BacAnn_df.loc[BacAnn_df['ID'] == subject_id]  # Get the row with matching ID
        uniprot_id = BacAnn_row['UniProt ID'].iloc[0]  # Get UniProt ID
        if str(uniprot_id) != 'nan':  # Check if UniProt ID is not NaN
            df = get_uniprot_features(uniprot_id)  # Fetch additional features from UniProt
            if df.empty:  # If no features found, create a placeholder DataFrame
                df = pd.DataFrame([['na','na','na','na','na']], columns=uniprot_features_headers)
            # Add additional columns to df from BacAnn_row
            df['query_id']       = row['query_id']
            df['subject_id']     = subject_id
            df['Name']           = BacAnn_row['Name'].iloc[0]
            df['StructuralGene'] = BacAnn_row['StructuralGene'].iloc[0]
            df['Class']          = BacAnn_row['Class'].iloc[0]
            df['SubClass']       = BacAnn_row['SubClass'].iloc[0]
            df['SubSubClass']    = BacAnn_row['SubSubClass'].iloc[0]
            df['Reference']      = BacAnn_row['Reference'].iloc[0]
            return df  # Return df with combined information
        else:
            return pd.DataFrame()  # Return an empty DataFrame if UniProt ID is NaN
    else:
        return pd.DataFrame()  # Return an empty DataFrame if subject_id is not found




######################################################################################################
##                                      Main                                                        ##
######################################################################################################

# load Bagel5 Bacteriocin annotation table
BacAnn_df = pd.read_table(args.progdir+'/tables/bagel5_bacteriocins_annotation.table', sep='\t')

# Perform blast
blast_df = blast(args.query)
blast_df.sort_values(by='query_id').to_csv(args.query+'.lanti.full.blastp', sep='\t', index=False)	

# Write as traditional -outfmt 6 blast table format
blast_df[blastp_headers].sort_values(by='query_id').to_csv(args.query+'.lanti.blastp', sep='\t', index=False)	


# Add UniProt lantibiotics features to blast subject(s)
LantiFeatures = []
for _, row in blast_df.iterrows():
    tmp_df = combine_Blast_BacAnn_UniProt(row)
    if not tmp_df.empty: LantiFeatures.append(tmp_df)

if LantiFeatures:
	LantiFeatures_df = pd.concat(LantiFeatures, ignore_index=True)	
	LantiFeatures_df[LantiFeatures_headers].sort_values(by='query_id').to_csv(args.query+'.lanti.features', sep='\t', index=False)	
	print('Output: '+args.query+'.lanti.features')


print('Output: '+args.query+'.lanti.full.blastp')
print('Output: '+args.query+'.lanti.blastp')


	