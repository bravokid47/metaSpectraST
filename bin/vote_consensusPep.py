#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Hao, Chunlin @ 20211115
import glob, argparse, textwrap
import pandas as pd
import numpy as np
from collections import Counter
from tqdm import tqdm

parser = argparse.ArgumentParser(\
	formatter_class=argparse.RawTextHelpFormatter,\
	description=textwrap.dedent('''\
		metaSpectraST (v0.0) by Hao, Chunlin.
		Reconcile any conflicting peptide sequences assigned to replicate spectra and their consensus spectrum by database search, open modification search, or de novo sequencing.
		Peptide sequences of every spectrum in one spectral cluster vote for the consensuspeptide sequences by majority rule.
								'''))
parser.add_argument('-d', nargs='?', action='store', dest='_df_db_raw',\
					required=True, help='''Result of database search of raw spectra.
Should be organized as a CSV file containing columns of \'RawSpectrum\', \'peptide\', and \'protein\'.''')

parser.add_argument('-t', nargs='?', action='store', dest='_df_tag_raw',\
					required=True, help='''Result of open modification of raw spectra.
Should be organized as a CSV file containing columns of \'RawSpectrum\', \'peptide\', and \'protein\'. Column of \'protein\' can be empty.''')

parser.add_argument('-n', nargs='?', action='store', dest='_df_denovo_raw',\
					required=True, help='''Result of de novo sequencing of raw spectra.
Should be organized as a CSV file containing columns of \'RawSpectrum\', \'peptide\', and \'protein\'. Column of \'protein\' can be empty.''')

parser.add_argument('-c', nargs='?', action='store', dest='_df_db_consensus',\
					required=True, help='''Result of database search of consensus spectra.
Should be organized as a CSV file containing columns of \'consensusSpec\', \'peptide\', and \'protein\'. Column of \'protein\' can be empty.''')

parser.add_argument('-i', nargs='?', action='store', dest='_df_idx',\
					required=False, default='consensusSpec_RawSpectra_idx.csv',\
					help='''Index of correspondence of raw spectrum and its consensus spectrum.
By default, \'consensusSpec_RawSpectra_idx.csv\' generated automatically by metaSpectraST_SC.py/metaSpectraST_SIn.py will be used.''')
argms = parser.parse_args()


try:
	df_db_raw = pd.read_csv(argms._df_db_raw, usecols=['RawSpectrum', 'peptide', 'protein'])
except ValueError as msg1:
	print(f'Files of database search of raw spectra: {msg1}')

try:
	df_taggraph_raw = pd.read_csv(argms._df_tag_raw, usecols=['RawSpectrum', 'peptide', 'protein'])
except ValueError:
	try:
		df_taggraph_raw = pd.read_csv(argms._df_tag_raw, usecols=['RawSpectrum', 'peptide'])
	except ValueError as msg2:
		print(f'Files of open modification search of raw spectra: {msg2}')
		
try:
	df_denovo_raw = pd.read_csv(argms._df_denovo_raw, usecols=['RawSpectrum', 'peptide', 'protein'])
except ValueError:
	try:
		df_denovo_raw = pd.read_csv(argms._df_denovo_raw, usecols=['RawSpectrum', 'peptide'])
	except ValueError as msg3:
		print(f'Files of de novo sequencing of raw spectra: {msg3}')

try:
	df_db_cons = pd.read_csv(argms._df_db_consensus, usecols=['consensusSpec', 'peptide', 'protein'])
except ValueError as msg4:
	print(print(f'Files of database search of consensus spectra: {msg4}'))

df_consID_rawID = pd.read_csv(argms._df_idx)

cons_list = df_consID_rawID['FullName'].unique()
df_voted_cons = []
for i in tqdm(cons_list, total=len(cons_list)):
	df_vote = []
	db_cons = df_db_cons[df_db_cons['consensusSpec'] == i].copy() # '_0317p_v_15661/2', _0309a_s_05380/3
	if db_cons.empty == False:
		db_cons.rename(columns={'consensusSpec':'RawSpectrum'}, inplace=True)
		for x in range(db_cons.shape[0]): # check _0317p_s_19259/3; df_db_cons has two rows of _0317p_s_19259/3, with different assumed charge.
			db_cons_single = db_cons.iloc[x, :].to_frame().T
			db_cons_single.insert(0, column='source', value = ['1_db'])
			df_vote.append(db_cons_single)

	cons_chunk = df_consID_rawID[df_consID_rawID['FullName'] == i].copy()
	raw_spectra_list = cons_chunk['RawSpectra'].to_list()
	
	for raw_spec in raw_spectra_list:
		db_raw = df_db_raw[df_db_raw['RawSpectrum'] == raw_spec]
		tag_raw = df_taggraph_raw[df_taggraph_raw['RawSpectrum'] == raw_spec]
		denovo_raw = df_denovo_raw[df_denovo_raw['RawSpectrum'] == raw_spec]

		if db_raw.empty:
			if tag_raw.empty:
				if denovo_raw.empty:
					continue
				else:
					denovo_raw.insert(0, column='source', value = ['3_denovo'])
					df_vote.append(denovo_raw)
			else:
				tag_raw.insert(0, column='source', value = ['2_open'])
				df_vote.append(tag_raw)
		else:
			db_raw.insert(0, column='source', value = ['1_db'])
			

	try:
		df_vote = pd.concat(df_vote)
	
		vote = Counter(df_vote['peptide'])
		winner_count = max(vote.values())
		winner_seq = [i for i in vote.keys() if vote[i] == winner_count][0]
		df_final = df_vote[df_vote['peptide']==winner_seq].sort_values(by=['source'])
		final = df_final.iloc[0, :].copy().rename(i, inplace=True)
		df_voted_cons.append(final)
	except ValueError:
		df_voted_cons.append(pd.Series({'source':None, 'RawSpectrum':None, 'peptide':None, 'protein':None}, name=i))
		
df = pd.concat(df_voted_cons, axis=1).T
df.to_csv('final_voted_consensus.csv', index_label='FullName')
print('\'final_voted_consensus.csv\' CREATED.')
print('Job completed!')

