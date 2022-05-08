#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Hao, Chunlin @ 20210315

from utilities import *
import time, datetime, itertools, argparse, textwrap
import pandas as pd
from tqdm import tqdm
# ______________________________________________________________
parser = argparse.ArgumentParser(\
	formatter_class=argparse.RawDescriptionHelpFormatter,\
	description=textwrap.dedent('''\
		metaSpectraST (v0.0) by Hao, Chunlin.
		Compute spectral counts (SC) of cnsensus spectra.
								'''))
parser.add_argument('-s', nargs='?', action='store', dest='sptxt',\
					required=False, default='grandConsensus.sptxt',\
					help='consensus spectra .sptxt file, grandConsensus.sptxt by default.')
# parser.add_argument('-m', nargs='+', action='store', dest='mgf',\
# 					required=True, help='raw spectra data sets in MGF format') # For SIn only
argms = parser.parse_args()

sptxt_file = argms.sptxt

print('_________ metaSpectraST (v0.0) by Hao, Chunlin _________')
print(f'Start running metaSpectraST @ {datetime.datetime.now()}')
p0 = time.time()

# ______________________________________________________________
# Read sptxt file
sptxt = sptxt_parser(sptxt_file)
df_consensus_lib = pd.DataFrame(sptxt)
df_consensus_lib.to_csv('grandConsensus_lib.csv', index=False)

# ______________________________________________________________
sample_NumRep_matrix = []
FullName_RawSpectra_idx = []
for row in tqdm(df_consensus_lib.itertuples(), total=df_consensus_lib.shape[0]):
	fullname = getattr(row, 'FullName')
	spectrum = {'FullName': fullname}
	spectrum_idx = {'FullName': fullname}

	# Get num of consensus replicates in each sample
	for sample in getattr(row, 'NumRep'):
		spectrum[sample.split(',')[0]] = float(sample.split(',')[1])
		# below use complete num of replicates (NOT num chosen by spectraST)
		# spectrum[sample.split(',')[0]] = float(sample.split(',')[2])
	sample_NumRep_matrix.append(spectrum)
	
	# Get consensus idx of raw spectra
	for raw_spec in getattr(row, 'RawSpectra'):
		spectrum_idx['RawSpectra'] = raw_spec
		FullName_RawSpectra_idx.append(spectrum_idx)
		spectrum_idx = {'FullName': fullname}
	
df_NumRep = pd.DataFrame(sample_NumRep_matrix)
df_FullName_RawSpectra_idx = pd.DataFrame(FullName_RawSpectra_idx)
# print(df_FullName_RawSpectra_idx['FullName'].nunique())
# print(df_NumRep.shape)

df_NumRep.to_csv('unnorm_consensusPep_SC.csv', index=False)
df_FullName_RawSpectra_idx.to_csv('consensusSpec_RawSpectra_idx.csv', index=False)
print('\'unnorm_consensusPep_SC.csv\' CREATED.')
print('\'consensusSpec_RawSpectra_idx.csv\' CREATED.')
print(f'Job completed in {time.time() - p0} sec.')
