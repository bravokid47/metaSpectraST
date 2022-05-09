#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Hao, Chunlin @ 20200516

from utilities import *
import time, datetime, itertools, argparse, textwrap
import pandas as pd
from tqdm import tqdm
# ______________________________________________________________
parser = argparse.ArgumentParser(\
	formatter_class=argparse.RawDescriptionHelpFormatter,\
	description=textwrap.dedent('''\
		metaSpectraST (v0.0) by Hao, Chunlin.
		Compute normalized spectral index (SIn) of consensus spectra.
								'''))
parser.add_argument('-s', nargs='?', action='store', dest='sptxt',\
					required=False, default='grandConsensus.sptxt',\
					help='consensus spectra .sptxt file, grandConsensus.sptxt by default.')
parser.add_argument('-m', nargs='+', action='store', dest='mgf',\
					required=True, help='raw spectra data sets in MGF format') # For SIn only
argms = parser.parse_args()

sptxt_file = argms.sptxt
mgf_file = argms.mgf # list of mgf files

print('_________ metaSpectraST (v0.0) by Hao, Chunlin _________')
print(f'Start running metaSpectraST @ {datetime.datetime.now()}')
print()
print('Loading consensus spectra library ...')
p0 = time.time()


# ______________________________________________________________
# Read sptxt file
sptxt = sptxt_parser(sptxt_file)
df_consensus_lib = pd.DataFrame(sptxt)

# Remove singly observed spectra
# rm_singleton_spectra = df_consensus_lib['RawSpectra'].map(len) > 1
# df_consensus_lib = df_consensus_lib[rm_singleton_spectra]

# Remove charge state is 0
rm_zero_charge = df_consensus_lib['MW'] > 0
df_consensus_lib = df_consensus_lib[rm_zero_charge]

# ______________________________________________________________
# Read mgf files
print('Loading mgf files ...')
mgf_holder = {}
for entity in tqdm(mgf_file, total=len(mgf_file)):
	base_name = entity.split('/')[-1].split('.')[0].lower()
	#print(f'Reading {base_name}.mgf ...')
	mgf_holder[base_name] = mgf_reader(entity)

# ______________________________________________________________
# Calculate SIn
print()
print('Calculating spectral index ...')
sample_SIn_matrix = []
FullName_RawSpectra_idx = []
for row in tqdm(df_consensus_lib.itertuples(), total=df_consensus_lib.shape[0]):
	fullname = getattr(row, 'FullName')
	spectrum = {'FullName':fullname}
	spectrum_idx = {'FullName': fullname}
	len_spectra = float(getattr(row, 'MW'))/110
	sorted_raw_spectra = sorted(getattr(row, 'RawSpectra'))
	#print(f'\rCalculating consensus spectra ... [{fullname}]', end='')

	grouped_raw_spectra = []	
	for _, i in itertools.groupby(sorted_raw_spectra, lambda x: x.partition('.')[0]):
		grouped_raw_spectra.append(list(i))	
		
	for file_of_replicates in grouped_raw_spectra:
		intensity_of_all_replicates = []
		file_name = file_of_replicates[0].partition('.')[0] # prefix is sample name (mgf file name)
		for replicate in file_of_replicates:
			intensity_of_single_replicate = []
			try:
				for mgf_scan in mgf_holder[file_name]:
					if replicate in mgf_scan['TITLE']:
						for consensus_mz in getattr(row, 'MZ'):
							lower_idx, upper_idx = match_peak(mgf_scan['mz'], float(consensus_mz), 0.4)
							try:
								intensity_of_single_replicate.append(max(mgf_scan['intensity'][lower_idx:upper_idx]))
							except ValueError:
								pass # when there is no matched mz
			except KeyError:
				#pass
				print()
				print(f'Cannot find {file_name}.mgf, please check')
				exit(1)
			intensity_of_all_replicates.append(sum(intensity_of_single_replicate))
			
		spectrum[file_name] = sum(intensity_of_all_replicates)/len_spectra

	# Get consensus idx of raw spectra
	for raw_spec in getattr(row, 'RawSpectra'):
		spectrum_idx['RawSpectra'] = raw_spec
		FullName_RawSpectra_idx.append(spectrum_idx)
		spectrum_idx = {'FullName': fullname}

	sample_SIn_matrix.append(spectrum)
	
# ______________________________________________________________
df = pd.DataFrame(sample_SIn_matrix)
df.to_csv('unnorm_consensusPep_SI.csv', index=False)

df_FullName_RawSpectra_idx = pd.DataFrame(FullName_RawSpectra_idx)
df_FullName_RawSpectra_idx.to_csv('consensusSpec_RawSpectra_idx.csv', index=False)

df = df.set_index('FullName')
df = df.fillna(0)
# remove all 0 rows
df = df[(df.T!=0).any()]
# normalized by total SIn of each sample
df = df.div(df.sum(axis=0), axis=1)
df.to_csv('consensusPep_SIn.csv', index=True)

print()
print('\'unnorm_consensusPep_SI.csv\' CREATED.')
print('\'consensusPep_SIn.csv\' CREATED.')
print('\'consensusSpec_RawSpectra_idx.csv\' CREATED.')
print(f'Job completed in {time.time() - p0} sec.')

			


