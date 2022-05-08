#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Hao, Chunlin @ 20200516
import re, itertools

def sptxt_parser(sptxt_file):
	# Parse SpectraST generated consensus spectra library (.sptxt)
	consensus_lib = []
	single_consensus_spectrum = {'MZ':[]}
	
	with open(sptxt_file, 'r') as sptxt:
		for line in sptxt:
			if line[0].isdigit(): # read mz-intensity pairs and store ONLY mz
				single_consensus_spectrum['MZ'].append(re.split(r'[\t,\s]\s*', line)[0])	
			elif line.startswith('MW'): # read precursor molecular weight 
				single_consensus_spectrum['MW'] = float(re.split(r':\s*',line.strip())[1])
			elif line.startswith('FullName'): # full name of the consensus spectrum
				single_consensus_spectrum['FullName'] = re.split(r':\s*', line.strip())[1].lower()
			elif line.startswith('Comment'): # read comment and find records of raw spectra (replicates)
				get_raw_spectra = [s.lower() for s in line.split(' ') if s.startswith('RawSpectr')][0].split('=')[1].split(',')
				single_consensus_spectrum['RawSpectra'] = get_raw_spectra
				
				get_num_rep = [s.lower() for s in line.split(' ') if s.startswith('Sample=')][0].split('=')[1].split('/')[1:]
				single_consensus_spectrum['NumRep'] = get_num_rep

			elif line.startswith('NumPeaks'): # num of peaks of consensus spectrum
				single_consensus_spectrum['NumPeaks'] = re.split(r':\s*', line.strip())[1]
			# Append dict of spectrum (n-1) at the beginning of spectrum (n)	
			elif line.startswith('Name') and len(single_consensus_spectrum) == 6:
				consensus_lib.append(single_consensus_spectrum)
				single_consensus_spectrum = {'MZ':[]}
				
	consensus_lib.append(single_consensus_spectrum) # Append the last spectrum
	
	return consensus_lib


'''
# Skipped num of replicates 
def sptxt_parser(sptxt_file):
	# Parse SpectraST generated consensus spectra library (.sptxt)
	consensus_lib = []
	single_consensus_spectrum = {'MZ':[]}
	
	with open(sptxt_file, 'r') as sptxt:
		for line in sptxt:
			if line[0].isdigit(): # read mz-intensity pairs and store ONLY mz
				single_consensus_spectrum['MZ'].append(re.split(r'[\t,\s]\s*', line)[0])	
			elif line.startswith('MW'): # read precursor molecular weight 
				single_consensus_spectrum['MW'] = float(re.split(r':\s*',line.strip())[1])
			elif line.startswith('FullName'): # full name of the consensus spectrum
				single_consensus_spectrum['FullName'] = re.split(r':\s*', line.strip())[1].lower()
			elif line.startswith('Comment'): # read comment and find records of raw spectra (replicates)
				get_raw_spectra = [s.lower() for s in line.split(' ') if s.startswith('RawSpectr')][0].split('=')[1].split(',')
				single_consensus_spectrum['RawSpectra'] = get_raw_spectra
			elif line.startswith('NumPeaks'): # num of peaks of consensus spectrum
				single_consensus_spectrum['NumPeaks'] = re.split(r':\s*', line.strip())[1]
			# Append dict of spectrum (n-1) at the beginning of spectrum (n)	
			elif line.startswith('Name') and len(single_consensus_spectrum) == 5:
				consensus_lib.append(single_consensus_spectrum)
				single_consensus_spectrum = {'MZ':[]}
				
	consensus_lib.append(single_consensus_spectrum) # Append the last spectrum
	
	return consensus_lib
'''
	
def sptxt_get_replicate_num(sptxt_file):
	# Get num of replicates of each consensus spectrum from spectra library (.sptxt)
	consensus_lib = []
	single_consensus_spectrum = {}
	
	with open(sptxt_file, 'r') as sptxt:
		for line in sptxt:
			if line.startswith('MW'): # read precursor molecular weight 
				single_consensus_spectrum['MW'] = float(re.split(r':\s*',line.strip())[1])
			elif line.startswith('FullName'): # full name of the consensus spectrum
				single_consensus_spectrum['FullName'] = re.split(r':\s*', line.strip())[1].lower()
			elif line.startswith('Comment'): # read comment and find records of raw spectra (replicates)
				get_raw_spectra = [s.lower() for s in line.split(' ') if s.startswith('Sample=')][0].split('=')[1].split('/')[1:]
				single_consensus_spectrum['NumRep'] = get_raw_spectra
			elif line.startswith('NumPeaks'): # num of peaks of consensus spectrum
				single_consensus_spectrum['NumPeaks'] = re.split(r':\s*', line.strip())[1]
			# Append dict of spectrum (n-1) at the beginning of spectrum (n)	
			elif line.startswith('Name') and len(single_consensus_spectrum) == 4:
				consensus_lib.append(single_consensus_spectrum)
				single_consensus_spectrum = {}
				
	consensus_lib.append(single_consensus_spectrum) # Append the last spectrum
	
	return consensus_lib

		
def mgf_reader(mgf_file):
	mgf = []
	mgf_spectra = {'mz':[], 'intensity':[]}
	with open(mgf_file, 'r') as read_mgf:
		for line in read_mgf:
			if line[0].isdigit():
				line = re.split(r'[\t,\s]\s*', line.strip())
				mgf_spectra['mz'].append(float(line[0]))
				mgf_spectra['intensity'].append(float(line[1]))
				#mgf_spectra['mzIntensity'][line[0]] = float(line[1])
			elif line.startswith('TITLE'):
				mgf_spectra['TITLE'] = re.split(r'[\t,\s]\s*', line.strip().split('=')[1])[0].lower()
			elif line.strip() == 'END IONS':
				mgf.append(mgf_spectra)
				mgf_spectra = {'mz':[], 'intensity':[]}
				
	return mgf
		
def anchor_left(a_, x_):
    lo = 0
    hi = len(a_)

    while lo < hi:
        mid = (lo + hi) // 2
        if a_[mid] < x_:
            lo = mid + 1
        else:
            hi = mid
    return lo

def anchor_right(a_, x_):
    lo = 0
    hi = len(a_)

    while lo < hi:
        mid = (lo + hi) // 2
        if a_[mid] > x_:
            hi = mid
        else:
            lo = mid + 1
    return lo
		
def match_peak(mz_list, mz_value, tolerance=0.8):
	# find elements fall within the range of value+/-tolerance from array
	min_value = mz_value - tolerance 
	max_value = mz_value + tolerance 

	lower_idx = anchor_left(mz_list, min_value)
	upper_idx = anchor_right(mz_list, max_value)
    
	return lower_idx, upper_idx
	

		
		
		
		
		
		
		
		
		
		
		
		