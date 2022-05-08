#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Hao, Chunlin @ 20210315

import pandas as pd
import numpy as np
import difflib, argparse, textwrap
import rpy2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

def divided_by_corrected_lib_size(scalingFactor, data_col):
	# scalingFactor, edgeR DataFrame, having column of lib.size and norm.factors.
	# data_col, pd.Series, SC or SIc values of different consensus spectra of a particular sample.
	sample = data_col.name
	return data_col/(df_scalingFactor.loc[sample, 'lib.size']*df_scalingFactor.loc[sample, 'norm.factors'])

# ______________________________________________________________________________________
parser = argparse.ArgumentParser(\
	formatter_class=argparse.RawDescriptionHelpFormatter,\
	description=textwrap.dedent('''\
		metaSpectraST (v0.0) by Hao, Chunlin.
		Normalize the SC/SIn data with TMM method (trimmed mean of M values).
								'''))
parser.add_argument('-u', nargs='?', action='store', dest='df',\
					required=True, help='Unnormalized SC or SIn data.')
argms = parser.parse_args()
data = argms.df
df_unnorm = pd.read_csv(data, index_col=0)
df_unnorm = df_unnorm.fillna(0)

# Calculate TMM scaling factor using R package edgeR
edgeR = importr('edgeR') # edgeR_3.34.0, limma_3.48.1
tmm_scalingFactor = ro.r('''
	function(unnorm_data){
	unnorm <- read.csv(file=unnorm_data)
	unnorm[is.na(unnorm)] <- 0
	unnorm <- unnorm[-c(1)] # remove column of FullName of consensus spectra
	matrix_data <- data.matrix(unnorm)
	y <- DGEList(counts=matrix_data, remove.zeros=TRUE)
	y <- calcNormFactors(y)
	y$samples
	write.csv(y$samples, file='tmm_scalingFactor.csv')
	}''')

tmm_scalingFactor(data)


df_scalingFactor = pd.read_csv('tmm_scalingFactor.csv', index_col=0)

# edgR sometimes may change the sample name
correct_sampleName = df_unnorm.columns
correct_edgR_index = {}
for sample in df_scalingFactor.index:
	correct_edgR_index[sample] = difflib.get_close_matches(sample, correct_sampleName, n=1)[0]

df_scalingFactor.rename(index=correct_edgR_index, inplace=True)

df_norm = df_unnorm.apply(lambda x: divided_by_corrected_lib_size(df_scalingFactor, x), axis=0)
df_norm.to_csv('tmmNorm_consensusPep.csv')
print('\'tmmNorm_consensusPep.csv\' CREATED.')
