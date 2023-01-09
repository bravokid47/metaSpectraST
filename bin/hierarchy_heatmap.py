#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Hao, Chunlin @20210322
# Hao, Chunlin @20230109, fix the TypeError of variable 'observe', line 39, change str to int

import pandas as pd
import numpy as np
import re, argparse, textwrap

import matplotlib.pyplot as plt
import seaborn as sns

import fastcluster
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage


import sys
sys.setrecursionlimit(100000)


parser = argparse.ArgumentParser(\
	formatter_class=argparse.RawTextHelpFormatter,\
	description=textwrap.dedent('''\
		metaSpectraST (v0.0) by Hao, Chunlin.
		Hierarchical clustering based on consensus peptide/spectra SC or SIn.
								'''))
parser.add_argument('-n', nargs='?', action='store', dest='df',\
					required=True, help='Normalized SC or SIn data.')
parser.add_argument('-r', nargs='?', action='store', dest='host',\
					required=False, help='*.txt file containing a single column of FullName of consensus spectra that need to be excluded.')
parser.add_argument('-o', nargs='?', action='store', dest='observation',\
					required=False, default=1, help='''The minimum number of samples that a consensus spectrum has to be observed in.
Used to filter out singly or rarely observed cosensus spectra. Default is 1.
						''')
argms = parser.parse_args()
data = argms.df
host = argms.host
observe = int(argms.observation)

df_norm = pd.read_csv(data, index_col=0)
df_norm.index.name = 'Consensus peptides'

# Remove host consensus spectra
if host is not None:
	df_host = pd.read_csv(host)
	df_norm = df_norm[~(df.index.isin(df_host.iloc[:,0]))]

# Remove singly or rarely observed consensus spectra
df_norm = df_norm[(df_norm==0).sum(axis=1)<=(df_norm.shape[1]-observe)]

# Missing value imputation
df_norm = df_norm.apply(lambda x: np.log2(x+df_norm[df_norm>0].min().min()))

df_data = df_norm

# Dendrogram
Z = linkage(df_data.T, method='average')
plt.figure(figsize=(9.7,6))
plt.xlabel('Distance')
dendrogram(Z,\
    leaf_rotation=0,  # rotates the x axis labels
    leaf_font_size=8.,  # font size for the x axis labels
    color_threshold=260.,
    labels=df_data.columns.tolist(),
    orientation='left',
    #truncate_mode='lastp', p=4
)
plt.savefig('dendrogram.png', bbox_inches='tight', dpi=300)

cmap = sns.cubehelix_palette(rot=-.2, gamma=1.2, hue=1, as_cmap=True)
sns.set(font_scale=0.6)
g = sns.clustermap(df_data, figsize=(6, 9), method='average', metric='euclidean',
				   cbar_kws={'orientation':'vertical', 'ticklocation':'left',
				   			 'label':r'$\log_2 SI_N$'},
				   cbar_pos=(0.1, 0.82, 0.02, 0.16),
				   cmap=cmap, dendrogram_ratio=(0.13,0.2),
				   #col_colors=[sample_color, time_color],
				   xticklabels=True, yticklabels=False,
				   row_cluster=True, col_cluster=True,
)
				   
# Hide the dendrogram
g.ax_row_dendrogram.set_visible(False)
plt.savefig('hierarchicalHeatmap.png', bbox_inches='tight', dpi=300)

print('Dendrogram: dendrogram.png, CREATED.')
print('Hierarchical heatmap: hierarchicalHeatmap.png, CREATED.')
