#!/usr/bin/env python3
# -*- coding:utf-8 -*-
'''
# File: variantSelection.py
# Created Date: Saturday July 4th 2020
# Author: Debora Antunes
# -----
# Last Modified: Friday, October 30th 2020, 12:52:47 pm
# -----
'''

import numpy as np
import pandas as pd
import gzip
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import VarianceThreshold

def loadData(name, case, control):
	"""Merges the data between the case and control datasets according to the variant, REF and ALT.
	Changes the format of the dataset to

	Samples		chrX:12568  chrX:12583  chrX:13587
	S1           0.0         0.0         0.0   
	S2           0.0         0.0         0.0   

	where rows are samples and columns are variants.

	Also creates a file INFO that provides the position, REF and ALT info for each variant.

	Arguments:
		name {string} -- Name of the chromosome.
		case {string} -- Path to the Cases dataset.
		control {string} -- Path to the Controls dataset.
	"""	
	print('>>> Loading datasets...')
	case = pd.read_csv(case, compression = 'gzip', sep='\t')
	var=list(case['VAR'])
	control = filterControls(control, var)
	dataset, info = mergeData(case,control)
	saveData(dataset, info, name)
	

def filterControls(control, var):
	"""Filters the variants of the controls dataset with the list of variants that exist in the cases dataset.

	Arguments:
		control {string} -- Path to controls dataset.
		var {list} -- List of variants in the cases dataset.

	Returns:
		pandas.Dataframe -- Filtered controls dataset.
	"""	
	df = pd.DataFrame()
	chunks=[]
	for chunk in pd.read_csv(control, quoting=3, sep='\t', chunksize=250000):
		chunks.append(chunk[chunk['VAR'].isin(var)].copy())
		data = pd.concat(chunks)
		df = pd.concat([df, data], ignore_index=True)
	return (df)

def mergeData(case,control):
	"""Merges the two dataset based in the position and ref and alt alleles. Only saves rows where these three \
		are the same.

	Arguments:
		case {pandas.Dataframe} -- Dataset for cases.
		control {pandas.Dataframe} -- Dataset for controls.

	Returns:
		tuple of pandas.Dataframe -- The merged dataset and the info for each variant saved.
	"""	
	print('>>> Merging...')
	dataset = pd.merge(case, control, on = ['VAR','REF','ALT'], how = 'inner')
	dataset = dataset.drop_duplicates()
	info = dataset[['VAR','REF','ALT']]
	dataset = dataset.drop(['REF','ALT'], axis=1).set_index('VAR')
	dataset.index.name = None
	dataset = dataset.T
	return (dataset, info)

def saveData(data, info, name):
	"""Saves the merged dataset in a csv.gz file. Also saves the info for each variant ina csv file.

	Arguments:
		data {pandas.Dataframe} -- Merged dataset
		name {string} -- The chromosome name
		info {pandas.Dataframe} -- The info dataset
	"""	
	df = data.to_csv(header=True, index=True, index_label='Samples')
	with gzip.open('../../data/datasets/chr/chr{}.csv.gz'.format(name), 'wt') as f:
		f.write(df)
		
	info.to_csv('../../data/datasets/chr/INFO_chr{}.csv'.format(name), header=True, index=False)
	
def mergeAll():
	"""Merges all chromosomes files in one, adds the label for cases and controls and saves a final merged csv.gz file.
	Changes the format of the dataset to

			chr1:12568  chr1:12583  chr1:13587  ... chrX:12568  chrX:12583  chrX:13587	labels
	0           0.0         0.0         0.0     ...     0.0         0.0         0.0       0
	1           0.0         0.0         0.0     ...     0.0         0.0         0.0       1

	where rows are samples and columns are variants, with the exception of the last one, labels, that
	characterizes the samples by case (1) and control (0).
	"""	
	region = [str(i) for i in range(1, 23)]
	region.extend(['X'])
	df = pd.DataFrame()
	for f in region:
		data = pd.read_csv('../../data/datasets/chr/chr{}.csv.gz'.format(f), compression = 'gzip')
		data = data.set_index('Samples')
		df = pd.concat([df, data], axis=1)

	for i in list(df.index):
		df.loc[i,'labels'] = '1' if i.startswith('Ex') else '0'
	df = df.reset_index().drop(['Samples'], axis=1)
	print(df)
	df.to_csv('../../data/datasets/merged_dataset.csv.gz', header=True, index=False, compression = 'gzip')

def cleanMissing(data, perc):
	"""Removes the columns that have more than a given percentage of missing data.

	Arguments:
		data {pandas.Dataframe} -- Dataset with no pre-processing.
		perc {float} -- Percentage of missing data chosen.

	Returns:
		pandas.Dataframe -- Cleaned dataset
	"""	
	cols = []
	for c in range(data.shape[1]-1):
		if data.iloc[:,c].isna().sum() / data.shape[0] * 100 > perc:
			cols.append(c)
	print('Number of deleted columns: ', len(cols))
	data = data.drop(data.columns[cols], axis=1)
	print(data)
	return data


def doImputation(data, strat):
	"""Imputation of NaN values using different stratagies

	Arguments:
		data {pandas.Dataframe} -- Dataset to impute
		strat {string} -- Imputation strategy

	Returns:
		pandas.Dataframe -- Imputated dataset
	"""	
	dataCol = np.array(data.columns)
	imp = SimpleImputer(missing_values=np.nan, strategy=strat)
	data=imp.fit_transform(data)
	data=pd.DataFrame(data, columns=dataCol)
	print(data)
	return data

