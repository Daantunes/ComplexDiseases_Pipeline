#!/usr/bin/env python3
# -*- coding:utf-8 -*-
'''
# File: vcfParser.py
# Created Date: Tuesday June 30th 2020
# Author: Debora Antunes
# -----
# Last Modified: Thursday, April 15th 2021, 10:26:32 am
# -----
'''

import gzip
import sys
import codecs
import pickle
import pandas as pd
import sys
import os
import csv
from operator import itemgetter


def main(dataList, samplesFile=None):
	"""Loads the vcf file and transformes it in a csv.gz file that can be manipulated in pandas. \
		Also loads this csv.gz file in pandas and performes a filtration and transformation of data, \
		saving it in a new csv.gz file.

	Arguments:
		dataList {list} -- List of arguments to be used: the type of data (cases or control); the \
		chromosome name; the path to the dataset.

	Keyword Arguments:
		samplesFile {string} -- A path to a tsv file that provides the samples name to filter. In case\
		no file is provided, all samples are used.  (default: {None})

	Returns:
		[string] -- When the arguments provided are not correct, an error message is returned. 
	"""	
	# Checking that the first arg in valid
	name = dataList[0]
	if name not in ['controls', 'cases']: return 'ERROR1'

	#Checking number of args
	dataList.pop(0)
	if len(dataList) < 1: return 'ERROR2'

	#Checking that the second arg is valid
	chr = dataList[0]
	region = [str(i) for i in range(1, 23)]
	region.extend(['X'])
	if chr not in region: return 'ERROR3'

	#Checking number of args	
	dataList.pop(0)
	if len(dataList) < 1: return 'ERROR4'

	samples = listSample(samplesFile)
	print('>>> Loading dataset...')
	if not os.path.isfile('../../data/vcf/{}/output_{}.csv.gz'.format(name,chr)):
		readFile(samples, dataList[0], name, chr)

	header=True
	for chunk in pd.read_csv('../../data/vcf/{}/output_{}.csv.gz'.format(name,chr), 
			dtype={'CHROM': str, 'POS': str, 'REF': str, 'ALT': str, 'INFO': str}, 
			quoting=3,
			chunksize=250000,
			sep='\t'):
			
		if name == 'controls': chunk = filterControls(chunk)
		else: chunk = filterCases(chunk)
		chunk = translateGT(chunk)
		print(chunk)
		df = chunk.to_csv(header=header, sep='\t', quoting=3, index=False)
		with gzip.open('../../data/vcf/{}/outputPandas_{}.csv.gz'.format(name,chr), 'at') as f:
			f.write(df)
		header=False


def listSample(tsv):
	"""Opens the TSV file and filter the sample's names.

	Arguments:
		tsv {string} -- Path to TSV file.

	Returns:
		list -- List of samples to be used, or None in case none was provided.
	"""	
	if tsv: 
		print('>>> Reading Samples File...')
		samples = pd.read_csv(tsv, sep='\t')
		samples = list(samples.iloc[:,0])
	else: 
		print('>>> No Samples File provided!')
		samples = None
	return samples

def readFile(samples, path, name, region):
	"""Parses the VCF file and saves the right rows and columns in a csv.gz file.

	Arguments:
		samples {list} -- List of samples to be used, or None in case none was provided.
		path {string} -- Path to VCF file.
		name {string} -- Type of dataset (cases or controls).
		region {string} -- Name of chromosome.
	"""	
	reader = gzip.GzipFile(path,'r')
	if sys.version > '3':
		reader = codecs.getreader('ascii')(reader)
	line = next(reader)

	#Ignores the initial informative line
	while line.startswith('##'): line = next(reader)

	#Selects the columns to be saved
	col_names = line[1:(len(line)-1)].split('\t')
	if not samples:
		samples = col_names[9:]	
	col = ['CHROM','POS','REF','ALT','INFO']
	samples = [s for s in samples if s in col_names]
	col.extend(samples)
	index = [i for i, val in enumerate(col_names) if val in set(col)]

	#Writes header in csv.gz file
	csvfile = gzip.open('../../data/vcf/{}/output_{}.csv.gz'.format(name,region), 'at')
	writer = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL, delimiter='\t')
	writer.writerow(col)

	flag = 0
	write = sys.stdout.write
	flush = sys.stdout.flush
	line = next(reader)

	#Writes the remaining lines, filtering the columns of interest
	while line != '##' :
		data = list(itemgetter(*index)(line[0:(len(line)-1)].split()))
		writer.writerow(data)
		write('\r')
		write(str(flag))
		flush()
		line = next(reader, '##')
		flag+=1
	print('Done!')
	csvfile.close()

def filterControls(df):
	"""Filters the controls dataset to standardise the data.

	Arguments:
		df {pandas.Dataframe} -- Initial controls dataset.

	Returns:
		pandas.Dataframe -- Filtered controls dataset.
	"""	
	print('>>> Filtering dataframe...')
	df.insert(loc=0, column='VAR', value=df[['CHROM', 'POS']].apply(lambda x: ':'.join(x), axis=1))
	df['VAR']='chr'+df['VAR']
	col=df.columns
	for c in col:
		df[c] = df[c].str.replace('|','/')

	index=[]
	idx=list(df.index)
	types=['VT=SNP','VT=INDEL'] #filter by type of variant
	info=df.loc[:,'INFO'].str.split(';')
	for i in idx:
		for j in info[int(i)]:
			if j in types:
					index.append(i)

	df = df.loc[index,:]
	df = df.drop(['CHROM', 'POS', 'INFO'], axis=1)
	return df

def filterCases(df):
	"""Filters the cases dataset to standardise the data.

	Arguments:
		df {pandas.Dataframe} -- Initial cases dataset.

	Returns:
		pandas.Dataframe -- Filtered cases dataset.
	"""	
	print('>>> Filtering dataframe...')
	df.insert(loc=0, column='VAR', value=df[['CHROM', 'POS']].apply(lambda x: ':'.join(x), axis=1))
	df['VAR']='chr'+df['VAR']

	for i in range(df.shape[1]-6):
		df.iloc[:,i+6] = df.iloc[:,i+6].str[0:3] #saves only the genotype info
	df.iloc[:,6:].replace(to_replace = '.:.', value = '.', inplace = True)
	df = df.drop(['CHROM', 'POS', 'INFO'], axis=1)
	return df

def translateGT(dataset):
	"""Translates the genotype information to a number that a ML algorithm can understand.

	Arguments:
		dataset {pandas.Dataframe} -- Dataset to me translated.

	Returns:
		[pandas.Dataframe] -- Translated dataset.
	"""	
	print('>>> Starting translation...')

	gt1=['0/0','0/1','1/1','0/2','1/2','2/2','0/3',
		'1/3','2/3','3/3','0/4','1/4','2/4','3/4',
		'4/4','0/5','1/5','2/5','3/5','4/5','5/5',
		'0/6','1/6','2/6','3/6','4/6','5/6','6/6','.','./.']
	gt2=['0/0','1/0','1/1','2/0','2/1','2/2','3/0',
		'3/1','3/2','3/3','4/0','4/1','4/2','4/3',
		'4/4','5/0','5/1','5/2','5/3','5/4','5/5',
		'6/0','6/1','6/2','6/3','6/4','6/5','6/6','.','./.']
	tr=list(range(0,28))
	tr.append(float('nan'))
	tr.append(float('nan'))

	dataset.iloc[:,3:].replace(to_replace = gt1, value = tr, inplace = True)
	dataset.iloc[:,3:].replace(to_replace = gt2, value = tr, inplace = True)
	return dataset
