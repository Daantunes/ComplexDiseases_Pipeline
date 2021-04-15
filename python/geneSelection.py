#!/usr/bin/env python3
# -*- coding:utf-8 -*-
'''
# File: geneSelection.py
# Created Date: Monday February 10th 2020
# Author: Debora Antunes
# -----
# Last Modified: Thursday, April 15th 2021, 10:37:21 am
# -----
'''
import pandas as pd 
import numpy as np 
import pickle 
import csv 

import plots 

import pyensembl

from rpy2.robjects.packages import importr
import rpy2.robjects as ro 

from scipy import stats

from progress.bar import Bar


def normalityTest(data):
	"""Test the normality of each variant; creates a plot
	
	Arguments:
		data {pandas.Dataframe} -- Dataset to test normality
	"""	
	print('>>> Testing Normality of Variants...')
	notNormal = 0
	Normal = 0
	plotData = {'chrNo': [], 'posNo': [], 'val': []} #Creates a dataset for the plot

	bar = Bar('Processing', max=len(data.columns[:-1]),suffix='%(percent)d%%')
	notnormal=[[],[]]
	for c in data.columns[:-1]:
		k2, p = stats.normaltest(data[c]) #P-value of the test
		
		plotData, hA = createPlotData(plotData, c, p)
		if hA: Normal += 1
		else: 
			notNormal +=1
			notnormal[0].append(c)
			notnormal[1].append(p)
		bar.next()
	bar.finish()

	plotData = pd.DataFrame(data=plotData)
	print('>>> Writting pickle...')
	pickle.dump(plotData, open('pickle/plotData_normality.p', 'wb'))
	plots.plotManhattan(plotData, '_normality')
	print('Normal: ', Normal)
	print('Not Normal: ', notNormal)
	print('>>> Writing txt...')
	result=zip(notnormal[0],notnormal[1])
	with open('../../data/variants/notnormal.csv', mode='w', newline='') as f:
		wr = csv.writer(f)
		wr.writerow(("Variant", "pval"))
		for item in result:
			wr.writerow(item)


def chiSquaredTest(data):
	"""Test if there is no association with the disease; creates a plot
	
	Arguments:
		data {pandas.Dataframe} -- Dataset to test
	"""	
	print('>>> Testing Significance of Variants...')
	notSig = 0
	Sig = 0
	plotData = {'chrNo': [], 'posNo': [], 'val': []} #Creates a dataset for the plot
	listSig=[[],[]]

	bar = Bar('Processing', max=len(data.columns[:-1]),suffix='%(percent)d%%')	

	for c in data.columns[:-1]:
		table = pd.crosstab(data[c], data['labels'])
		f_obs = np.array(table.values)
		chi2, p = stats.chi2_contingency(f_obs, correction=False)[0:2]

		listSig[0].append(c)
		listSig[1].append(p)
		plotData, hA = createPlotData(plotData, c, p)
		if hA == True:
			Sig += 1
		else:
			notSig +=1
		bar.next()
	bar.finish()

	plotData = pd.DataFrame(data=plotData)
	print('>>> Writting pickle...')
	pickle.dump(plotData, open('pickle/plotData_chi2.p', 'wb'))
	plots.plotManhattan(plotData, '_chi2')
	print('Significant: ', Sig)
	print('Not Significant: ', notSig)
	print('>>> Writing txt...')
	result=zip(listSig[1],listSig[0])
	with open('../../data/variants/sigVars.csv', mode='w', newline='') as f:
		wr = csv.writer(f)
		wr.writerow(("pval", "Variants"))
		for item in result:
			wr.writerow(item)
	
def createPlotData(plotData, var, p):
	"""Prepares data to use in the Manhattan plot

	Args:
		plotData (dict): Dictionary with the information
		var (string): variant
		p (float): p-value of the variant

	Returns:
		tuple: The updated dictionary, and the result for the alternative hypothesis
	"""	
	alpha = 0.05
	if p <= 0:
		p = 1e-40
		hA = True
	elif p < alpha:
		hA = True
	elif p >= alpha:
		hA = False

	chrNo, posNo = var.split(':') #To populate dataset
	chrNoRes = makeInteger(chrNo[3:])
	if not chrNoRes:
		if chrNo[3:] == 'X':
			plotData['chrNo'].append(23)
			plotData['posNo'].append(int(float(posNo)))
			plotData['val'].append(p)
		elif chrNo[3:] == 'Y':
			plotData['chrNo'].append(24)
			plotData['posNo'].append(int(float(posNo)))
			plotData['val'].append(p)
	else: 
		plotData['chrNo'].append(chrNoRes)
		plotData['posNo'].append(int(float(posNo)))
		plotData['val'].append(p)

	return plotData, hA

def makeInteger(num):
	"""Tries to change the type of value to integer
	
	Arguments:
		num {value} -- value to be changed
	
	Returns:
		integer -- value changed
		boolean -- warning that value type has not changed
	"""	
	try: return int(num)
	except: return False
	else: return False

def geneTranslation(var):
	"""Translates variants to genes
	
	Arguments:
		var {list} -- List of chromosome regions
	"""	
	data = pyensembl.Genome(reference_name='GRCh37',
				annotation_name='my_genome_features', 
				gtf_path_or_url='../../data/datasets/Original/ensembl_v37.gtf')
	data.index()
	genes, variants = [], []
	bar = Bar('Processing', max=len(var),suffix='%(percent)d%%')
	for v in var:
		locus = v.split(':')
		chr, pos = locus[0][3:], int(float(locus[1]))
		gene = data.genes_at_locus(chr,pos)
		if len(gene) == 0:
			genes.append('None')
			variants.append(v)
		elif gene[0].biotype == 'protein_coding':
			genes.append(gene[0].gene_name)
			variants.append(v)
		else:
			genes.append('None')
			variants.append(v)
		bar.next()
	bar.finish()
	result=zip(genes,variants) # Creates a list of (gene,variant region)
	print('>>> Writing txt...')
	with open('../../data/genes/geneList.csv', mode='w', newline='') as myfile:
		wr = csv.writer(myfile)
		wr.writerow(("Genes", "Variants"))
		for row in result:
			wr.writerow(list(row))


def addGenes(data, genes, col):
	"""[Creates a row with corresponding values]
	
	Arguments:
		data {pandas.Dataframe} -- Dataset with variants in columns
		genes {pandas.Dataframe} -- Dataset with variants and corresponding values
		col {string} -- Name of column to add

	Returns:
		[pandas.Dataframe] -- Dataset the added corresponding values
	"""	
	data.loc[data.shape[0]] = 'NONE' #Creates a last line with the string 'NONE'

	bar = Bar('Processing', max=len(genes.index))


	for i in range(len(genes.index)):
		g = genes.iloc[i,1]
		data.loc[data.shape[0]-1, g] = genes.iloc[i,0]
		bar.next()
	bar.finish()
	print(data)
	return data
	



