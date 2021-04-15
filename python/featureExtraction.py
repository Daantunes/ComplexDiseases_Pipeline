#!/usr/bin/env python3
# -*- coding:utf-8 -*-
'''
# File: featureSelection.py
# Created Date: Wednesday January 29th 2020
# Author: Debora Antunes
# -----
# Last Modified: Thursday, April 15th 2021, 12:09:55 pm
# -----
'''

import pandas as pd
import numpy as np
import pickle
import csv

import geneSelection
import models
import plots

from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from sklearn.decomposition import NMF
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.manifold import LocallyLinearEmbedding
from sklearn.manifold import TSNE


from progress.bar import Bar

def reduceFeatures(data, parameter='pval'):
	"""For \'pval\': Selects the averaged p value per gene minor than 0.05;
	For \'risk\': Selects risk genes from a list that also are present in the data
	For \'network\': Selects the top genes in the centrality measures
	
	Arguments:
		data {pandas.Dataframe} -- Dataset with genes
		parameter {string} -- performing different reductions according to the type of genes
	"""
	# Select genes by significance
	if parameter=='pval':
		try:
			data_simple = pickle.load(open('pickle/dataPval.p', 'rb'))

		except:
			print('>>> File not found, creating file')
			geneList = pd.read_csv('../../data/variants/sigVars.csv')
			data_simple = geneSelection.addGenes(data.copy(), geneList, 'pval') #Add p-values
			index = data_simple.iloc[:,:-1].T
			# print(data)
			print(index)
			# index = index[index[178]<1e-9]
			# print(list(index.index))

			#Select gene and pval rows and transform in a new dataset
			data_simple = data_simple.iloc[-2:,:].reset_index().rename(index={0: 'Genes', 1:'pval'}).T.iloc[1:-1,:]
			data_simple = data_simple.astype({'Genes':str,'pval':float})

			print('>>> Writting pickle...')
			pickle.dump(data_simple, open('pickle/dataPval.p', 'wb'))

		finally:
			
			#Calculates the average pvalue/ number of variants, per gene (if only average: aggfunc=np.mean)
			# data_simple = data_simple[data_simple['pval']>=1e-19]
			pvalues = data_simple.groupby('Genes')['pval'].apply(lambda x: x.mean())
			par = list(pvalues.index[pvalues<=0.05]) #alpha choosen to reduce features
	
	# Select genes by centrality
	elif parameter == 'network':
		par = pd.read_csv('../../data/proteins/r_genes.csv')
		netB = list(par.sort_values(by=['betweenness'], ascending = False).iloc[0:100,0])
		netC = list(par.sort_values(by=['closeness'], ascending = False).iloc[0:100,0])
		netD = list(par.sort_values(by=['degree_distribuition'], ascending = False).iloc[0:100,0])
		par = list((set(netB) & set(netC)) & set(netD))
		

	# Select genes by know risk
	riskGenes = pd.read_csv('../../data/genes/riskGenes.csv', sep=';')
	risk = list(riskGenes['Locus'])
	if parameter == 'risk':
		par = list(riskGenes['Locus'])
		risk=[]

	# Select genes in dataset
	genes = []
	geneGroup = data.iloc[:-1,:-1].groupby(data.iloc[-1,:-1], axis=1)
	gene = list(geneGroup.indices)
	if 'None' in gene:
		gene.remove('None')
	for i in gene:
		for j in set(par):
			if j == i:
				genes.append(i)
	
	# Save List of genes
	with open('../../data/genes/{}.csv'.format(parameter), mode='w') as f:
		f.write("\nTotal genes of interest: {}".format(len(par)))
		f.write("\nGenes of interest in dataset: {}".format(len(genes)))
		for g in genes:
			f.write('\n')
			f.write(g)
	
	# Extract features
	features = pd.DataFrame()
	for g in genes:
		features = extractFeatures(data, geneGroup, g, features)
	features['labels'] = data.iloc[:,-1]
	print(features)
	print('>>> Writing csv...')
	features.to_csv('../../data/datasets/reduced_dataset_{}.csv.gz'.format(parameter), index=False, compression = 'gzip')

def extractFeatures(data, geneGroup, gene, features):
	"""Creates new features (PDA, LDA, mean and variance) per gene
	
	Arguments:
		data {pandas.Dataframe} -- Dataset with genes
		geneGroup {pandas.DataFrame.groupby} -- Groups of genes composed with variants
		gene {string} -- Name of the gene
		features {pandas.Dataframe} -- Dataset with the reduced features
	
	Returns:
		pandas.Dataframe -- Dataset with the added reduced features
	"""	
	varList = list(geneGroup.groups[gene])
	for v in varList:
		data[v].iloc[:-1] = pd.to_numeric(data[v].iloc[:-1])
	pca = PCA(n_components=1)
	X_pca = pca.fit(data[varList].iloc[:-1,:]).transform(data[varList].iloc[:-1,:])
	X_pca = np.squeeze(np.asarray(X_pca))
	tsne = TSNE(n_components=1)
	X_tsne = tsne.fit_transform(data[varList].iloc[:-1,:])
	X_tsne = np.squeeze(np.asarray(X_tsne))
	mean = data[varList].iloc[:-1,:].mean(axis=1)
	var = data[varList].iloc[:-1,:].var(axis=1)
	var = var.fillna(0)
	features[str(gene) + '_PCA'] = X_pca
	features[str(gene) + '_TSNE'] = X_tsne
	features[str(gene) + '_mean'] = mean
	features[str(gene) + '_var'] = var
	return(features)
	
def selectTopFeatures(data, name, value):
	"""Using the repetion of tree models, selects the top genes to be used in the models.
	
	Arguments:
		data {pandas.Dataframe} -- Reduced dataset, with the new features
		name {string} -- Name for file
		value {int} -- Number of top features
	"""	
	try:
		sorted_top100 = pd.read_csv('../../data/features/top100_{}.csv'.format(name), header=None)
	except:
		print('>>> File not found, creating file')
		final_list=[]
		k=1000
		bar = Bar('Processing', max=1000,suffix='%(percent)d%%')
		#trains tree 1000 times
		
		while k > 0:
			X_train, X_test, y_train, y_test = train_test_split(data.iloc[:,:-1], data.iloc[:,-1], test_size = 0.4)
			result = models.useTree(data.iloc[:,:-1], X_train , y_train)
			final_list.extend(result)
			k-=1
			bar.next()
		bar.finish()

		final_list = np.array(final_list)
		feature, counts = np.unique(final_list, return_counts=True)
		top100 = dict(zip(feature, counts*100/1000))
		# ordered list of features and respective counts
		sorted_top100 = sorted(top100.items(), key=lambda kv: kv[1], reverse=True)

		
		print('>>> Writing csv...')
		with open('../../data/features/top100_{}.csv'.format(name), mode='w', newline='') as f:
				writer = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC)
				for l in sorted_top100:
					writer.writerow(l)
		sorted_top100 = pd.read_csv('../../data/features/top100_{}.csv'.format(name), header=None)
	finally:
		#selecting the best features
		top_names = list(sorted_top100.iloc[0:int(value),0])
		top_values = list(sorted_top100.iloc[0:int(value),1])

		# Save List of genes
		with open('../../data/genes/top100_{}.csv'.format(name), mode='w') as f:
			for g in top_names:
				f.write('\n')
				f.write(g)

		plots.plotFeatures(top_names,top_values, name)
		top_names.append('labels')
		topgenes = data[top_names]
		print('>>> Writing csv...')
		topgenes.to_csv('../../data/datasets/top_dataset_{}.csv.gz'.format(name), index=False, compression = 'gzip')