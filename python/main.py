#!/usr/bin/env python3
# -*- coding:utf-8 -*-
'''
# File: main.py
# Created Date: Wednesday March 25th 2020
# Author: Debora Antunes
# -----
# Last Modified: Thursday, April 15th 2021, 11:12:53 am
# -----
'''

import pandas as pd
import argparse
import pickle

import vcfParser
import variantSelection
import geneSelection
import featureExtraction
import models
import plots
import network



def main():
	arg_parser = argparse.ArgumentParser(description = 'Feature Selection')
	arg_parser.add_argument('-v','--vcf', type = str, nargs=3, 
				help = 'First argument must be \'cases\' or \'controls\'; Second argument must be the \
					number of the chromosome (1 to 22) or \'X\'; \
					Third argument is the path to the VCF file to be processed')
	arg_parser.add_argument('-s','--samples', type = str, 
				help = 'TSV samples file from 1000 Genome')			
	arg_parser.add_argument('-m','--merge', type = str, nargs='*',
				help = 'Merge datasets from a specific chromosome (1 to 22 and X) of cases and controls \
					for the construction of dataset; First arg is a the chromosome, Second is the cases \
					dataset, and Third is the controls dataset. If the argument \'all\' is provided, merges all \
					choromosome datasets in one but all need to be available.')
	arg_parser.add_argument('-w','--workspace', type = str, 
				help = 'Saved workspace')
	arg_parser.add_argument('-d','--dataset', type = str, 
				help = 'Input file')
	arg_parser.add_argument('-cV','--cleanVar', type = float, 
				help = 'Selects the columns to delete according to a percentage')
	arg_parser.add_argument('-iV','--imputeVar', type = str, 
				help = 'Imputation strategy: \'mean\', \'median\',\'most_frequent\'')
	arg_parser.add_argument('-c','--classifier', type = str, nargs=1, 
				help = 'Select one model: \'svm\',\'tree\',\'knn\',\'log\',\'rf\',\'nb\'')
	arg_parser.add_argument('-g','--grid_search', type = str, nargs='*', 
				help = 'Chosen hyperparameters in format \'parameter:value\'')
	arg_parser.add_argument('-t','--test', type = str, 
				help = 'Apply test per variant: \'normality\', \'chi2\'')
	arg_parser.add_argument('-tG','--translation_genes', type = str, 
				help = 'Path to list of variants to translate into genes')
	arg_parser.add_argument('-aG','--add_genes', type = str, 
				help = 'Path to list of genes and variants')
	arg_parser.add_argument('-nF','--new_features', type = str, 
				help = 'Group by gene according to the given parameter (pval, risk or network) and apply metric (pca, tsme, mean and variance)')
	arg_parser.add_argument('-tF','--top_features', type = float, 
				help = 'Selects the top <value> features')
	arg_parser.add_argument('-n','--network', type = str, nargs='*', 
				help = 'Performs the integration into the network: the first argument is the path to the interactions file;\
					the second argument is the path to a Rdata file with the network saved')
	args = arg_parser.parse_args()

	if args.vcf:
		if args.samples: 
			flag = vcfParser.main(args.vcf, args.samples)
		else: 
			flag = vcfParser.main(args.vcf)
		if flag == 'ERROR1':
			print('>>> First argument is not \'cases\' or \'controls\', please add a correct label.')
		elif flag =='ERROR2':
			print('>>> No mode provided.')
		elif flag =='ERROR2':
			print('>>> Mode not valid, make sure is \'region\' or one chromossome (1 to 22, or X)')
		elif flag == 'ERROR4':
			print('>>> No VCF file provided.')

	elif args.merge:
		if args.merge[0] == 'all':
			variantSelection.mergeAll()
		else:
			variantSelection.loadData(args.merge[0],args.merge[1],args.merge[2])

	elif args.workspace:
		print('>>> Loading workspace...')
		dataset = pickle.load(open(args.workspace, 'rb'))
		print(dataset)
		name = args.workspace[7:]

	elif args.dataset:
		print('>>> Loading dataset...')
		dataset = pd.read_csv(args.dataset, compression = 'gzip')#.iloc[:,1:]
		print('>>> Writting pickle...')
		pickle.dump(dataset, open('pickle/data.p', 'wb'))
		print(dataset)
		name = args.dataset[20:]
	
	else:
		print('>>> No dataset or workspace provided!')

	if args.cleanVar:
		print('>>> Cleaning dataset...')
		dataset=variantSelection.cleanMissing(dataset, args.cleanVar)
		print('>>> Writing csv...')
		dataset.to_csv('../../data/datasets/cleaned_dataset.csv.gz', index=False, compression = 'gzip')
		print('>>> Writting pickle...')
		pickle.dump(dataset, open('pickle/dataCln.p', 'wb'))
	
	if args.imputeVar:
		print('>>> Imputation...')
		dataset=variantSelection.doImputation(dataset, args.imputeVar)
		print('>>> Writing csv...')
		dataset.to_csv('../../data/datasets/imputed_dataset.csv.gz', index=False, compression = 'gzip')
		print('>>> Writting pickle...')
		pickle.dump(dataset, open('pickle/dataImp.p', 'wb'))

	if args.classifier:
		X=dataset.iloc[:,:-1].values
		y=pd.to_numeric(dataset.iloc[:,-1].values.ravel())

		# Use correct tag
		name = name[12:-7]
		if name != 'pval' and name != 'network':
			name = name[4:]
			if name == 'pval':
				name = 'all_pval'
			elif name == 'network':
				name = 'all_network'
				
		models.featuresSel(dataset.iloc[:,:-1],pd.to_numeric(dataset.iloc[:,-1]),name)
		if args.classifier[0] == 'svm': plots.plotRocCurve(X, y, 'svm', name, args.grid_search)
		elif args.classifier[0] == 'tree': plots.plotRocCurve(X, y, 'tree', name, args.grid_search)
		elif args.classifier[0] == 'knn': plots.plotRocCurve(X, y, 'knn', name, args.grid_search)
		elif args.classifier[0] == 'log': plots.plotRocCurve(X, y, 'log', name, args.grid_search)
		elif args.classifier[0] == 'rf': plots.plotRocCurve(X, y, 'rf', name, args.grid_search)
		elif args.classifier[0] == 'nb': plots.plotRocCurve(X, y, 'nb', name)
		else: print('>>> The classifier chosen is not valid!')

	if args.test:
		if args.test == 'normality': geneSelection.normalityTest(dataset)
		elif args.test == 'chi2': geneSelection.chiSquaredTest(dataset)
		else: print('>>> The test chosen is not valid!')


	if args.translation_genes:
		with open(args.translation_genes, 'r') as f:
			f = f.read().replace('\n', ', ').split(',')
			varList = [i for i in f if i.startswith('chr')]
		geneSelection.geneTranslation(varList)

	elif args.add_genes:
		geneList = pd.read_csv(args.add_genes)
		print(geneList)
		data = geneSelection.addGenes(dataset, geneList, 'Genes')
		print('>>> Writing csv...')
		data.to_csv('../../data/datasets/dataset_with_genes.csv.gz', index=False, compression = 'gzip')
		print('>>> Writting pickle...')
		pickle.dump(data, open('pickle/dataGenes.p', 'wb'))

	if args.new_features: 
		featureExtraction.reduceFeatures(dataset, args.new_features)

	if args.top_features:
		name = name[16:-7]
		if args.top_features <= 100 and args.top_features >= 0: 
			featureExtraction.selectTopFeatures(dataset, name, args.top_features)
		else: 
			print('>>> The value must be between 0 and 100!')
	
	if args.network:
		network.editGenesFile()
		if len(args.network) == 1:
			rData = False
		else: rData = args.network[1]
		network.createNetwork(args.network[0], rData)



	print('>>> END')

if __name__ == '__main__':
	main()
