#!/usr/bin/env python3
# -*- coding:utf-8 -*-
'''
# File: network.py
# Created Date: Thursday, January 11th 2021
# Author: Debora Antunes
# -----
# Last Modified: Thursday, April 15th 2021, 11:43:55 am
# -----
'''

import pandas as pd 
import numpy as np

from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
import logging
rpy2_logger.setLevel(logging.ERROR)

def editGenesFile():
	"""Transforms values from the original file into values acepted by the R package
	"""
	print('>>> Transforming genes file...')
	geneList = pd.read_csv('../../data/genes/geneList.csv')
	varList = pd.read_csv('../../data/variants/sigVars.csv')
	result = pd.merge(geneList,varList, on=['Variants'])
	result = result.groupby(['Genes'])['pval'].min().reset_index() ## more significant
	result.replace(1, 0.99999, inplace=True) ### PREVENT FUTURE ERROR
	p_vals = result[result['pval']>=1e-16]
	num_min = p_vals.min(axis=0)
	for i in range(len(result['pval'])):
		if result.loc[i,'pval']<1e-16:
			result.loc[i,'pval']= num_min['pval']
	result.to_csv('../../data/proteins/genes.csv',index=False, columns=['Genes' , 'pval'], header=['gene','weight'])

def createNetwork(interaction, rData = False):
	pandas2ri.activate()
	r_script = ro.r
	r_script['source'](r'../R/network.R')
	r_script.fun(interaction, rData)

