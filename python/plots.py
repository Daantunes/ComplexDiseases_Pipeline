#!/usr/bin/env python3
# -*- coding:utf-8 -*-
'''
# File: plots.py
# Created Date: Wednesday March 25th 2020
# Author: Debora Antunes
# -----
# Last Modified: Thursday, April 15th 2021, 11:28:42 am
# -----
'''
import numpy as np
import pandas as pd

import models

from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
import logging
rpy2_logger.setLevel(logging.ERROR)

from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import confusion_matrix, roc_curve, auc, plot_roc_curve
from scipy.interpolate import interpn

import matplotlib.pyplot as plt
import seaborn as sns



def plotManhattan(data, name):
	"""Uses a R script to create a manhattan plot
	
	Arguments:
		data {pandas.Dataframe} -- Dataset used to plot
		name {string} -- Sufix for plot file
	"""	
	print('>>> Creating figure...')
	pandas2ri.activate()
	r_script = ro.r
	r_script['source'](r'../R/plotManhattan.R')
	r_script.fun(data, name)

def plotConfusionMatrix(real,prev,dir,model):
	"""Image for the confusion matrix
	
	Arguments:
		real {list} -- Real labels
		prev {list} -- Predicted labels
		dir {string} -- Name of the directory
		model {string} -- add the model name to the file name
	"""	
	print('>>> Creating figure...')

	fig = plt.figure()
	plt.plot([2,2,4]) 
	ax0 = plt.subplot(2, 2, 1)
	ax0 = sns.heatmap(confusion_matrix(real[0],prev[0]), annot=True, cbar = False, cmap=sns.light_palette('#dea369'))
	ax0.set(xlabel='Predicted', ylabel='Real')
	ax0.title.set_text('Confusion Matrix Fold 1')
	ax1 = plt.subplot(2, 2, 2) 
	ax1 = sns.heatmap(confusion_matrix(real[1],prev[1]), annot=True, cbar = False, cmap=sns.light_palette('#dea369'))
	ax1.set(xlabel='Predicted', ylabel='Real')
	ax1.title.set_text('Confusion Matrix Fold 2')
	ax2 = plt.subplot(2, 2, 3) 
	ax2 = sns.heatmap(confusion_matrix(real[2],prev[2]), annot=True, cbar = False, cmap=sns.light_palette('#dea369'))
	ax2.set(xlabel='Predicted', ylabel='Real')
	ax2.title.set_text('Confusion Matrix Fold 3')
	ax3 = plt.subplot(2, 2, 4) 
	ax3 = sns.heatmap(confusion_matrix(real[3],prev[3]), annot=True, cbar = False, cmap=sns.light_palette('#dea369'))
	ax3.set(xlabel='Predicted', ylabel='Real')
	ax3.title.set_text('Confusion Matrix Fold 4')

	fig.tight_layout()
	print('>>> Saving figure...')
	plt.savefig('../../data/figures/{}/matrix_{}.png'.format(dir, model))

def plotRocCurve(X,y, name, dir, best=False):
	"""Image for the ROC curve
	
	Arguments:
		X {numpy.ndarray} -- Dataset to train
		y {numpy.ndarray} -- Labels for the dataset
		name {string} -- Name for the file
		dir {string} -- Name of the directory
		best {list} -- Hyperparameters
	"""	
	print(dir, name)
	skf = StratifiedKFold(n_splits=5, shuffle = True)
	tprs, aucs, real, prev = [], [], [], []
	mean_fpr = np.linspace(0, 1, 100)

	# Parse hyperparameters if provided
	if best:
		hyper = {}
		for i in best:
			key, val = i.split(':')[0], i.split(':')[1]
			try: 
				if val == str(float(val)):  val = float(val)
				else: val = int(float(val))
			except: pass
			finally: hyper[key]=val
		best = hyper
		print(best)

	# Cross-Validation
	fig, ax = plt.subplots()
	for i, (train, test) in enumerate(skf.split(X, y)):
		if name == 'svm': model, p = models.trainSvm(X[train], y[train], X[test], y[test], best)
		elif name == 'tree': model, p = models.trainTree(X[train], y[train], X[test], y[test], best)
		elif name == 'knn': model, p, best = models.trainKnn(X[train], y[train], X[test], y[test], best)
		elif name == 'log': model, p = models.trainLog(X[train], y[train], X[test], y[test], best)
		elif name == 'rf': model, p, best = models.trainRf(X[train], y[train], X[test], y[test], best)
		elif name == 'nb': model, p = models.trainNb(X[train], y[train], X[test], y[test])
		prev.append(p)
		real.append(list(y[test]))
		viz = plot_roc_curve(model, X[test], y[test],
                         name='ROC fold {}'.format(i+1),
                         alpha=0.3, lw=1, ax=ax)
		interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
		interp_tpr[0] = 0.0
		tprs.append(interp_tpr)
		aucs.append(viz.roc_auc)

	# ROC Curves
	print('>>> Creating figure...')
	ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r', label='Chance', alpha=.8)

	mean_tpr = np.mean(tprs, axis=0)
	mean_tpr[-1] = 1.0
	mean_auc = auc(mean_fpr, mean_tpr)
	std_auc = np.std(aucs)
	ax.plot(mean_fpr, mean_tpr, color='b',
			label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
			lw=2, alpha=.8)

	std_tpr = np.std(tprs, axis=0)
	tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
	tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
	ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
					label=r'$\pm$ 1 std. dev.')

	ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05], title="Receiver operating characteristic")
	ax.legend(loc="lower right")
	print('>>> Saving figure...')
	plt.savefig('../../data/figures/{}/roc_auc_{}.png'.format(dir, name))
	plt.close()

	plotConfusionMatrix(real,prev, dir, name)



def plotFeatures(labels, values, name):
	"""Image for the top features frequencies
	
	Arguments:
		labels {list} -- Selected features
		values {list} -- Frequence of the selected features
		name {string} -- add the dataset name to the file name
	"""	
	print('>>> Creating figure...')
	fig, ax = plt.subplots()
	ax.bar(labels, values)
	ax.set_xticklabels(labels, rotation = 90)
	fig.set_size_inches(20, 12)
	print('>>> Saving figure...')
	ax.figure.savefig('../../data/figures/features_{}.png'.format(name))
	plt.close()