#!/usr/bin/env python3
# -*- coding:utf-8 -*-
'''
# File: models.py
# Created Date: Wednesday February 12th 2020
# Author: Debora Antunes
# -----
# Last Modified: Thursday, April 15th 2021, 10:42:36 am
# -----
'''
import pickle
import pandas as pd
import numpy as np
import csv

import plots
import collections

from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from feature_selector import FeatureSelector
from sklearn.metrics import accuracy_score
from sklearn.metrics import f1_score

from matplotlib import pyplot as plt

from progress.bar import Bar

def trainSvm(train_X,train_y,test_X,test_y,best_params):
	"""Uses a SVM model to fit the data.
	
	Arguments:
		train_X {numpy.ndarray} -- 0.8% of the original dataset for training
		train_y {numpy.ndarray} -- labels for the training data
		test_X {numpy.ndarray} -- 0.2% of the original dataset for testing
		test_y {numpy.ndarray} -- labels for the testing data
		best_params {dictionary} -- hyperparametes
	"""	
	if not best_params:
		print('>>> Starting grid search...')
		parameters = {'kernel':['linear'], 
					'C':[0.25,0.4,0.5,0.55,0.75,1], 
					'tol':[1e-3,1e-4,1e-5], 
					'gamma':[25,50,75,100,150,'auto'], 
					'degree':[1,2,3,5,10]}
		svm_model = SVC()
		grid = GridSearchCV(svm_model, parameters, cv = 5, scoring='f1')
		grid_result = grid.fit(train_X, train_y)
		best_params = grid_result.best_params_
		print(best_params)
	print('>>> Creating model...')
	svm_model = SVC(kernel=best_params["kernel"],
					C=best_params["C"], 
					tol=best_params["tol"], 
					gamma=best_params["gamma"], 
					degree=best_params["degree"])
	svm_model = svm_model.fit(train_X, train_y)
	prev = svm_model.predict(test_X)
	print(accuracy_score(test_y,prev))
	print(f1_score(test_y,prev, average='macro'))
	return svm_model, list(prev)

def trainKnn(train_X,train_y,test_X,test_y,best_params):
	"""Uses a KNN model to fit the data.
	
	Arguments:
		train_X {numpy.ndarray} -- 0.8% of the original dataset for training
		train_y {numpy.ndarray} -- labels for the training data
		test_X {numpy.ndarray} -- 0.2% of the original dataset for testing
		test_y {numpy.ndarray} -- labels for the testing data
		best_params {dictionary} -- hyperparametes
	"""
	if not best_params:	
		print('>>> Starting grid search...')
		parameters = {'n_neighbors':[3,5,7,9,11],
					'leaf_size':[1,3,5,10,15,20,30],
					'algorithm':['auto', 'kd_tree', 'ball_tree', 'brute']}
		knn_model = KNeighborsClassifier()
		grid = GridSearchCV(knn_model, parameters, cv = 5, scoring='f1')
		grid_result = grid.fit(train_X, train_y)
		best_params = grid_result.best_params_
		print(best_params)
	print('>>> Creating model...')
	knn_model = KNeighborsClassifier(n_neighbors = best_params['n_neighbors'],
									leaf_size = best_params['leaf_size'],
									algorithm = best_params['algorithm'])
	knn_model = knn_model.fit(train_X, train_y)
	prev = knn_model.predict(test_X)
	print(accuracy_score(test_y,prev))
	print(f1_score(test_y,prev, average='macro'))
	return knn_model, list(prev)

def trainLog(train_X,train_y,test_X,test_y,best_params):
	"""Uses a Logarithmic model to fit the data.
	
	Arguments:
		train_X {numpy.ndarray} -- 0.8% of the original dataset for training
		train_y {numpy.ndarray} -- labels for the training data
		test_X {numpy.ndarray} -- 0.2% of the original dataset for testing
		test_y {numpy.ndarray} -- labels for the testing data
		best_params {dictionary} -- hyperparametes
	"""	
	if not best_params:
		print('>>> Starting grid search...')
		parameters = {'penalty': ['l1', 'l2'],
					'C':[0.25,0.4,0.5,0.55,0.75,1], 
					'tol':[1e-3,1e-4,1e-5],
					'solver':['liblinear']}
		log_model = LogisticRegression()
		grid = GridSearchCV(log_model, parameters, cv = 5, scoring='f1')
		grid_result = grid.fit(train_X, train_y)
		best_params = grid_result.best_params_
		print(best_params)	
	print('>>> Creating model...')
	log_model = LogisticRegression(penalty = best_params['penalty'],
									C = best_params['C'], 
									tol = best_params['tol'],
									solver = best_params['solver'])
	log_model = log_model.fit(train_X, train_y)
	prev = list(log_model.predict(test_X))
	print(accuracy_score(test_y,prev))
	print(f1_score(test_y,prev, average='macro'))
	return log_model, list(prev)

def trainNb(train_X,train_y,test_X, test_y):
	"""
	Arguments:
		train_X {numpy.ndarray} -- 0.8% of the original dataset for training
		train_y {numpy.ndarray} -- labels for the training data
		test_X {numpy.ndarray} -- 0.2% of the original dataset for testing
		test_y {numpy.ndarray} -- labels for the testing data
	"""	
	print('>>> Creating model...')
	nb_model = GaussianNB()
	nb_model = nb_model.fit(train_X, train_y)
	prev = nb_model.predict(test_X)
	print(accuracy_score(test_y,prev))
	print(f1_score(test_y,prev, average='macro'))
	return nb_model, list(prev)

def trainRf(train_X,train_y,test_X,test_y,best_params):
	"""Uses a Random forest model to fit the data.
	
	Arguments:
		train_X {numpy.ndarray} -- 0.8% of the original dataset for training
		train_y {numpy.ndarray} -- labels for the training data
		test_X {numpy.ndarray} -- 0.2% of the original dataset for testing
		test_y {numpy.ndarray} -- labels for the testing data
		best_params {dictionary} -- hyperparametes
	"""	
	if not best_params:
		print('>>> Starting grid search...')
		parameters = {'n_estimators':[20,50,75,100], 
				'criterion':['entropy','gini'], 
				'min_samples_leaf':[1,2,3,5,10], 
				'min_samples_split':[2,4,5,8,10], 
				'max_leaf_nodes':[2,20,50,75,100]}
		rf_model = RandomForestClassifier()
		grid = GridSearchCV(rf_model, parameters, cv = 5, scoring='f1')
		grid_result = grid.fit(train_X, train_y)
		best_params = grid_result.best_params_
		print(best_params)	
	print('>>> Creating model...')
	rf_model = RandomForestClassifier(criterion=best_params["criterion"], 
				max_leaf_nodes=best_params["max_leaf_nodes"], 
				min_samples_leaf=best_params["min_samples_leaf"], 
				min_samples_split=best_params["min_samples_split"], 
				n_estimators=best_params["n_estimators"])
	rf_model = rf_model.fit(train_X, train_y)
	prev = rf_model.predict(test_X)
	print(accuracy_score(test_y,prev))
	print(f1_score(test_y,prev, average='macro'))
	return rf_model, list(prev)


def trainTree(train_X,train_y,test_X,test_y,best_params):
	"""Uses a tree model to fit the data.
	
	Arguments:
		train_X {numpy.ndarray} -- 0.8% of the original dataset for training
		train_y {numpy.ndarray} -- labels for the training data
		test_X {numpy.ndarray} -- 0.2% of the original dataset for testing
		test_y {numpy.ndarray} -- labels for the testing data
		best_params {dictionary} -- hyperparametes
	"""	
	if not best_params:
		print('>>> Starting grid search...')
		parameters = {'n_estimators':[20,50,75,100], 
					'criterion':['entropy','gini'], 
					'min_samples_leaf':[1,2,3,5,10], 
					'min_samples_split':[2,4,5,8,10], 
					'max_leaf_nodes':[2,20,50,75,100]}
		tree_model = ExtraTreesClassifier()
		grid = GridSearchCV(tree_model, parameters, cv = 5, scoring='f1')
		grid_result = grid.fit(train_X, train_y)
		best_params = grid_result.best_params_
		print(best_params)	
	print('>>> Creating model...')
	tree_model = ExtraTreesClassifier(criterion=best_params["criterion"], 
				max_leaf_nodes=best_params["max_leaf_nodes"], 
				min_samples_leaf=best_params["min_samples_leaf"], 
				min_samples_split=best_params["min_samples_split"], 
				n_estimators=best_params["n_estimators"])
	tree_model = tree_model.fit(train_X, train_y)
	prev = tree_model.predict(test_X)
	print(accuracy_score(test_y,prev))
	print(f1_score(test_y,prev, average='macro'))
	return tree_model, list(prev)


def useTree(data,train_X,train_y):
	"""Uses an Extremely Randomized Trees ensemble to fit the data.
	
	Arguments:
		data {pandas.Dataframe} -- Original dataset
		train_X {numpy.ndarray} -- Original dataset for training
		train_y {numpy.ndarray} -- Labels for the training data
	"""
	if train_X.shape[1] < 100: thr = train_X.shape[1]-1
	else: thr = 99

	tree_model = ExtraTreesClassifier(100, n_jobs = -1)
	tree_model = tree_model.fit(train_X, train_y)
	importance = np.array(tree_model.feature_importances_)
	select = (importance > sorted(importance, reverse=True)[thr])
	dataset_sel = data.iloc[:, select]
	return list(dataset_sel.columns)

def featuresSel(train, train_labels, name):
	"""Plots the curve for the importantant features
	
	Arguments:
		train {pandas.Dataframe} -- Dataset
		train_labels {numpy.ndarray} -- Labels for the dataset
		name {string} -- Name for file
	"""
	print('>>> Feature Selection...')
	fs = FeatureSelector(data = train, labels = train_labels)
	fs.identify_zero_importance(task = 'classification', 
                            eval_metric = 'auc', 
                            n_iterations = 10, 
                             early_stopping = True)
	plt.figure(figsize=(15,15))
	fs.plot_feature_importances(threshold = 0.99, plot_n = 50, name = name)
	plt.savefig('../../data/figures/rank_{}.png'.format(name))
	plt.close()

	