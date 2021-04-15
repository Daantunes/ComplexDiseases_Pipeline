#!/usr/bin/env python3
# -*- coding:utf-8 -*-
'''
# File: biogridParser.py
# Created Date: Thursday, April 15th 2021
# Author: Debora Antunes
# -----
# Last Modified: Thursday, April 15th 2021, 11:02:06 am
# -----
'''
import pandas as pd


dt = pd.read_csv('human.txt', sep='\t', na_values='-', index_col=False)
print(dt)
dt.columns = ['ID(s) interactor A', 'ID(s) interactor B', 'Alt. ID(s) interactor A',
	   'Alt. ID(s) interactor B', 'Al_A',
	   'Al_B', 'Interaction detection method(s)',
	   'Publication 1st author(s)', 'Publication Identifier(s)',
	   'Taxid_A', 'Taxid_B', 'Interaction type(s)',
	   'Source database(s)', 'Interaction identifier(s)',
	   'Confidence value(s)', 'Expansion method(s)',
	   'Biological role(s) interactor A', 'Biological role(s) interactor B',
	   'Exp_A',
	   'Exp_B', 'Type(s) interactor A',
	   'Type(s) interactor B', 'Xref(s) interactor A', 'Xref(s) interactor B',
	   'Interaction Xref(s)', 'Annotation(s) interactor A',
	   'Annotation(s) interactor B', 'Interaction annotation(s)',
	   'Host organism(s)', 'Interaction parameter(s)', 'Creation date',
	   'Update date', 'Checksum(s) interactor A', 'Checksum(s) interactor B',
	   'Interaction Checksum(s)', 'Negative', 'Feature(s) interactor A',
	   'Feature(s) interactor B', 'Stoichiometry(s) interactor A',
	   'Stoichiometry(s) interactor B', 'Identification method participant A',
	   'Identification method participant B']
dt.replace('', np.nan, inplace=True)

dt = dt.loc[dt['Taxid_A'] == 'taxid:9606(human)|taxid:9606(Homo sapiens)']
dt = dt.loc[dt['Taxid_B'] == 'taxid:9606(human)|taxid:9606(Homo sapiens)']


A=dt.loc[:,'Al_A'].str.split('|')
B=dt.loc[:,'Al_B'].str.split('|')
idx=list(dt.index)
Symb_A, Symb_B =[],[]
for i in idx:
	for j in A[int(i)]:
		if j.startswith('psi-mi:') and j.endswith('(display_short)'):
			if j[7]=='"':
				Symb_A.append(None)
			else:
				Symb_A.append(j[7:-15])
			break
		elif j == A[int(i)][-1]:
			Symb_A.append(None)
	for j in B[int(i)]:
		if j.startswith('psi-mi:') and j.endswith('(display_short)'):
			if j[7]=='"':
				Symb_B.append(None)
			else:
				Symb_B.append(j[7:-15])
			break
		elif j == B[int(i)][-1]:
			Symb_B.append(None)

dt.loc[:,'Al_A']=Symb_A
dt.loc[:,'Al_B']=Symb_B

dt = dt.dropna(subset=['Al_A','Al_B','Taxid_A','Taxid_B','Exp_A','Exp_B'])[
	['Al_A','Al_B','Taxid_A','Taxid_B','Exp_A','Exp_B']]

	
dt.to_csv(r'interactions_simple.csv',
index=False, columns=['Al_A','Al_B'])