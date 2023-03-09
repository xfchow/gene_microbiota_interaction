#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 09:44:46 2023

@author: zhou
"""

############Spearman correlation analysis###########
import pandas as pd
LUSC_fpkm = pd.read_csv(r'df_LUSC_fpkm.txt',sep='\t',index_col = 0)
sig_microbe = pd.read_csv(r'LUSC_sig_micriobe.csv',index_col = 0)
#'LUSC_sig_micriobe' is the abundance data of microbiome obtained after LEfSe analysis

relation = []
correlation = []
pvalue = []
import scipy
for i in list(LUSC_fpkm.columns):
    for j in list(sig_microbe.columns):
        k = scipy.stats.spearmanr(LUSC_fpkm.loc[:,i], sig_microbe.loc[:,j])
        r = (i,j)
        relation.append(r)
        correlation.append(k[0])
        pvalue.append(k[1])
print(len(relation));print(len(correlation));print(len(pvalue))
gene = []
microbe = []
for i in range(len(relation)):
    gene.append(relation[i][0])
    microbe.append(relation[i][1])

import numpy as np
dic = {
    "microbe": microbe,  
    "correlation": correlation,  
    "pvalue": pvalue
}  
data = pd.DataFrame(dic, index=gene) 