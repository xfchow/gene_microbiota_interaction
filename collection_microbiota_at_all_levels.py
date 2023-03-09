#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 09:04:11 2023

@author: zhou
"""

#############Collection of microbial data at all levels##########
import pandas as pd
import numpy as np
import re
LUAD_microbe = pd.read_csv('LUAD_microbe.csv',index_col = 0)
LUAD_microbe = LUAD_microbe.T
microbe = list(LUAD_microbe.columns)
phylum=[]; order=[]; Class=[]; family = []

for i in microbe:
    if re.findall(r”\.p__(.+?)\.”,i):
        phylum.append(re.findall(r”\.p__(.+?)\.”,i)[0])
    else:
        phylum.append(0)

for i in microbe:
    if re.findall(r”\.o__(.+?)\.”,i):
        order.append(re.findall(r”\.o__(.+?)\.”,i)[0])
    else:
        order.append(0)

for i in microbe:
    if re.findall(r"\.f__(.+?)\.",i):
        Class.append(re.findall(r"\.f__(.+?)\.",i)[0])
    else:
        Class.append(0)

for i in microbe:
    if re.findall(r"\.f__(.+?)\.",i):
        family.append(re.findall(r"\.f__(.+?)\.",i)[0])
    else:
        family.append(0)

LUAD_phylum = LUAD_microbe.T
LUAD_phylum[‘phylum’] = phylum
LUAD_phylum = LUAD_phylum.set_index(‘phylum’,drop = True)
LUAD_phylum = LUAD_phylum.drop(0,axis=1)
df_phylum = LUAD_phylum.groupby([‘phylum’]).sum()

LUAD_order = LUAD_microbe.T
LUAD_order[‘order’] = order
LUAD_order = LUAD_order.set_index(‘order’,drop = True)
LUAD_order = LUAD_order.drop(0,axis=1)
df_order = LUAD_order.groupby([‘order’]).sum()

LUAD_class = LUAD_microbe.T
LUAD_class[‘class’] = Class
LUAD_class = LUAD_class.set_index(‘class’,drop = True)
LUAD_class = LUAD_class.drop(0,axis=1)
df _class= LUAD_class.groupby([‘class’]).sum()

LUAD_family = LUAD_microbe.T
LUAD_family['family'] = family
LUAD_family = LUAD_family.set_index('family',drop = True)
LUAD_family = LUAD_family.drop(0,axis=1)
df_family = LUAD_family.groupby([‘family’]).sum()

############Spearman correlation analysis###########
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



