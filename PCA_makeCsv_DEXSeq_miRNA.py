#!/usr/bin/env python
#This program makes a master CSV file for the PCA analysis from the DAVID KEGG info and the DEXSeq results
import pandas as pd
from pandas import DataFrame
import numpy as np
import os, sys

genomicFile = 'DEXSeq_pval_less_0.1_miRNAs_rmNan.csv'
DEXseq = 'DEXSeqCounts_run182_183.SB_ETOH_SALINE_All_Lanes_unstranded-B_pval0.1.csv'

flag = []
from_csv = pd.read_csv(genomicFile, sep=None, header=0, engine='python')
flag = from_csv[['Gene stable ID','miRBase ID',]] #reads header columns
print flag.head()

flag_d = []
from_csv_d = pd.read_csv(DEXseq, sep=None, header=0, engine='python')
flag_d = from_csv_d[['Exons','saline_1','saline_2','ETOH_1','ETOH_2','log2fold_SALINE_ETOH']] #reads in fourth header column and fifth header column
print flag_d.head()
#replace 'Exons' with 'Genes' in DESeq gene level analysis

genes = []
pathways = []

for index, row in flag.iterrows():
    print "ROW:", row[0], row [1]
    if ',' in row[1]:
        for j in range(0, len(row[1].split(','))):
            genes.append(row[0].strip()[:-1])
            pathways.append((row[1].split(',')[j]).strip())
            #print row[0].strip()[:-1],(row[1].split(',')[j]).strip()
    else:
        genes.append(row[0].strip())
        pathways.append(row[1].strip())
        #print row[0].strip()[:-1],row[1].strip()

exons = []
log2fold_SALINE_ETOH = []
countDataSaline_1 = []
countDataSaline_2 = []
countDataETOH_1 = []
countDataETOH_2 = []
genes_d = []
pathways_d = []

#dict(zip(genes, pathways)) #combines two lists into a dictionary

for index, row in flag_d.iterrows():
    for k in range(0, len(genes)):
        if genes[k] in row[0]:
            print 'genes:', genes[k], row[0], row[1], row[2], row[3], row[4], row[5]
            genes_d.append(genes[k])
            pathways_d.append(pathways[k])
            exons.append(row[0]) #replace with DESeq gene ID in gene level analysis
            log2fold_SALINE_ETOH.append(row[5])
            countDataSaline_1.append(row[1])
            countDataSaline_2.append(row[2])
            countDataETOH_1.append(row[3])
            countDataETOH_2.append(row[4])

df = pd.DataFrame({'GeneID': genes_d, 'miRBase ID': pathways_d, 'Exons': exons, 'log2fold_SALINE_ETOH': log2fold_SALINE_ETOH, 'countDataSaline_1': countDataSaline_1, 'countDataSaline_2': countDataSaline_2, 'countDataETOH_1': countDataETOH_1, 'countDataETOH_2': countDataETOH_2}) #creates a dataframe
print df.head()

df.to_csv('DEXSeq_pval_less_0.1_miRNAs_PCA_file.csv')
