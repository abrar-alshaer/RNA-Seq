#!/usr/bin/env python
#This program makes a master CSV file for the PCA analysis from the DAVID KEGG info and the DEXSeq results
import pandas as pd
from pandas import DataFrame
import numpy as np
import os, sys

genomicFile = 'DAVID_TEST_DEXSeqResults_neg-posFC_padj0.1_unstr-B_rmNA_filterKEGGs.csv'
DEXseq = 'DEXSeqResults_run182_183.SB_ETOH_SALINE_All_Lanes_rmNA_0.1unstranded-B.csv'

flag = []
from_csv = pd.read_csv(genomicFile, sep=None, header=0, engine='python')
flag = from_csv[['ENSEMBL_GENE_ID','KEGG_PATHWAY',]] #reads header columns
print flag.head()

flag_d = []
from_csv_d = pd.read_csv(DEXseq, sep=None, header=0, engine='python')
flag_d = from_csv_d[['Exons','exonBaseMean','ETOH','SALINE','log2fold_SALINE_ETOH','countData.saline_1','countData.saline_2','countData.ETOH_1','countData.ETOH_2']] #reads in fourth header column and fifth header column
print flag_d.head()
#replace 'Exons' with 'Genes' in DESeq gene level analysis

genes = []
pathways = []

for index, row in flag.iterrows():
    #print "ROW:", row[0], row [1]
    if ',' in row[1]:
        for j in range(0, len(row[1].split(','))):
            genes.append(row[0].strip()[:-1])
            pathways.append((row[1].split(',')[j]).strip())
            #print row[0].strip()[:-1],(row[1].split(',')[j]).strip()
    else:
        genes.append(row[0].strip()[:-1])
        pathways.append(row[1].strip())
        #print row[0].strip()[:-1],row[1].strip()

exons = []
exonBaseMean = []
ETOH = []
SALINE = []
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
            print 'genes:', genes[k], row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8]
            genes_d.append(genes[k])
            pathways_d.append(pathways[k])
            exons.append(row[0]) #replace with DESeq gene ID in gene level analysis
            exonBaseMean.append(row[1])
            ETOH.append(row[2])
            SALINE.append(row[3])
            log2fold_SALINE_ETOH.append(row[4])
            countDataSaline_1.append(row[5])
            countDataSaline_2.append(row[6])
            countDataETOH_1.append(row[7])
            countDataETOH_2.append(row[8])

df = pd.DataFrame({'GeneID': genes_d, 'KEGG_PATHWAY': pathways_d, 'Exons': exons, 'exonBaseMean': exonBaseMean, 'ETOH': ETOH, 'SALINE': SALINE, 'log2fold_SALINE_ETOH': log2fold_SALINE_ETOH, 'countDataSaline_1': countDataSaline_1, 'countDataSaline_2': countDataSaline_2, 'countDataETOH_1': countDataETOH_1, 'countDataETOH_2': countDataETOH_2}) #creates a dataframe
print df.head()

df.to_csv('DAVID_TEST_DEXSeqResults_neg-posFC_padj0.1_unstr-B_rmNA_filterKEGGs_1E1P.csv')
