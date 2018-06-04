#!/usr/bin/env python
#This program takes the DEXSeq output from the column groupID (gene IDs for exons) and parses the gene IDs to put them into an individual list then get a unique set of that list.

#genomicFile = 'DEXSeq-run182_183.SB_ETOH_SALINE_All_Lanes_RemovedNA_unstr-B_groupID_PosFC.txt' #'TEST-negFC_ALLGenes.txt'
genomicFile = 'DAVID_TEST-KEGGs_DEXSeqResults_neg-posFC_padj0.1_unstr-B_rmNA.txt'

with open(genomicFile) as lines:
	genes = lines.readlines()
lines.close()

IDs = []
for i in genes:
    if ',' in i:
        for j in range(0, len(i.split(','))):
            IDs.append((i.split(',')[j]).strip())
    else:
        IDs.append(i.strip())

unique = list(set(IDs))

print "Total length:", len(IDs)
print len(set(IDs))
print len(unique)
print unique

thefile = open('groupKEGGs_UNIQUE_DEXSeqResults_neg-posFC_padj0.1_unstr-B_rmNA.txt', 'w')
for item in unique:
  thefile.write("%s\n" % item) #target.write('%r\n%r\n%r\n' % (line1, line2, line3))
'''
#######################################################################################
genomicFile2 = 'DEXSeq-run182_183.SB_ETOH_SALINE_All_Lanes_RemovedNA_unstr-B_groupID_NegFC.txt' #'TEST-posFC_ALLGenes.txt'
#'groupID_DEXSeqResults_negFC_padj0.1.txt'

with open(genomicFile2) as lines2:
	genes2 = lines2.readlines()
lines2.close()

IDs2 = []
for i in genes2:
    if '+' in i:
        for j in range(0, len(i.split('+'))):
            IDs2.append((i.split('+')[j]).strip())
    else:
        IDs2.append(i.strip())

unique2 = list(set(IDs2))

print 'Total length:', len(IDs2)
print len(set(IDs2))
print len(unique2)

print "DIFFERENCE", len(list(set(unique)-set(unique2)))

afile = open('groupID_UNIQUE_DEXSeqResults_negFC_padj0.1_unstr-B.txt', 'w')
for item in unique2:
  afile.write("%s\n" % item)
'''
