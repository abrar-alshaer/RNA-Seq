rm(list=ls())
library(dplyr)
library(DEXSeq)
library(BiocParallel)
library(GenomicRanges)
library(data.table)
source("load_SubreadOutput.R") #load subread (featurecounts) to DEXSeq R code
suppressPackageStartupMessages(library("DEXSeq"))
setwd("D:/OneDrive - University of North Carolina at Chapel Hill/scripts/Smith Lab/Wei Sha Analysis/Abrar Analysis/Illumia Data/RERUN/Galgal5Run/")

#MAKE SURE TO ALWAYS INITIALIZE CONTROLS FIRST (in code below) THEN CASES (SALINE, ETOH) & MAKE SURE IN YOUR FEATURE COUNTS OUTPUT YOU HAVE YOUR CONTROLS FIRST THEN YOUR CASE COUNTS (SALINE COUNTS FIRST THEN ETOH COUNTS)
#If you initialize your cases before your controls in DEXSeq code below but it's the reverse order in your counts files you will get the incorrect counts assigned to your conditions.
#If you initialize your cases before your controls in both your DEXSeq code below and your featurecounts output file then the interpretation of your fold change will be different. 
samp <- data.frame(row.names = c("saline_1","saline_2","ETOH_1","ETOH_2"), condition = rep(c("SALINE","ETOH"),each=2))
dxd.fc <- DEXSeqDataSetFromFeatureCounts("ALL_RUN-182-183.SB_ETOH_SALINE_EXONS-DEXSeq.counts", flattenedfile = "featurecounts.gtf",sampleData = samp)
#Run DEXSeq#
dxr = DEXSeq(dxd.fc)
#ALL_RUN-182-183.SB_ETOH_SALINE_EXONS-DEXSeq_unstranded-B.counts

#DON'T RUN IF RAN LINE: dxr = DEXSeq(dxd.fc)
#The lines below are all the substeps DEXSeq does when running the mass command above
###############################################
dxd.fc = estimateSizeFactors(dxd.fc)
dxd.fc = estimateDispersions(dxd.fc)
plotDispEsts(dxd.fc)

dxd.fc = testForDEU(dxd.fc)
dxd.fc = estimateExonFoldChanges(dxd.fc, fitExpToVar="condition")
dxr1 = DEXSeqResults(dxd.fc)
results_annotation <- mcols(dxr1)$description
#################################################

#From this object, we can ask how many exonic regions are significant with a false discovery rate of 10%
table (dxr$padj < 0.1)
#We may also ask how many genes are affected
table (tapply(dxr$padj < 0.1, dxr$groupID, any))

#plots the logarithm of fold change versus average normalized count per exon and marks by red colour 
#the exons which are considered significant; here, the exons with an adjusted p values of less than 0.1
plotMA( dxr, cex=0.8 )

#order results (res) by pvalue
res <- dxr[order(dxr$padj), ]
filter <- na.omit(dxr) #remove genes with NA values
resCutoff <- filter[(filter$padj)<0.1, ]

#################### NORMALIZED COUNTS
df <- data.frame(resCutoff[,c(0,7,10)]) #initializing a dataframe with the rows as the first column, adjusted pvalue, and log2FC
df_rows_1 <- setDT(df, keep.rownames = TRUE)[] #converting rows to first column

counts <- counts(dxr, normalized=TRUE) #retrieving normalized DEXSeq counts
results_counts <- as.data.frame(counts) #converting the counts to a dataframe
counts_rmNA <- na.omit(results_counts) #remove genes with NA values
row_sub = apply(counts_rmNA, 1, function(row) all(row !=0 )) #removing zeroes 
counts_noZero <- counts_rmNA[row_sub,] #removing zeroes 
counts_noZero_rows_1 <- setDT(counts_noZero, keep.rownames = TRUE)[] #setting the row names as first column

counts_noZero$padj <- df_rows_1$padj[match(counts_noZero_rows_1$rn, df_rows_1$rn)] #adding the padj column to the normalized counts columns
counts_noZero$log2fold_SALINE_ETOH <- df_rows_1$log2fold_SALINE_ETOH[match(counts_noZero_rows_1$rn, df_rows_1$rn)] #adding the log2FC column to the normalized counts columns

counts_rmNA_2 <- na.omit(counts_noZero) #remove genes with NA values
countsCutoff <- counts_rmNA_2[(counts_rmNA_2$padj)<0.1, ] #retrieving all the exons with an adjusted pval less than 0.1
#output normalized counts with log2FC and p-value information
write.csv(countsCutoff, file="DEXSeqCounts_run182_183.SB_ETOH_SALINE_All_Lanes_unstranded-B_pval0.1.csv")
############################# NORMALIZED COUNTS

#BELOW the results outputted into the files DO NOT include normalized counts
#DEXSeq Results in dataframe
results_data <- as.data.frame(res)
resdata_rmNA <- na.omit(results_data) #remove genes with NA values

#write to CSV #BELOW the results outputted into the files DO NOT include normalized counts
write.csv(results_data, file="DEXSeqResults_run182_183.SB_ETOH_SALINE_All_Lanes_unstranded-B.csv")
write.csv(resdata_rmNA, file="DEXSeqResults_run182_183.SB_ETOH_SALINE_All_Lanes_RemovedNA_unstranded-B.csv")
write.csv(resCutoff, file="DEXSeqJUNC_run182_183.SB_ETOH_SALINE_All_Lanes_rmNA_0.1Pvalcut_unstranded-B.csv")

#initialize a high resolution graph
png("ENSGALG00000036474+ENSGALG00000003959.png", units="in", width=11, height=9, res=600)

#Visualization
plotDEXSeq(dxr, "ENSGALG00000036474+ENSGALG00000003959", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

#release PNG file graph into directory
dev.off() 
