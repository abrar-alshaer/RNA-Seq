setwd("D:/OneDrive - University of North Carolina at Chapel Hill/scripts/Smith Lab/Wei Sha Analysis/Abrar Analysis/Illumia Data/RERUN/Galgal5Run/miRNA database files/DAVID_R/")
rm(list=ls())
require(rJava)
library(dplyr)
library(stringi)
library("RDAVIDWebService") 
#If rJava doesn't load then download correct Java (x32 or x64) & set it to the correct JRE folder name and reload the rJava library
#Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre1.8.0_161') 

#genome <- read.csv("Galgal5.0-geneIDs_parse3.txt", header = FALSE) #full chicken genome from ensembl
genome_TS <- read.csv("mart_export_galgal_targetScan_geneID.txt", header = TRUE) #UTR gene targets from TargetScan Chicken Genome

david <- DAVIDWebService(email='email@unc.edu', url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
BG <- addList(david, genome_TS$Gene.stable.ID, idType="ENSEMBL_GENE_ID", listName="all", listType="Background")

datalist = list()
for(i in 1:200) #loop the number of times you want to randomly generate KEGG annotations
{
print(i)
str <- stri_rand_strings(1, 10, pattern = "[A-Za-z0-9]")
#To randomly sample rows from a dataframe
genome2 <- sample_n(genome_TS, 800)
FG <- addList(david, genome2$Gene.stable.ID, idType="ENSEMBL_GENE_ID", listName=str, listType="Gene")
setAnnotationCategories(david, c("KEGG_PATHWAY")) #setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
FuncAnnotChart <- getFunctionalAnnotationChart(david)
datalist[[i]] <- FuncAnnotChart
}

big_data = do.call(rbind, datalist) #combines all previous dataframes from for loop
#"gga04810:Regulation of actin cytoskeleton", "gga04530:Tight junction", "gga04520:Adherens junction", "gga04540:Gap junction"
pathway <- big_data[is.element(big_data$Term, c("gga04510:Focal adhesion", "gga04514:Cell adhesion molecules (CAMs)")),] #place specific pathway(s)
pathway_shh <- big_data[is.element(big_data$Term, c("gga04340:Hedgehog signaling pathway")),] #place specific pathway(s)
hist(pathway$PValue) #histogram 
hist(pathway_shh$PValue) #histogram 
dotchart(pathway$PValue) #dot plot

write.csv(big_data, "Random_KEGGs_miR-3533_5.csv")

#getFunctionalAnnotationChartFile(david, "comboChart.tsv")

#I ran this program 5 times over 5 days (not to exceed DAVID API limit) and then combined all 5 files to visualize the pvalues for 1000 iterations across all the files (200 iterations for 5 files = 1000). 
k1 <- read.csv("Random_KEGGs_miR-3533_1.csv", header = TRUE)
k2 <- read.csv("Random_KEGGs_miR-3533_2.csv", header = TRUE)
k3 <- read.csv("Random_KEGGs_miR-3533_3.csv", header = TRUE)
k4 <- read.csv("Random_KEGGs_miR-3533_4.csv", header = TRUE)
k5 <- read.csv("Random_KEGGs_miR-3533_5.csv", header = TRUE)

ktotal <- rbind(k1, k2, k3, k4, k5)
#"gga04810:Regulation of actin cytoskeleton", "gga04530:Tight junction", "gga04520:Adherens junction", "gga04540:Gap junction"
pathway_k <- ktotal[is.element(ktotal$Term, c("gga04510:Focal adhesion", "gga04514:Cell adhesion molecules (CAMs)")),] #place specific pathway(s)
hist(pathway_k$PValue, main = "Randomly Generated miRNA Gene Targets - Cell Adhesions \n800 Targets", xlab = "P-values")


pathway_shh_k <- ktotal[is.element(ktotal$Term, c("gga04340:Hedgehog signaling pathway")),]
hist(pathway_shh_k$PValue, main = "Randomly Generated miRNA Gene Targets - SHH \n800 Targets", xlab = "P-values")

myData <- read.csv("Cell_Adhesions_KEGGs_All_miRNAs_myData.csv", header = TRUE)
hist(myData$PValue, main = "miRNA Gene Targets RNASeq Dataset - Cell Adhesions", xlab = "P-values")
