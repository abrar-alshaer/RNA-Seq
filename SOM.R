setwd("C:/Users/Abrar/Dropbox/UNC_OneDrive/scripts/Smith Lab/Wei Sha Analysis/Abrar Analysis/Illumia Data/RERUN/Galgal5Run/DEXSeq PCA_Kmeans_Code+Output/")
#Code was adapted from: https://www.r-bloggers.com/self-organising-maps-for-customer-segmentation-using-r/ 
library("FactoMineR")
library(ggplot2)
library("factoextra")
library(cluster)
library(kohonen)
rm(list=ls())

exonData <- read.csv("DAVID_TEST_DEXSeqResults_neg-posFC_padj0.1_unstr-B_rmNA_1E1P_PCA_FINAL-UP.csv", header = TRUE, fill = TRUE) #row.names = 1
exonData = exonData[!duplicated(exonData$Exons),] #FILTER AND REMOVE DUPLICATES
nams = exonData[,1] #place 1 to label by exon IDs and 11 to label by pathway
rownames(exonData) = make.names(nams, unique=TRUE) #to avoid duplicate names error

exonData.active <- exonData[,c(2:5, 7:8)]

#Self Organizing Maps (SOMs)
## training grid prior to training the SOM. Hexagonal and Rectangular are possible
#aim for 5-10 samples per node 
#we are making a 15*16 map which is 240 nodes with about 10 samples per node given 2420 nodes
som_grid <- somgrid(xdim = 15, ydim=10, topo="hexagonal") #the x*y should be roughly equal to your n observations

# Train the SOM, options for the number of iterations,
# the learning rates, and the neighbourhood are available
som_model <- som(as.matrix(exonData.active), 
                 grid=som_grid, rlen=15000, #rlen = epochs (iterations that will learn the data and try to fit the observations in a model of the attributes (predictors))
                 alpha=c(0.05,0.01), #learning rate = how much displacement between each epoch will the neurons have
                 keep.data = TRUE)

#As the SOM training iterations progress, the distance from each node's weights to the samples 
#represented by that node is reduced. #Ideally, this distance should reach a minimum plateau. 
#This plot option shows the progress over time. 
#If the curve is continually decreasing, more iterations are required. 
png("TrainingSOM_DUMMY.png", units="in", width=9, height=7, res=600)
plot(som_model, type="changes", main = "Training Progress - Downregulated Exons")
dev.off()

#allows us to visualise the count of how many samples are mapped to each node on the map
#This metric can be used as a measure of map quality - ideally the sample distribution is relatively uniform. 
png("SOM_CountsPlot_AllExons.png", units="in", width=9, height=7, res=600)
plot(som_model, type="count")
dev.off()

#Often referred to as the "U-Matrix", this visualisation is of the distance between each node and its neighbours
#areas of low neighbour distance indicate groups of nodes that are similar
#Areas with large distances indicate the nodes are more dissimilar & indicate natural boundaries btwn node clusters
plot(som_model, type="dist.neighbours")

plot(som_model, type="codes")

plot(som_model, type = "property", property = som_model$codes[,1], main=names(som_model$data)[1])

df_c <- data.frame(som_model$codes)

png("SOM_predictor_heatmaps_All_Exons.png", units="in", width=9, height=7, res=600)
par(mfrow=c(3,2))
plot(som_model, type = "property", property = df_c[,1], main = "ETOH Counts 1")
plot(som_model, type = "property", property = df_c[,2], main = "ETOH Counts 2")
plot(som_model, type = "property", property = df_c[,3], main = "Saline Counts 1")
plot(som_model, type = "property", property = df_c[,4], main = "Saline Counts 2")
plot(som_model, type = "property", property = df_c[,5], main = "ETOH Exon Usage Coefficient")
plot(som_model, type = "property", property = df_c[,6], main = "Saline Exon Usage Coefficient")
dev.off()

#hierarchal clustering
#get optimal # of clusters
png("SOM_optimal_clusters_All_Exons3.png", units="in", width=9, height=7, res=600)
fviz_nbclust(df_c, kmeans, method = "wss") + geom_vline(xintercept = 3, linetype = 2) + labs(subtitle = "Elbow Method - Kmeans Clustering on SOM") #k = 3
fviz_nbclust(df_c, kmeans, method = "silhouette") + labs(subtitle = "Silhouette Method on SOM") #k=2
fviz_nbclust(df_c, kmeans, method = "gap_stat", nboot = 500) + labs(subtitle = "Gap Statistic Method - K-means Clustering on SOM, nboot = 500") #k=4
dev.off()

## use hierarchical clustering to cluster the codebook vectors
som_cluster <- cutree(hclust(dist(df_c)), 2)
# plot the results:
png("SOM_hierarchal_clustering_k2_All_Exons.png", units="in", width=9, height=7, res=600)
plot(som_model, type="mapping", main = "Clusters of All Exons \n Self Organizing Map")
add.cluster.boundaries(som_model, som_cluster)
dev.off()
