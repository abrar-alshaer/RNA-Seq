setwd("D:/OneDrive - University of North Carolina at Chapel Hill/scripts/Smith Lab/Wei Sha Analysis/Abrar Analysis/Illumia Data/RERUN/Galgal5Run/miRNA DEXSeq+PCA_Kmeans/")
#setwd("C:/Users/Abrar/OneDrive - University of North Carolina at Chapel Hill/scripts/Smith Lab/Wei Sha Analysis/Abrar Analysis/Illumia Data/RERUN/Galgal5Run/miRNA DEXSeq+PCA_Kmeans/")
library("FactoMineR")
library(ggplot2)
library("factoextra")
library(NbClust)
library("corrplot")
library(clustertend)
library(randomcoloR)
library(zoom)
library(ggfortify)
library(cluster)
rm(list=ls())

#DEXSeq_pval_less_0.1_miRNAs_PCA_file_ALL_FC
#DEXSeq_pval_less_0.1_miRNAs_PCA_file_PosFC
#DEXSeq_pval_less_0.1_miRNAs_PCA_file_NegFC
exonData <- read.csv("DEXSeq_pval0.1_miRexons_PCA_ALL_FC.csv", header = TRUE, fill = TRUE) #row.names = 1
exonData <- exonData[-2,] #to remove 1 row from the dataframe
nams = exonData[,12] #change to 12 for miR ID to be row names
rownames(exonData) = make.names(nams, unique=TRUE) #to avoid duplicate names error

#Defining active variables for PCA
exonData.active <- exonData[,c(2:5,7:8)]

#running PCA
res.PCA <- PCA(exonData.active, graph=FALSE)

#getting eigenvalues
eig.val <- get_eigenvalue(res.PCA)
eig.val

#initialize a high resolution graph
png("scree plot miRNAs exclude3064_highRes.png", units="in", width=7, height=7, res=600) 

#Scree plot
fviz_eig(res.PCA, addlabels = TRUE)

#release PNG file graph into directory
dev.off()

#provides a list of matrices containing all the results for the active variables 
#(coordinates, correlation between variables and axes, squared cosine and contributions)
var <- get_pca_var(res.PCA)
var
#coordinates
var$coord
#Cos2:Quality on the factor map
var$cos2
#Contributions to the principal components
var$contrib

#Visualize correlation between the variables (columns) and principal componenet (PC)
fviz_pca_var(res.PCA, col.var = "black", repel = TRUE)

#visualize cos2
corrplot(var$cos2, is.corr = FALSE)
#fviz_cos2(res.PCA, choice = "var", axes = 1:2) #bar plot of above

#color by cos2 values: quality on the factor map
fviz_pca_var(res.PCA, col.var = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) #repel avoids overlapping labels

#correlation plot of variable contributions (variables/columns with contributions in last dimensions should be removed)
corrplot(var$contrib, is.corr = FALSE)
#bar graph of contributions to PC1
fviz_contrib(res.PCA, choice = "var", axes = 1)
#bar graph of contributions to PC2
fviz_contrib(res.PCA, choice = "var", axes = 2)
#bar graph of contributions to PC1-2
fviz_contrib(res.PCA, choice = "var", axes = 1:2)
#contributions on correlation plot
fviz_pca_var(res.PCA, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) #repel avoids overlapping labels

#Coloring variables (columns) by Kmean clustering
#determine how many groups to use for kmeans algorithm
exonData.keggs <- exonData[,c(2,3,4,5,6,7,8)]
names = exonData[,9]
rownames(exonData.keggs) = make.names(names, unique=TRUE)
fviz_nbclust(exonData.keggs, kmeans, method = "wss") + geom_vline(xintercept = 4, linetype = 2)
#optimal kmeans clusters = 4
res.km <- kmeans(var$coord, centers = 4, nstart = 25)
grp <- as.factor(res.km$cluster)
fviz_pca_var(res.PCA, col.var = grp, palette = c("purple", "blue", "red", "green"), legend.title = "Cluster", repel = TRUE)

#Coloring variables (columns) by KEGG Pathway groups ####DOES NOT WORK####
#grp1 <- as.factor(unique(exonData$KEGG_PATHWAY))
#colors <- distinctColorPalette(k = 23) #make 23 random colors
#fviz_pca_var(res.PCA, col.var = grp1, palette = colors, legend.title = "Cluster", repel = TRUE)

#Statistical corr & pvalue of dimension descriptions
res.desc <- dimdesc(res.PCA, axes = c(1,2), proba = 0.05)
res.desc$Dim.1
res.desc$Dim.2

ind <- get_pca_ind(res.PCA)
ind

head(ind$coord)
head(ind$cos2)
head(ind$contrib)

dev.off()
#All individuals
fviz_pca_ind(res.PCA, geom = "point")
#All cos2
fviz_pca_ind(res.PCA, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),  geom = "point")
#selecting top 50 according to the cos2
fviz_pca_ind(res.PCA, geom = "point", select.ind = list(cos2 = 50))
#selecting top 50 according to the contribution
fviz_pca_ind(res.PCA, geom = "point", select.ind = list(contrib = 50))
#selecting top 10 according to the contribution
fviz_pca_ind(res.PCA, repel = TRUE, select.ind = list(contrib = 10)) #label by pathway????

fviz_cos2(res.PCA, choice = "ind", select.ind = list(contrib = 10)) #select didn't work..?

#initialize a high resolution graph
png("PCA allFC miRNAs exclude3064_highRes.png", units="in", width=9, height=7, res=600) 

#PCA GRAPH WORKS
ind.p <- fviz_pca_ind(res.PCA, geom.ind = c("point", "text"), labelsize = 3, pointshape = 20, pointsize = 3, repel = TRUE, col.ind = exonData$miRBase.ID, mean.point = FALSE, addEllipses = FALSE, legend.title = "miRNA IDs")
ggpubr::ggpar(ind.p, title = "Principal Component Analysis", subtitle = "All miRNAs excluding gga-miR-3064 Exon 1", xlab = "PC1 (93.8%)", ylab = "PC2 (3.3%)", legend.title = "miRNA IDs", legend.position = "top", ggtheme = theme_classic()) #palette = "simpsons

dev.off()

#initialize a high resolution graph
png("2-PCA allFC miRNAs exclude3064_highRes.png", units="in", width=9, height=7, res=600) 

#PCA GRAPH WORKS Log2 FC coloring
ind.p <- fviz_pca_ind(res.PCA, geom.ind = c("point", "text"), labelsize = 3, pointshape = 20, pointsize = 3, repel = TRUE, col.ind = exonData$log2fold_SALINE_ETOH, mean.point = FALSE, addEllipses = FALSE, legend.title = "Log2 Fold Change")
ggpubr::ggpar(ind.p, title = "Principal Component Analysis", subtitle = "All miRNAs excluding gga-miR-3064 Exon 1", xlab = "PC1 (93.8%)", ylab = "PC2 (3.3%)", legend.title = "Log2 Fold Change", legend.position = "top", ggtheme = theme_classic()) #simpsons

dev.off()

#PCA GRAPH WORKS
scale_exons <- scale(exonData.active)
myPCA <- princomp(scale_exons)
summary(myPCA)
plot(myPCA$scores[,1], myPCA$scores[,2])
autoplot(prcomp(exonData.active))
#Without labels 
autoplot(prcomp(exonData.active, scale. = TRUE), data = exonData, colour = 'Fold.Change', label = TRUE, label.repel=TRUE, label.size = 3)
#With labels
autoplot(prcomp(exonData.active, scale. = TRUE), data = exonData, colour = 'miRBase.ID', label = TRUE, label.repel=TRUE, label.size = 3, title = "miRNA Upregulated in ETOH")#label.repel=TRUE

#plotting a biplot
fviz_pca_biplot(res.PCA, repel = TRUE, select.ind = list(contrib = 20), col.var = "#2E9FDF", col.ind = "#696969")

#HCPC 
res.hcpc <- HCPC(res.PCA, graph=FALSE)
#step below takes a long time - may skip
fviz_dend(res.hcpc, cex = 0.7, palette = "jco", rect = TRUE, rect_fill = TRUE, rect_border = "jco", show_labels = TRUE)

#Creating a factor map
fviz_cluster(res.hcpc, show.clust.cent = TRUE, palette = "jco", ggtheme = theme_minimal(), main = "Factor Map", geom = "point")
fviz_cluster(res.hcpc, show.clust.cent = TRUE, palette = "jco", ggtheme = theme_minimal(), main = "Factor Map", geom = c("point", "text"), repel = TRUE)
#DOESN'T WORK# fviz_dend(res.hcpc, show_labels = FALSE)

plot(res.hcpc, choice = "3D.map")
write.csv(res.hcpc$data.clust, "PCA__HCPC_DEXSeq_NormCounts_miRNA.csv")

#K-means clustering
#xintercept specifies where the line is drawn #USE ELBOW METHOD FOR K-MEANS!#
png("k-means clusters allFC miRNAs exclude 3064_highRes.png", units="in", width=9, height=7, res=600) 

fviz_nbclust(exonData.active, kmeans, method = "wss") + geom_vline(xintercept = 3, linetype = 2) + labs(subtitle = "Elbow Method - Kmeans Clustering") #k = 3

dev.off()

fviz_nbclust(exonData.active, kmeans, method = "silhouette") + labs(subtitle = "Silhouette Method") #k=2
fviz_nbclust(exonData.active, kmeans, method = "gap_stat", nboot = 500) + labs(subtitle = "Gap Statistic Method - K-means Clustering, nboot = 500") #k=4 with nboot=50, k=5 with nboot=500

png("k-means allFC miRNAs exclude3064_highRes.png", units="in", width=9, height=7, res=600)
#optimal kmeans clusters = 3
km <- kmeans(exonData.active, 3, nstart = 50)
fviz_cluster(km, data = exonData.active, ellipse.type = "euclid", star.plot = TRUE, repel = TRUE, ggtheme = theme_minimal(), title = "K-means Clustering", subtitle = "All miRNAs excluding gga-mir-3064 exon 1")
dev.off()

autoplot(km, data = exonData.active, label = FALSE, label.repel = TRUE, label.size = 3) #label = TRUE, label.size = 3

#assessing cluster tendency statistics (measures probability that a given dataset is generated by a uniform distribution - tests for spatial randomness in the data)
hopkins(exonData.active, n = nrow(exonData.active)-1)
#h = 0.2104383 - H=0.5 is the threshold, the closer to 0.5 and above the higher probability that the clustering is due to random chance, and that your dataset is uniformly clustered

#compute hopkins stat for random dataset
random <- data.frame(replicate(7,sample(5:4000,22,rep=TRUE)))
hopkins(random, n = nrow(random)-1)
#h=0.4994303 - expected to be near 0.5, since this is a random dataset

#Visual assessment of cluster tendency - VAT (visually seeing the statistics)
png("miRNA RNA-Seq Exon Data exclude3064-VAT_highRes2.png", units="in", width=9, height=7, res=600)
fviz_dist(dist(exonData.active), show_labels = TRUE) + labs(title = "miRNA RNA-Seq Exon Data")
dev.off()

png("Random Data Matrix-VAT_highRes.png", units="in", width=9, height=7, res=600)
fviz_dist(dist(random), show_labels = TRUE) + labs(title = "Random Data Matrix") #RANDOM DATASET
dev.off()

########################Output isn't great
#CLARA clustering
#optimal clusters = 2
fviz_nbclust(exonData.active, clara, method = "wss") + geom_vline(xintercept = 3, linetype = 2) #+ theme_classic()
clara.res <- clara(exonData.active,3,samples = 100, pamLike = TRUE)
autoplot(clara.res, frame = TRUE)
########################Output isn't great

#Generating bar graphs of the pathways in each HCPC cluster
up <- read.csv("PCA_Exons_HCPC_UP.csv", header = TRUE, fill = TRUE)
clust <- down[down$clust == 5,]
tab <- table(clust$KEGG_PATHWAY)
op <- par(mar = c(10,4,4,2) + 0.1) #TO EXTEND MARGINS
barplot(tab[order(tab, decreasing = T)], cex.names=0.6, las=2, main = "Cluster 5 - Increasing in ETOH", ylab = "Number of Exons")
dotchart(tab) #cool dot plot
#HCPC_DW-Exon_Clust3
