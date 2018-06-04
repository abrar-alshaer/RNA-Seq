setwd("D:/OneDrive - University of North Carolina at Chapel Hill/scripts/Smith Lab/Wei Sha Analysis/Abrar Analysis/Illumia Data/RERUN/Galgal5Run/DEXSeq PCA_Kmeans_Code+Output/")
library("FactoMineR")
library(ggplot2)
library("factoextra")
library("corrplot")
library(randomcoloR)
library(zoom)
library(ggfortify)
library(cluster)
rm(list=ls())

#DAVID_TEST_DEXSeqCounts_neg-posFC_padj0.1_unstr-B_rm0_filterKEGGs_1E1P #NOT NORMALIZED FILE!!!
#DAVID_TEST_DEXSeqResults_neg-posFC_padj0.1_unstr-B_rmNA_1E1P_PCA_FINAL-UP.csv #FOR UPREGULATED EXONS, DW for DOWNREGULATED
exonData <- read.csv("DAVID_TEST_DEXSeqResults_neg-posFC_padj0.1_unstr-B_rmNA_1E1P_PCA_FINAL-UP.csv", header = TRUE, fill = TRUE) #row.names = 1
exonData = exonData[!duplicated(exonData$Exons),] #FILTER AND REMOVE DUPLICATES
nams = exonData[,1] #place 1 to label by exon IDs and 11 to label by pathway
rownames(exonData) = make.names(nams, unique=TRUE) #to avoid duplicate names error

#Defining active variables for PCA
#exonData[, 2:5] <- log(exonData[2:5], 2)
exonData.active <- exonData[,c(2:5, 7:8)]

#running PCA
res.PCA <- PCA(exonData.active, graph=FALSE) #scale.unit = TRUE IS THE DEFAULT!!!

#getting eigenvalues
eig.val <- get_eigenvalue(res.PCA)
eig.val


#initialize a high resolution graph
png("PCA-DW_ScreePlot_highRes.png", units="in", width=7, height=7, res=600) 

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

#initialize a high resolution graph
png("corr plot downregulated exons_highRes.png", units="in", width=7, height=7, res=600) 

#color by cos2 values: quality on the factor map
fviz_pca_var(res.PCA, col.var = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) #repel avoids overlapping labels

#release PNG file graph into directory
dev.off() 

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
names = exonData[,11]
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
#plot top 50 contributions
fviz_contrib(res.PCA, choice = "ind", axes = 1:2, top = 50)

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
png("PCA-AllFC_col_Unique50-highRes.png", units="in", width=11, height=9, res=600) 

#PCA GRAPH WORKS (USE THIS ONE TO SELECT TOP INDIVIDUALS FOR PCA)
#colors <- distinctColorPalette(k = 21) #make 23 random colors
####NO COLORING####ind.p <- fviz_pca_ind(res.PCA, geom.ind = c("point", "text"), labelsize = 3, pointshape = 20, pointsize = 3, repel = TRUE, select.ind = list(contrib = 50), mean.point = FALSE, addEllipses = FALSE) #NO COLORING
ind.p <- fviz_pca_ind(res.PCA, geom.ind = c("point", "text"), labelsize = 3, pointshape = 20, pointsize = 3, repel = TRUE, select.ind = list(contrib = 50), col.ind = exonData$KEGG_PATHWAY, mean.point = FALSE, addEllipses = FALSE, legend.title = "KEGG PATHWAYS") # palette = "jco",  # + coord_cartesian(xlim = c(-250000, 50), ylim=c(-9000,9000))
ggpubr::ggpar(ind.p, title = "Principal Component Analysis - Downregulated Exons", subtitle = "Top 50 exons contributing to variance of principal components", xlab = "PC1 (69.2%)", ylab = "PC2 (30%)", legend.title = "KEGG Pathways", legend.position = "top", ggtheme = theme_classic(), palette = "simpsons") #simpsons

#color by log2 FC
ind.p <- fviz_pca_ind(res.PCA, geom.ind = c("point", "text"), labelsize = 3, pointshape = 20, pointsize = 3, repel = TRUE, select.ind = list(contrib = 50), col.ind = exonData$log2fold_SALINE_ETOH, mean.point = FALSE, addEllipses = FALSE, legend.title = "Log2 Fold Change") # palette = "jco",  # + coord_cartesian(xlim = c(-250000, 50), ylim=c(-9000,9000))
ggpubr::ggpar(ind.p, title = "Principal Component Analysis - All Exons", subtitle = "Top 50 exons contributing to variance of principal components", xlab = "PC1 (69.6%)", ylab = "PC2 (29.5%)", legend.title = "Log2 Fold Change", legend.position = "top", ggtheme = theme_classic()) #simpsons

#release PNG file graph into directory
dev.off() 

#PCA GRAPH WORKS
scale_exons <- scale(exonData.active)
myPCA <- princomp(scale_exons)
summary(myPCA)
plot(myPCA$scores[,1], myPCA$scores[,2])
autoplot(prcomp(exonData.active), label = TRUE)
#Without labels
autoplot(prcomp(exonData.active, scale. = TRUE), data = exonData, colour = 'KEGG_PATHWAY')
#With labels
autoplot(prcomp(exonData.active, scale. = TRUE), data = exonData, colour = 'KEGG_PATHWAY', label = TRUE, label.size = 3)#label.repel=TRUE

#plotting a biplot
fviz_pca_biplot(res.PCA, repel = TRUE, select.ind = list(contrib = 50), col.var = "#2E9FDF", col.ind = "#696969")

#HCPC 
res.hcpc <- HCPC(res.PCA, graph=FALSE)
#step below takes a long time - may skip
fviz_dend(res.hcpc, cex = 0.7, palette = "jco", rect = TRUE, rect_fill = TRUE, rect_border = "jco", show_labels = FALSE)

#Creating a factor map
fviz_cluster(res.hcpc, show.clust.cent = TRUE, palette = "jco", ggtheme = theme_minimal(), main = "Factor Map", geom = "point")
fviz_cluster(res.hcpc, show.clust.cent = TRUE, palette = "jco", ggtheme = theme_minimal(), main = "Factor Map", geom = c("point", "text"), repel = TRUE)
#DOESN'T WORK# fviz_dend(res.hcpc, show_labels = FALSE)

#initialize a high resolution graph
png("HCPC Upregulated exons-highRes.png", units="in", width=9, height=7, res=600)

plot(res.hcpc, choice = "3D.map", ind.names = FALSE)
write.csv(res.hcpc$data.clust, "PCA_Exons_Counts_HCPC_UP_exonIDs.csv")

#release PNG file graph into directory
dev.off() 

#K-means clustering
fviz_nbclust(exonData.active, kmeans, method = "wss") + geom_vline(xintercept = 4, linetype = 2)
#optimal kmeans clusters = 4
km <- kmeans(exonData.active, 3, nstart = 25) #could also use exonData.active (same thing)
autoplot(km, data = exonData.active, label = TRUE, label.size = 3) #label = TRUE, label.size = 3

########################Output isn't great
#CLARA clustering
#optimal clusters = 2
fviz_nbclust(exonData.active, clara, method = "silhouette") + theme_classic()
clara.res <- clara(exonData.active,2,samples = 100, pamLike = TRUE)
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
