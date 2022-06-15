# http://compbio.cs.toronto.edu/SNF/SNF/Software.html

# Please conduct spectral clustering on each omics data, and the SNF
# approach on multi-omics data. Summary the similarity matrices with the 
# heatmap similar to the following figure, and the clustering results (including 
# the number of clusters and size of each cluster).


library(SNFtool)

Dist1 <- (dist2(as.matrix(CNA.reduction), as.matrix(CNA.reduction)))^0.5
Dist2 <- (dist2(as.matrix(mRNA.reduction), as.matrix(mRNA.reduction)))^0.5

# construct similarity graphs
K = 20  #number of neighbors, usually (10~30)
alpha = 0.5   #hyperparameter, usually (0.3~0.8)
T = 20  #Number of Iterations, usually (10~50)

W1 <- affinityMatrix(Dist1, K, alpha)
W2 <- affinityMatrix(Dist2, K, alpha)

# Next, we fuse all the graphs, then the overall matrix can be computed by
W = SNF(list(W1,W2), K, T)


# With this unified graph W of size n x n, 
# you can do either spectral clustering or Kernel NMF. 
# If you need help with further clustering, please let us know. 

# You can display clusters in the data by the following function
# where C is the number of clusters.
C = 3

# You can get cluster labels for each data point by spectral clustering
labels = spectralClustering(W, C)

# On each omics data
labels_CNA <- spectralClustering(W1, C)
labels_mRNA <- spectralClustering(W2, C)

truelabel = c(matrix(1,1100/2,1),matrix(2,1100/2,1))
M_label=cbind(labels,truelabel)
colnames(M_label)=c("spectralClustering","TrueLabel")
                                                                                                                      colors=c("orange","cyan")))

displayClustersWithHeatmap(W1, labels_CNA)
displayClustersWithHeatmap(W2, labels_mRNA)

displayClustersWithHeatmap(W, labels)
displayClusters(W, labels)
