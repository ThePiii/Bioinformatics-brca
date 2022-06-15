library(data.table)
library(factoextra)
library(FactoMineR)

mRNA <- fread("data/mRNA.csv", header=T) # 读的飞快
CNA <- fread("data/CNA.csv", header=T)

# 计算缺失值数量
sum(is.na(mRNA))
sum(is.na(CNA))

myfun <- function(x){
  cons <- mean(x,na.rm=TRUE)
  x <- nafill(x,'const',fill=cons)
  x
}


#CNA中有全部缺失的，先删掉
new.CNA <- subset(CNA,select=apply(CNA, 2, function(x) !all(is.na(x))))

# 到此CNA中没有全为na的列了，我检查过了

new.CNA <- cbind(new.CNA[,1],apply(new.CNA[,2:19827],2,myfun)) #检查过了，没有缺失值  495650


#mRNA也来一遍
new.mRNA <- subset(mRNA,select=apply(mRNA, 2, function(x) !all(is.na(x))))
#看看还有没有
sum(is.na(new.mRNA)) # 没有了，但是没想到这个时候mRNA的列数比CNA要少...


# PCA
CNA.pca <- PCA(new.CNA[,2:length(new.CNA)], ncp=500, graph=F)
CNA.var <- get_pca_var(CNA.pca)
res <- CNA.var$cor
head(CNA.var$cos2)
res <- CNA.pca$ind
View(res$dist)


CNA.pca <- prcomp(new.CNA[,2:length(new.CNA)], center=T, scale.=T)
eigs <- (CNA.pca$sdev)^2  # 各成分的方差
var_cum <- cumsum(eigs)/sum(eigs) # 方差比例
num_pc <- which.max(var_cum >= var_cum[500])
CNA.pcs <- CNA.pca$x[, 1:num_pc] # 降维后CNA 

mRNA.pca <- prcomp(new.mRNA[,2:length(new.mRNA)], center=T, scale.=T)
mRNA.pcs <- mRNA.pca$x[, 1:500] # 降维后mRNA
