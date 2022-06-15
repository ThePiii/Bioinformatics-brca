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


# 计算方差降维
CNA.var <- apply(new.CNA[,2:length(new.CNA)], 2, var)
CNA.reduction <- new.CNA %>%
  select(order(-CNA.var[1:500]))  # order默认升序，需要加个负号

mRNA.var <- apply(new.mRNA[,2:length(new.mRNA)], 2, var)
mRNA.reduction <- new.mRNA %>%
  select(order(-mRNA.var)[1:500])
