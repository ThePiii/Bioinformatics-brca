library(data.table)

mRNA <- fread("data/mRNA.csv", header=T) # 读的飞快
CNA <- fread("data/CNA.csv", header=T)

# 计算缺失值数量
sum(is.na(mRNA))
sum(is.na(CNA))

