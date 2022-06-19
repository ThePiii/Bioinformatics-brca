
# 第五题
#CNA data
library(igraph)
ADMM <- function(x){
  
  
  # 输入是要研究的数据表
  # 输出是一个list[theta1，data（邻接矩阵），
  # graph（这里的graph是挑选了某一个社区的数据进行可视化的）]
  
  S = cov(x)
  
  # lam rho 是自己调的，调到最后边数有100-200个
  lam = 0.2
  rho = 0.2
  
  theta1 = matrix(0,500,500)
  Z1 = matrix(0,500,500)
  U1 = matrix(0,500,500)
  r = matrix(1,500,500)
  s = matrix(1,500,500)
  epsilon.abs <- 1e-4
  epsilon.rel <- 1e-4
  epsilon.pri <- 0
  epsilon.dual <- 0
  i <- 0
  
  #迭代开始
  while (!(norm(r,'2')<=epsilon.pri & norm(s,'2')<=epsilon.dual)) {
    i <- i+1
    theta0 <- theta1
    Z0 <- Z1
    U0 <- U1
    rho <- 6.6*log(i)+0.2   #这里让他随迭代变大，在第20次迭代，大概达到20的样子；对偶系数没调
    
    
    decompo <- svd(rho*(Z0-U0)-S)
    Q <- decompo$u
    theta.hat <- (decompo$d+sqrt(decompo$d**2+4*rho))/2/rho
    
    theta1 <- Q%*%diag(theta.hat)%*%t(Q)
    Z1 <- (1-lam/(rho *  max(0, norm(theta1+U0,type='2')) ) )*(theta1+U0)
    U1 <- U0+(theta1-Z1)
    r <- theta1 - Z1
    s <- -rho*(Z1 - Z0)
    epsilon.pri <- sqrt(500)*epsilon.abs+epsilon.rel*max(norm(theta1,'2'),norm(Z1,'2'))
    epsilon.dual <- sqrt(500)*epsilon.abs+epsilon.rel*norm(rho*U1,'2')  
    if(i==22){break}
  }
  
  #迭代结束
  #到这里大概迭代了22次，我看网上说ADMM收敛比较慢，迭代个十几次就差不多行了
  
  
  
  
  # 看看留下150条边，阈值要取多少
  thres <- quantile(theta1[upper.tri(theta1)],c(1-150/250/499),na.rm=TRUE)
  
  
  #画图，画图之前记得把没有没有边的赋值成0
  
  
  
  data <- theta1
  colnames(data) <- colnames(x)
  diag(data) <- 0
  data[data <= thres ] <- 0
  data <- data[apply(data, 1, function(x) !all(x==0)),]
  data <- data[,apply(data, 2, function(x) !all(x==0))]
  network <- graph_from_adjacency_matrix(data,mode='undirected',weighted = TRUE)
  graph <- graph_from_adjacency_matrix(data[which(colnames(data) %in% edge.betweenness.community(network)[[1]]),which(colnames(data) %in% edge.betweenness.community(network)[[1]])],mode='undirected',weighted = TRUE)
  return(list(data,network,graph))
}

CNA.net <- ADMM(CNA.reduction)
mRNA.net <- ADMM(mRNA.reduction)


plot(mRNA.net[[3]])





# test
# 这里试了半天，发现在筛选变量的时候调大阈值会导致后面网络的节点个数增多；但不管怎么调，公共结点都是0 

# div <- function(x){
#   c <- quantile(x,c(0,0.25,0.5,0.75,1))
#   a <- c[2:5] - c[1:4]
#   ifelse(max(a)<1e6*min(a),0,1)  # 0 no prob   
# }
# 
# v <- apply(new.mRNA[,2:ncol(new.mRNA)],2,div)
# mRNA.red <- new.mRNA[,2:ncol(new.mRNA)][,!as.logical(v),with=FALSE]
# mRNA.v <- apply(mRNA.red, 2, var)
# mRNA.red <- mRNA.red %>%
#   select(order(-mRNA.v)[1:500])
# 
# mRNA.n <- ADMM(mRNA.red)
# intersection(CNA.graph,mRNA.n[[2]],keep.all.vertices = FALSE)
# 
# V(mRNA.n[[2]])$name %in% V(CNA.graph)$name
# 
# V(mRNA.n[[2]])


# (1)
CNA.graph <- CNA.net[[2]]
mRNA.graph <- mRNA.net[[2]]



CNA.graph.report <- cbind(get.edgelist(CNA.graph),E(CNA.graph)$weight)
colnames(CNA.graph.report) <- c('start','end','weight')
CNA.graph.report

mRNA.graph.report <- cbind(get.edgelist(mRNA.graph),E(mRNA.graph)$weight)
colnames(mRNA.graph.report) <- c('start','end','weight')
mRNA.graph.report

# (2)
# shared edges
CNA.edge <- E(CNA.graph)
mRNA.edge <- E(mRNA.graph)
c <- c()
for (i in length(CNA.edge)) {
  if(CNA.edge[i] %in% mRNA.edge){c <- c(c,i)}
}
common.edge <- intersection(CNA.graph,mRNA.graph,keep.all.vertices = FALSE)
length(common.edge)



CNA.vertex <- V(CNA.graph)
mRNA.vertex <- V(mRNA.graph)
sum(CNA.vertex$name %in% mRNA.vertex$name) 
CNA.vertex[CNA.vertex$name %in% mRNA.vertex$name]

# (3)
CNA.close <- closeness(CNA.graph)
CNA.deg <- degree(CNA.graph)
CNA.close[order(-CNA.close)[1:5]]
CNA.deg[order(-CNA.deg)[1:5]]

mRNA.close <- closeness(mRNA.graph)
mRNA.deg <- degree(mRNA.graph)
mRNA.close[order(-mRNA.close)[1:5]]
mRNA.deg[order(-mRNA.deg)[1:5]]

