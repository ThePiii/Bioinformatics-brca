library(glmgraph)
library(dplyr)

clinical <- read.csv("data/clincal.csv")

response <- clinical %>%
  select(X, HER2_CENT17_RATIO) %>%
  left_join(cbind(mRNA[,1], mRNA.reduction), by=c("X" = "V1")) %>%
  filter(!(is.na(HER2_CENT17_RATIO)))

# glmgraph(X, Y, L, family=c("gaussian","binomial"), penalty=c("MCP","lasso") ,
#          mcpapproach=c("mmcd", "adaptive", "original"),gamma=8,
#          lambda1,nlambda1=100,lambda2=c(0, 0.01 * 2^(0:7)),eps=1e-3,max.iter=2000,
#          dfmax=round(ncol(X)/2),penalty.factor=rep(1,ncol(X)),standardize=TRUE,warn=FALSE,...)

X <- response[, 3:length(response)]
Y <- response$HER2_CENT17_RATIO

L <- diag(500)

glmgraph(X = X, Y = Y, L = L, family = "gaussian", penalty = "lasso")
