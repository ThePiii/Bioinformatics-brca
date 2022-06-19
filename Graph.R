library(glmgraph)
library(dplyr)

clinical <- read.csv("data/clincal.csv")

brca_related <- select(mRNA, V1, BRCA1, BRCA2, PALB2, CHEK2, CDH1, PTEN, STK11, TP53, 
                       ATM, BARD1, BRIP1, CASP8, CTLA4, CYP19A1, FGFR2, H19, MAP3K1, NBN, RAD51, TERT) 

response <- clinical %>%
  select(`X`, HER2_CENT17_RATIO) %>%
  left_join(brca_related, by=c("X" = "V1")) %>%
  filter(!(is.na(HER2_CENT17_RATIO)))


# glmgraph(X, Y, L, family=c("gaussian","binomial"), penalty=c("MCP","lasso") ,
#          mcpapproach=c("mmcd", "adaptive", "original"),gamma=8,
#          lambda1,nlambda1=100,lambda2=c(0, 0.01 * 2^(0:7)),eps=1e-3,max.iter=2000,
#          dfmax=round(ncol(X)/2),penalty.factor=rep(1,ncol(X)),standardize=TRUE,warn=FALSE,...)

X <- as.matrix(response[, 3:length(response)])  # 234*520
Y <- as.matrix(response$HER2_CENT17_RATIO)  # 234*1


# 构建拉普拉斯矩阵，主要围绕BRCA1和BRCA2瞎整

p <- 20  # X的个数
p1 <- 8  # 强相关
A <- matrix(rep(0, p*p), p, p)
A[1:2,1:p] <- 1
A[1:p, 1:2] <- 1
A[1:p1, 1:p1] <- 1
diag(A) <- 0
D <- apply(A,1,sum)
L <- - A
diag(L) <- D  # L = D - A

cv.glmgraph(X, Y, L, type.measure = "mse", nfolds=5, trace=TRUE)  # 用交叉验证来进行参数选择

obj <- glmgraph(X = X, Y = Y, L = L, family = "gaussian", penalty = "lasso", lambda2 = 1.28, lambda1 = 0.094498, warn=TRUE)

beta <- as.data.frame(obj$betas)

# ---- 考虑500个 ----

brca_all <- select(mRNA, V1, BRCA1, BRCA2, PALB2, CHEK2, CDH1, PTEN, STK11, TP53, 
                   ATM, BARD1, BRIP1, CASP8, CTLA4, CYP19A1, H19, MAP3K1, NBN, RAD51)  %>%
  cbind(mRNA.reduction) %>%
  unique()

response1 <- clinical %>%
  select(`X`, HER2_CENT17_RATIO) %>%
  left_join(brca_all, by=c("X" = "V1")) %>%
  filter(!(is.na(HER2_CENT17_RATIO)))

X <- as.matrix(response1[, 3:length(response1)])  # 234*520
Y <- as.matrix(response1$HER2_CENT17_RATIO)  # 234*1


# 构建拉普拉斯矩阵，主要围绕BRCA1和BRCA2瞎整

n <- 518
p <- 20  # X的个数
p1 <- 8  # 强相关
A <- matrix(rep(0,n*n), n, n)
A[1:2,1:n] <- 1
A[1:n, 1:2] <- 1
A[1:p, 1:p] <- 1
diag(A) <- 0
D <- apply(A,1,sum)
L <- - A
diag(L) <- D  # L = D - A

cv.glmgraph(X, Y, L, type.measure = "mse", nfolds=5, trace=TRUE)  # 用交叉验证来进行参数选择

obj <- glmgraph(X = X, Y = Y, L = L, family = "gaussian", penalty = "lasso", lambda2 = 1.28 , lambda1 = 0.031257, warn=TRUE)

beta <- as.data.frame(obj$betas)  %>%
  filter(beta$X0.0313 != 0)
