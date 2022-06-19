library(PMA)
out <- MultiCCA(list(CNA.reduction,mRNA.reduction))

CNA.nonzero <- colnames(subset(CNA.reduction,select=out$ws[[1]] != 0.0))
CNA.cano.covar <- out$ws[[1]][out$ws[[1]] != 0.0]

CNA.cano <- rbind(gene=CNA.nonzero,covar=CNA.cano.covar)
CNA.cano


mRNA.nonzero <- colnames(subset(mRNA.reduction,select=out$ws[[2]] != 0.0))
mRNA.cano.covar <- out$ws[[2]][out$ws[[2]] != 0.0]

mRNA.cano <- rbind(gene=mRNA.nonzero,covar=mRNA.cano.covar)
mRNA.cano