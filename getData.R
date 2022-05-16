# pack <- "D:/Download/cgdsr_1.3.0.tar.gz"
# install.packages(pack, repos=NULL, type="source")

library(cgdsr)

# Create CGDS object
mycgds <- CGDS("http://www.cbioportal.org/")

# Test the CGDS endpoint URL using a few simple API tests
test(mycgds) 

# Get list of cancer studies at server
all_TCGA_studies <- getCancerStudies(mycgds)

# Get available case lists (collection of samples) for a given cancer study  
mycancerstudy <- getCancerStudies(mycgds)[50,1]
mycaselist <- getCaseLists(mycgds,mycancerstudy)[8,1]

# Get available genetic profiles
profiles <- getGeneticProfiles(mycgds,mycancerstudy)
mRNAid <-  getGeneticProfiles(mycgds,mycancerstudy)[14,1]
CNAid <- getGeneticProfiles(mycgds,mycancerstudy)[7,1]

# Get data slices for a specified list of genes, genetic profile and case list
# mRNA <- getProfileData(mycgds,c('BRCA1', 'BRCA2'),mRNAid,mycaselist) 
# CNA <- getProfileData(mycgds,c('BRCA1', 'BRCA2'),CNAid,mycaselist) 

# Get clinical data for the case list
clinicaldata <- getClinicalData(mycgds,mycaselist)  # 临床数据直接这里下应该没问题

geneslist <- read.table("data/geneslist.txt")

# 循环，500一轮读入所有基因的数据 
i <- 1
while(i <= length(geneslist[,1])){
  if(i == 1){
    mRNA <- matrix()
    CNA <- matrix()
  }
  k <- i+499
  temp1 <- getProfileData(mycgds, geneslist[i:k,], mRNAid, mycaselist)
  temp2 <- getProfileData(mycgds, geneslist[i:k,], CNAid, mycaselist)
  message(paste("done for genes -->",i," to ",k))
  if(dim(temp1)[1] != 0){
    mRNA <- cbind(mRNA, temp1)
    CNA <- cbind(CNA, temp2)
  }
  i <- k+1
}

mRNA <- mRNA[-1]
CNA <- CNA[-1]

# 保存本地 
write.csv(clinicaldata, "data/clincal.csv")
write.csv(mRNA, "data/mRNA.csv")
write.csv(CNA, "data/CNA.csv")

