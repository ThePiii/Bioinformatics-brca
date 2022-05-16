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
mRNA <- getProfileData(mycgds,c('BRCA1', 'BRCA2'),mRNAid,mycaselist)
CNA <- getProfileData(mycgds,c('BRCA1', 'BRCA2'),CNAid,mycaselist)

# Get clinical data for the case list
clinicaldata <- getClinicalData(mycgds,mycaselist)

# save data
library(writexl)
write_xlsx(mRNA, 'data/mRNA.xlsx')
write_xlsx(CNA, 'data/CNA.xlsx')
write_xlsx(clinicaldata, 'data/clinical.xlsx')
