# RCSM
RCSM(Recommended Connectivity-map Scoring Methods)
# Description
Connecting disease gene expression signatures to
compound-induced gene expression profiles is valuable for understanding
drug Mechanism of Action (MoA) and applying an existing therapeutic to
a new disease indication. RCSM is an R package for measuring those 
connectivities. RCSM implemented seven methods that are reported in 
published papers.
# Details of seven methods
KSScore(refMatrix, queryUp, queryDown, permuteNum = 10000,
  pAdjMethod = "BH", mcCore = 1)

XCosScore(refMatrix, query, topN = 500, permuteNum = 10000,
  pAdjMethod = "BH", mcCore = 1)
  
XSumScore(refMatrix, queryUp, queryDown, topN = 500,
  permuteNum = 10000, pAdjMethod = "BH", mcCore = 1)
  
GSEAweight0Score(refMatrix, queryUp, queryDown, permuteNum = 10000,
  pAdjMethod = "BH", mcCore = 1)
  
GSEAweight1Score(refMatrix, queryUp, queryDown, permuteNum = 10000,
  pAdjMethod = "BH", mcCore = 1)

GSEAweight2Score(refMatrix, queryUp, queryDown, permuteNum = 10000,
  pAdjMethod = "BH", mcCore = 1)
  
ZhangScore(refMatrix, queryUp, queryDown, permuteNum = 10000,
  pAdjMethod = "BH", mcCore = 1)

# How to install this package?
install.packages("devtools")

library(devtools)

install_github("Jasonlinchina/RCSM")

library(RCSM)

# If there is any bug or suggestion, please contact kequanlin@gmail.com

# FAQ
1, If you encounter the error saying "readRDS(dest) XXX" when installing
the RCSM, you could first remove the old version of RCSM with "remove.package(RCSM)",
and then restart the R or Rstudio.
