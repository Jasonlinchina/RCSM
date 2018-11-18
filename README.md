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
cmapScore(refMatrix, queryUp, queryDown, permuteNum = 10000,
  pAdjMethod = "BH", mcCore = 1)

eXtremeCosScore(refMatrix, query, topN = 500, permuteNum = 10000,
  pAdjMethod = "BH", mcCore = 1)
  
eXtremeSumScore(refMatrix, queryUp, queryDown, topN = 500,
  permuteNum = 10000, pAdjMethod = "BH", mcCore = 1)
  
GSEAweight0Score(refMatrix, queryUp, queryDown, permuteNum = 10000,
  pAdjMethod = "BH", mcCore = 1)
  
GSEAweight1Score(refMatrix, queryUp, queryDown, permuteNum = 10000,
  pAdjMethod = "BH", mcCore = 1)

GSEAweight2Score(refMatrix, queryUp, queryDown, permuteNum = 10000,
  pAdjMethod = "BH", mcCore = 1)
  
zhangScore(refMatrix, queryUp, queryDown, permuteNum = 10000,
  pAdjMethod = "BH", mcCore = 1)

# How to install this package?
install.packages("devtools")

library(devtools)

install_github("Jasonlinchina/RCSM")

library(RCSM)
