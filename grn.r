setwd(dir = "/Users/farihatanjin/Desktop/lab/REIN/data")
library(GENIE3)

grn <- function(condition_name) {
  
  #read single cell expression data into a matrix
  filename <- paste("condition_data/scRecover+scImpute_",condition_name,"_condition.csv", sep="")
  exprMatr <- read.csv(filename, sep=",", row.names = 1)
  exprMatr <- as.matrix(exprMatr) 
  
  set.seed(123) # For reproducibility of results
  
  #set selected genes as regulators 
  reg <- read.csv("genelist.csv", sep=",", row.names = 1)
  regulators <- c(row.names(reg))
  
  weightMat <- GENIE3(exprMatr, regulators=regulators)
  weightMat <- as.matrix(weightMat)
  
  #threshold is the weight of the regulatory link between Pou5f1 and Nanog  
  threshold <- weightMat["Pou5f1","Nanog"]
  
  #list of gene interactions with its weight 
  linkList <- getLinkList(weightMat, threshold=threshold)
  
  #save list to csv
  filename <- paste("/Users/farihatanjin/Desktop/lab/REIN/results_genie3/interactions_",condition_name,".csv",sep="")
  write.csv(linkList,filename, row.names=FALSE)
}

#list of conditions
data <- c("2iLif","CiTC","LB")

for (x in data) {
  grn(x)
}




