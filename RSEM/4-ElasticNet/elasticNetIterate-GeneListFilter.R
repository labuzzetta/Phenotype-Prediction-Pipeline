directories <- list.dirs(recursive=FALSE)

library(glmnet)
library(AUC)

filter <- read.table("iLOO/Gene-List-without-Outliers.iLOO.txt", header=TRUE)

phenotypes <- c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1)

predictionList <- c()
labels <- c()

for(i in 1:length(directories)){
  
  setwd(directories[i])
  
  tableTRAIN <- read.table("EBSeq-TRAIN-selected.genes.txt", sep=" ", header=TRUE)
  tableTEST <- read.table("EBSeq-TEST-selected.genes.txt", sep=" ", header=TRUE)
  
  matrixTRAIN <- data.matrix(tableTRAIN)
  matrixTEST <- data.matrix(tableTEST)
  
  matrixTRAIN <- matrixTRAIN[,colnames(matrixTRAIN)%in%filter$x]
  matrixTEST <- matrixTEST[,colnames(matrixTEST)%in%filter$x]
  
  TRAINCV <- cv.glmnet(matrixTRAIN, phenotypes, nfolds=19, type.measure="deviance")
  
  prediction <- predict(TRAINCV, matrixTEST)
  
  predictionList <- c(predictionList, prediction)
  
  labels <- c(labels, 0, 1)
  
  setwd("..")
  
}

print(predictionList)

labels <- factor(labels)

print("Area Under ROC:")
print(auc(roc(predictionList, labels)))
print("Accuracy:")
print(auc(accuracy(predictionList, labels)))
print("Sensitivity:")
print(auc(sensitivity(predictionList, labels)))
print("Specificity:")
print(auc(specificity(predictionList, labels)))

plot(roc(predictionList, labels))
