t.test.select.randomForest <- function(expTable, ctrlTable, p_value) {
  
  #expTable <- read.table(exp_path, header = TRUE)
  #ctrlTable <- read.table(ctrl_path, header = TRUE)
  
  exp <- data.matrix(expTable)
  ctrl <- data.matrix(ctrlTable)
  
  genes <- c()
  labels <- c()
  predictionList <- c()
  total_features = 0
  
  for(i in 1:ncol(exp)){
    
    exp_minus_one = data.matrix(exp[,-i])
    exp_test = data.matrix(exp[,i])
    
    for(j in 1:ncol(ctrl)){
      
      ctrl_minus_one = data.matrix(ctrl[,-j])
      ctrl_test = data.matrix(ctrl[,j])
      
      train_matrix = cbind(exp_minus_one, ctrl_minus_one)
      test_matrix = cbind(exp_test, ctrl_test)
      
      #The following block of code is used only in the Mouse Data due to the fact that there was
      #low read coverage. Either genes with any sample having 0 read count in any sample, or genes with
      #more than one 0 per phenotype are dropped due to the sample size of 8.
      
      row_sub = apply(train_matrix, 1, function(row) all(row !=0 ))
      #row_sub = data.matrix(apply(train_matrix, 1, function(row) (sum(row[1:3] ==0 )<=1 & sum(row[4:6] == 0)<=1)))
      train_matrix <- train_matrix[row_sub,]
      test_matrix <- test_matrix[row_sub,]
      
      t_test <- data.matrix(apply(train_matrix,1,function(x){
          obj<-try(t.test(x[1:(ncol(exp)-1)],x[(ncol(exp)):((ncol(exp)-1)+(ncol(ctrl)-1))]), silent=TRUE)
          if (is(obj, "try-error")) return(NA)
          else return(obj$p.value)
      }))
      
      train_matrix <- t(data.matrix(train_matrix[t_test[,1] < p_value & !is.na(t_test[,1]),]))
      test_matrix <- t(data.matrix(test_matrix[t_test[,1] < p_value & !is.na(t_test[,1]),]))
    
      library(randomForest)
      library(pROC)
      library(stringr)
      
      phenotypes <- c(rep(0,ncol(exp)-1), rep(1,ncol(ctrl)-1))
      
      genes <- c(genes, colnames(train_matrix))
      
      total_features = total_features + ncol(test_matrix)
      
      TRAINCV <- randomForest(train_matrix, phenotypes)
      
      prediction <- predict(TRAINCV, test_matrix)
      
      predictionList <- c(predictionList, prediction)
      
      labels <- c(labels, 1, 2)
      
    }
    
  }
  
  print(predictionList)
  
  roc <- roc(labels, predictionList, plot=TRUE)
  print(roc)
  print(ci(roc$auc))
  
  detach("package:pROC", unload=TRUE)
  library(AUC)
  
  labels <- factor(labels)
  
  print("Area Under ROC:")
  print(auc(roc(predictionList, labels)))
  print("Accuracy:")
  print(auc(accuracy(predictionList, labels)))
  print("Sensitivity:")
  print(auc(sensitivity(predictionList, labels)))
  print("Specificity:")
  print(auc(specificity(predictionList, labels)))
  print("Mean Features Selected:")
  print(total_features/(ncol(exp) * ncol(ctrl)))
  
  return(genes)
}
