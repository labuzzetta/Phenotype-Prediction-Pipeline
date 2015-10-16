t.test.select.randomForest.combine <- function(expTable_Gene, ctrlTable_Gene, expTable_Iso, ctrlTable_Iso, nonZero, p_value) {
  
  #expTable <- read.table(exp_path, header = TRUE)
  #ctrlTable <- read.table(ctrl_path, header = TRUE)
  
  expG <- data.matrix(expTable_Gene)
  ctrlG <- data.matrix(ctrlTable_Gene)
  
  expI <- data.matrix(expTable_Iso)
  ctrlI <- data.matrix(ctrlTable_Iso)
  
  genes <- c()
  labels <- c()
  predictionList <- c()
  total_features = 0
  
  for(i in 1:ncol(expG)){
    
    exp_minus_one_G = data.matrix(expG[,-i])
    exp_test_G = data.matrix(expG[,i])
    
    exp_minus_one_I = data.matrix(expI[,-i])
    exp_test_I = data.matrix(expI[,i])
    
    for(j in 1:ncol(ctrlG)){
      
      ctrl_minus_one_G = data.matrix(ctrlG[,-j])
      ctrl_test_G = data.matrix(ctrlG[,j])
      
      ctrl_minus_one_I = data.matrix(ctrlI[,-j])
      ctrl_test_I = data.matrix(ctrlI[,j])
      
      train_matrix_G = cbind(exp_minus_one_G, ctrl_minus_one_G)
      test_matrix_G = cbind(exp_test_G, ctrl_test_G)
      
      train_matrix_I = cbind(exp_minus_one_I, ctrl_minus_one_I)
      test_matrix_I = cbind(exp_test_I, ctrl_test_I)
      
      #The following two blocks of code are used only in the Mouse Data due to the fact that there was
      #low read coverage. Either genes with any sample having 0 read count in any sample, or genes with
      #more than one 0 per phenotype are dropped due to the sample size of 8.
      
      row_sub = apply(train_matrix_G, 1, function(row) all(row !=0 ))
      #row_sub = data.matrix(apply(train_matrix, 1, function(row) (sum(row[1:3] ==0 )<=1 & sum(row[4:6] == 0)<=1)))
      train_matrix_G <- train_matrix_G[row_sub,]
      test_matrix_G <- test_matrix_G[row_sub,]
      
      row_sub = apply(train_matrix_I, 1, function(row) all(row !=0 ))
      #row_sub = data.matrix(apply(train_matrix, 1, function(row) (sum(row[1:3] ==0 )<=1 & sum(row[4:6] == 0)<=1)))
      train_matrix_I <- train_matrix_I[row_sub,]
      test_matrix_I <- test_matrix_I[row_sub,]
      
      t_test_G <- data.matrix(apply(train_matrix_G,1,function(x){
        obj<-try(t.test(x[1:(ncol(expG)-1)],x[(ncol(expG)):((ncol(expG)-1)+(ncol(ctrlG)-1))]), silent=TRUE)
        if (is(obj, "try-error")) return(NA)
        else return(obj$p.value)
      }))
      
      t_test_I <- data.matrix(apply(train_matrix_I,1,function(x){
        obj<-try(t.test(x[1:(ncol(expI)-1)],x[(ncol(expI)):((ncol(expI)-1)+(ncol(ctrlI)-1))]), silent=TRUE)
        if (is(obj, "try-error")) return(NA)
        else return(obj$p.value)
      }))
      
      train_matrix_I <- t(data.matrix(train_matrix_I[t_test_I[,1] < p_value & !is.na(t_test_I[,1]),]))
      test_matrix_I <- t(data.matrix(test_matrix_I[t_test_I[,1] < p_value & !is.na(t_test_I[,1]),]))
         
      train_matrix_G <- t(data.matrix(train_matrix_G[t_test_G[,1] < p_value & !is.na(t_test_G[,1]),]))
      test_matrix_G <- t(data.matrix(test_matrix_G[t_test_G[,1] < p_value & !is.na(t_test_G[,1]),]))

      train_matrix <- cbind(train_matrix_G, train_matrix_I)
      test_matrix <- cbind(test_matrix_G, test_matrix_I)
      
      
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
