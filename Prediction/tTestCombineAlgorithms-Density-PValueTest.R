t.test.select.AllAlgorithms <- function(expTable, ctrlTable, p_value, quant) {
  
  #expTable <- read.table(exp_path, header = TRUE)
  #ctrlTable <- read.table(ctrl_path, header = TRUE)
  
  exp <- data.matrix(expTable)
  ctrl <- data.matrix(ctrlTable)
  
  genes <- c()
  labels <- c()
  predictionListElasticNet <- c()
  predictionListRandomForest <- c()
  predictionListSPLS <- c()
  total_features = 0
  
  library(extremevalues)
  
  for(i in 1:ncol(exp)){
    
    exp_minus_one = data.matrix(exp[,-i])
    exp_test = data.matrix(exp[,i])
    
    for(j in 1:ncol(ctrl)){
      
      ctrl_minus_one = data.matrix(ctrl[,-j])
      ctrl_test = data.matrix(ctrl[,j])
      
      train_matrix = cbind(exp_minus_one, ctrl_minus_one)
      test_matrix = cbind(exp_test, ctrl_test)

      complete_test <- cbind(train_matrix, test_matrix)
      train_matrix <- train_matrix[complete.cases(complete_test),]
      test_matrix <- test_matrix[complete.cases(complete_test),]
      

      row_sub = apply(train_matrix, 1, function(row) (all(row != 0)))
      #row_sub = data.matrix(apply(train_matrix, 1, function(row){(sum(row[1:(ncol(exp)-1)] ==0 )<=((ncol(exp)-1)/5) & sum(row[(ncol(exp)):((ncol(exp)-1)+(ncol(ctrl)-1))] == 0)<=((ncol(ctrl)-1)/5))}))
      train_matrix <- train_matrix[row_sub,]
      test_matrix <- test_matrix[row_sub,]
      
      row_sub = data.matrix(apply(test_matrix, 1, function(row) (all(row != 0))))
      train_matrix <- train_matrix[row_sub,]
      test_matrix <- test_matrix[row_sub,]
      
      t_test <- data.matrix(apply(train_matrix,1,function(x){
        obj<-try(t.test(x[1:(ncol(exp)-1)],x[(ncol(exp)):((ncol(exp)-1)+(ncol(ctrl)-1))]), silent=TRUE)
        if (is(obj, "try-error")) return(NA)
        else return(obj$p.value)
      }))
      
      train_matrix <- data.matrix(train_matrix[t_test[,1] < p_value & !is.na(t_test[,1]),])
      test_matrix <- data.matrix(test_matrix[t_test[,1] < p_value & !is.na(t_test[,1]),])
      row.names(test_matrix) <- row.names(train_matrix)
      
      input <- cbind(train_matrix, test_matrix)
      
      returned <- apply(input,1,try(function(row){
        curve <- density(row[1:(ncol(exp)-1)])
        test1 <- (ncol(exp)-1)+(ncol(ctrl)-1)+1
        test2 <- (ncol(exp)-1)+(ncol(ctrl)-1)+2
        data <- c(mean(curve$x), var(curve$x),dnorm(row[test1], mean(curve$x), sd(curve$x)),dnorm(row[test2], mean(curve$x), sd(curve$x)))
      }))
      
      returned2 <- apply(input,1,try(function(row){
        curve <- density(row[(ncol(exp)):((ncol(exp)-1)+(ncol(ctrl)-1))])
        test1 <- (ncol(exp)-1)+(ncol(ctrl)-1)+1
        test2 <- (ncol(exp)-1)+(ncol(ctrl)-1)+2
        data <- c(mean(curve$x), var(curve$x),dnorm(row[test1], mean(curve$x), sd(curve$x)),dnorm(row[test2], mean(curve$x), sd(curve$x)))
      }))
      
      #if(is.vector(returned)){returned <- do.call(cbind, returned)}
      #if(is.vector(returned2)){returned2 <- do.call(cbind, returned2)}
      
      diffp <- apply(t(rbind(returned[3:4,],returned2[3:4,])),1,try(function(col){
        data <- c(abs(col[1] - col[3]) * abs(col[2] - col[4]))
      }))
    
      diffp <- data.matrix(diffp)
      

      
        #returned <- returned[,intersect(colnames(returned), colnames(returned2))]
        #returned2 <- returned2[,intersect(colnames(returned), colnames(returned2))]
      
        #complete_test <- t(rbind(returned, returned2))
        #returned <- returned[,complete.cases(complete_test)]
        #returned2 <- returned2[,complete.cases(complete_test)]
  
      #print(sort(abs(returned[1,] - returned2[1,])/(returned[2,] + returned2[2,]),decreasing = T)[1:20])
      
      #plot(density(train_matrix[,1:10]));lines(density(train_matrix[,11:19]))
      
      #test_matrix <- t(test_matrix[sort(names(sort(abs(returned[1,] - returned2[1,])/(returned[2,] + returned2[2,]),decreasing = T)[1:nSelected])),])
      #train_matrix <- t(train_matrix[sort(names(sort(abs(returned[1,] - returned2[1,])/(returned[2,] + returned2[2,]),decreasing = T)[1:nSelected])),])
      
      #print(quantile(abs(returned[1,] - returned2[1,])/(returned[2,] + returned2[2,]),0.9))
      
      #Squared
      #test_matrix <- t(test_matrix[(returned[1,] - returned2[1,])^2/(returned[2,] + returned2[2,])>quantile((returned[1,] - returned2[1,])^2/(returned[2,] + returned2[2,]),quant),])
      #train_matrix <- t(train_matrix[(returned[1,] - returned2[1,])^2/(returned[2,] + returned2[2,])>quantile((returned[1,] - returned2[1,])^2/(returned[2,] + returned2[2,]),quant),])
      
      
      #Bhattachayya
      #print(quantile((1/4)*log((1/4)*((returned[2,]/returned2[2,]) + (returned[2,]/returned2[2,]) + 2)) + (1/4)*((returned[1,] - returned2[1,])^2/(returned[2,] + returned2[2,])),.9))
      #test_matrix <- t(test_matrix[(1/4)*log((1/4)*((returned[2,]/returned2[2,]) + (returned[2,]/returned2[2,]) + 2)) + (1/4)*((returned[1,] - returned2[1,])^2/(returned[2,] + returned2[2,]))>quantile((1/4)*log((1/4)*((returned[2,]/returned2[2,]) + (returned[2,]/returned2[2,]) + 2)) + (1/4)*((returned[1,] - returned2[1,])^2/(returned[2,] + returned2[2,])),quant),])
      #train_matrix <- t(train_matrix[(1/4)*log((1/4)*((returned[2,]/returned2[2,]) + (returned[2,]/returned2[2,]) + 2)) + (1/4)*((returned[1,] - returned2[1,])^2/(returned[2,] + returned2[2,]))>quantile((1/4)*log((1/4)*((returned[2,]/returned2[2,]) + (returned[2,]/returned2[2,]) + 2)) + (1/4)*((returned[1,] - returned2[1,])^2/(returned[2,] + returned2[2,])),quant),])
      
      #Hellinger
      #print(quantile(1 - sqrt(  (2*sqrt(returned[2,])*sqrt(returned2[2,])) / (returned[2,] + returned2[2,]) )*exp(-(1/4)*((returned[1,] - returned2[1,])^2/(returned[2,] + returned2[2,]))),0.9))
      #test_matrix <- t(test_matrix[(1 - sqrt(  (2*sqrt(returned[2,])*sqrt(returned2[2,])) / (returned[2,] + returned2[2,]) )*exp(-(1/4)*((returned[1,] - returned2[1,])^2/(returned[2,] + returned2[2,]))))>quantile((1/4)*log((1/4)*((returned[2,]/returned2[2,]) + (returned[2,]/returned2[2,]) + 2)) + (1/4)*((returned[1,] - returned2[1,])^2/(returned[2,] + returned2[2,])),quant),])
      #train_matrix <- t(train_matrix[(1 - sqrt(  (2*sqrt(returned[2,])*sqrt(returned2[2,])) / (returned[2,] + returned2[2,]) )*exp(-(1/4)*((returned[1,] - returned2[1,])^2/(returned[2,] + returned2[2,]))))>quantile((1/4)*log((1/4)*((returned[2,]/returned2[2,]) + (returned[2,]/returned2[2,]) + 2)) + (1/4)*((returned[1,] - returned2[1,])^2/(returned[2,] + returned2[2,])),quant),])
      
      test_matrix <- t(test_matrix[abs(returned[1,] - returned2[1,])/(returned[2,] + returned2[2,])>quantile(abs(returned[1,] - returned2[1,])/(returned[2,] + returned2[2,]),quant) & diffp[,1] > quantile(diffp[,1], 0),])
      train_matrix <- t(train_matrix[abs(returned[1,] - returned2[1,])/(returned[2,] + returned2[2,])>quantile(abs(returned[1,] - returned2[1,])/(returned[2,] + returned2[2,]),quant) & diffp[,1] > quantile(diffp[,1], 0),])
      
      
#     library(caret)
#     descrCorr <- cor(train_matrix)
#     highCorr <- findCorrelation(descrCorr, 0.9)
#     train_matrix <- train_matrix[,-highCorr]
#     test_matrix <- test_matrix[,-highCorr]
      
     #train_matrix <- train_matrix[,-nearZeroVar(train_matrix)]
     #test_matrix <- test_matrix[,-nearZeroVar(train_matrix)]
              
      heatmap(cbind(t(train_matrix), t(test_matrix)))
      
      library(randomForest)
      library(pROC)
      library(stringr)
      
      phenotypes <- c(rep(0,ncol(exp)-1), rep(1,ncol(ctrl)-1))
      
      total_features = total_features + ncol(test_matrix)
      
      RandomForestCV <- randomForest(train_matrix, phenotypes)
      
      predictionRandomForest <- predict(RandomForestCV, test_matrix)
      
      predictionListRandomForest <- c(predictionListRandomForest, predictionRandomForest)
      
      print(predictionRandomForest)
      
      labels <- c(labels, 0, 1)
      
      library(glmnet)
      library(pROC)
      library(stringr)

      ElasticNetCV <- cv.glmnet(train_matrix, phenotypes, nfolds=nrow(train_matrix), type.measure="deviance")
      
      predictionElasticNet <- predict(ElasticNetCV, test_matrix)
      
      print(predictionElasticNet)
      
      predictionListElasticNet <- c(predictionListElasticNet, predictionElasticNet)
      
      library(spls)
      library(pROC)
      library(stringr)
      
      SPLSCV <- spls(train_matrix, phenotypes, eta = 0.4, K = 3)
      
      predictionSPLS <- predict.spls(SPLSCV, test_matrix)
      
      print(predictionSPLS)
      
      predictionListSPLS <- c(predictionListSPLS, predictionSPLS)
      
    }
    
  }
  
  print(predictionListRandomForest)
  
  roc <- roc(labels, predictionListRandomForest, plot=TRUE)
  print(roc)
  print(ci(roc$auc))
  
  detach("package:pROC", unload=TRUE)
  library(AUC)
  
  labels <- factor(labels)
  
  print("RandomForest")
  print("Area Under ROC:")
  try(print(auc(roc(predictionListRandomForest, labels))))
  print("Accuracy:")
  try(print(auc(accuracy(predictionListRandomForest, labels))))
  print("Sensitivity:")
  try(print(auc(sensitivity(predictionListRandomForest, labels))))
  print("Specificity:")
  try(print(auc(specificity(predictionListRandomForest, labels))))
  print("Mean Features Selected:")
  try(print(total_features/(ncol(exp) * ncol(ctrl))))
  
  library(pROC)
  
  print(predictionListElasticNet)
  
  roc <- roc(labels, predictionListElasticNet, plot=TRUE)
  print(roc)
  print(ci(roc$auc))
  
  detach("package:pROC", unload=TRUE)
  library(AUC)
  
  labels <- factor(labels)
  
  print("ElasticNet")
  print("Area Under ROC:")
  try(print(auc(roc(predictionListElasticNet, labels))))
  print("Accuracy:")
  try(print(auc(accuracy(predictionListElasticNet, labels))))
  print("Sensitivity:")
  try(print(auc(sensitivity(predictionListElasticNet, labels))))
  print("Specificity:")
  try(print(auc(specificity(predictionListElasticNet, labels))))
  print("Mean Features Selected:")
  try(print(total_features/(ncol(exp) * ncol(ctrl))))
  
  library(pROC)
  
  print(predictionListSPLS)
  
  roc <- roc(labels, predictionListSPLS, plot=TRUE)
  print(roc)
  print(ci(roc$auc))
  
  detach("package:pROC", unload=TRUE)
  library(AUC)
  
  labels <- factor(labels)
  
  print("SPLS")
  print("Area Under ROC:")
  try(print(auc(roc(predictionListSPLS, labels))))
  print("Accuracy:")
  try(print(auc(accuracy(predictionListSPLS, labels))))
  print("Sensitivity:")
  try(print(auc(sensitivity(predictionListSPLS, labels))))
  print("Specificity:")
  try(print(auc(specificity(predictionListSPLS, labels))))
  print("Mean Features Selected:")
  try(print(total_features/(ncol(exp) * ncol(ctrl))))
  
}
