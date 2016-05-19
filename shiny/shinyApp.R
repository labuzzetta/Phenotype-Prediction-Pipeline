##Phenotype Prediction Pipeline shinyApp.R##

library(shiny)
library(shinydashboard)

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Predict Unknowns", tabName = "unknowns", icon = icon("magic"),
           badgeLabel = "new", badgeColor = "green"),
    menuItem("Cross Validation", tabName = "validation", icon = icon("line-chart"),
           badgeLabel = "new", badgeColor = "green"),
    menuItem("Clone via Github", icon = icon("github"), 
           href = "https://github.com/clabuzze/Phenotype-Prediction-Pipeline.git"),
    fileInput(inputId = "data", label="Upload Expression Table")
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "unknowns",
            h2("Predict Phenotype of Unknown Samples")

    ),
    
    tabItem(tabName = "validation",
            h2("Cross Validate Machine Learning Methods"),
            box(
              title = "Inputs", status = "warning", solidHeader = TRUE,
              "Select the range of columns in the uploaded expression file corresponding to each phenotype with the sliders.",
              sliderInput("colPheno1", "Phenotype 1:", min = 1, max = 30, value = c(1,15), ticks = FALSE),
              sliderInput("colPheno2", "Phenotype 2:", min = 1, max = 30, value = c(16,30), ticks = FALSE),
              "Input a p-value to identify differentially expressed transcripts.",
              textInput("pValue", label = NULL, value = 0.05),
              radioButtons("SelFil",label = "Select Filtering Method",choices = list("None", "MVP"),inline = TRUE),
              actionButton("run.validate", "Click to run validation")
            ),
            box(
              title = "Random Forest ROC Curve Analysis", status = "primary", solidHeader = TRUE,
              collapsible = TRUE,
              plotOutput("rocRF.plot", height = 180),
              textOutput("rocRF.mla"),
              textOutput("rocRF.roc")
            ),
            box(
              title = "Elastic Net ROC Curve Analysis", status = "primary", solidHeader = TRUE,
              collapsible = TRUE,
              plotOutput("rocEN.plot", height = 180),
              textOutput("rocEN.mla"),
              textOutput("rocEN.roc")
            )
    )
  )
)


ui <- dashboardPage(
  dashboardHeader(title = "MVP Pipeline"), sidebar, body
)

server <- function(input, output) {
  
  validate <- observeEvent(input$run.validate, {
    
    tableIn <- read.table(input$data$datapath, header=T)
    
    expTable <- data.matrix(tableIn[,input$colPheno1[1]:input$colPheno1[2]])
    ctrlTable <- data.matrix(tableIn[,input$colPheno2[1]:input$colPheno2[2]])
    
    p_value <- input$pValue

    exp <- data.matrix(expTable)
    ctrl <- data.matrix(ctrlTable)
    
    genes <- c()
    labels <- c()
    predictionListElasticNet <- c()
    predictionListRandomForest <- c()
    predictionListSPLS <- c()
    total_features = 0
    quant = 0
    
    MVPq <- FALSE
    if(input$SelFil == "MVP"){
      MVPq <- TRUE
    }
    
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
        
        if(MVPq == TRUE){
          row_sub = apply(train_matrix, 1, function(row) (all(row != 0)))
          train_matrix <- train_matrix[row_sub,]
          test_matrix <- test_matrix[row_sub,]
          
          row_sub = data.matrix(apply(test_matrix, 1, function(row) (all(row != 0))))
          train_matrix <- train_matrix[row_sub,]
          test_matrix <- test_matrix[row_sub,]
          
          quant = 0.9
        }
        
        train_matrix <- train_matrix * 100
        test_matrix <- test_matrix * 100
        train_matrix <- round(train_matrix,10)
        test_matrix <- round(test_matrix,10)
        
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
          data <- (c(mean(curve$x), var(curve$x)))
        }))
        
        returned2 <- apply(input,1,try(function(row){
          curve <- density(row[(ncol(exp)):((ncol(exp)-1)+(ncol(ctrl)-1))])
          test1 <- (ncol(exp)-1)+(ncol(ctrl)-1)+1
          test2 <- (ncol(exp)-1)+(ncol(ctrl)-1)+2
          data <- (c(mean(curve$x), var(curve$x)))
        }))
        
        test_matrix <- t(test_matrix[abs(returned[1,] - returned2[1,])/(returned[2,] + returned2[2,])>quantile(abs(returned[1,] - returned2[1,])/(returned[2,] + returned2[2,]),quant),])
        train_matrix <- t(train_matrix[abs(returned[1,] - returned2[1,])/(returned[2,] + returned2[2,])>quantile(abs(returned[1,] - returned2[1,])/(returned[2,] + returned2[2,]),quant),]) 

        test_matrix <- test_matrix[,apply(train_matrix,2,var)>0.1e-50]
        train_matrix <- train_matrix[,apply(train_matrix,2,var)>0.1e-50]
        
        #heatmap(cbind(t(train_matrix), t(test_matrix)))
        
        library(randomForest)
        library(pROC)
        library(stringr)
        
        phenotypes <- c(rep(0,ncol(exp)-1), rep(1,ncol(ctrl)-1))
        
        total_features = total_features + ncol(test_matrix)
        
        RandomForestCV <- randomForest(train_matrix, phenotypes)
        
        predictionRandomForest <- predict(RandomForestCV, test_matrix)
        
        predictionListRandomForest <- c(predictionListRandomForest, predictionRandomForest)
        
        labels <- c(labels, 0, 1)
        
        library(glmnet)
        library(pROC)
        library(stringr)
        
        ElasticNetCV <- cv.glmnet(train_matrix, phenotypes, nfolds=nrow(train_matrix), type.measure="deviance")
        
        predictionElasticNet <- predict(ElasticNetCV, test_matrix)
    
        predictionListElasticNet <- c(predictionListElasticNet, predictionElasticNet)
        
      }
      
    }
    
      rocRF <- roc(labels, predictionListRandomForest, plot=FALSE)
      output$rocRF.plot <- renderPlot({plot.roc(rocRF)})
      output$rocRF.mla <- renderPrint("Random Forest")
      output$rocRF.roc <- renderPrint(rocRF$auc)
      
      rocEN <- roc(labels, predictionListElasticNet, plot=FALSE)
      output$rocEN.plot <- renderPlot({plot.roc(rocEN)})
      output$rocEN.mla <- renderPrint("Elastic Net")
      output$rocEN.roc <- renderPrint(rocEN$auc)
    
  })
  
}

shinyApp(ui, server)