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
            h2("Predict Phenotype of Unknown Samples"),
            box(
              title = "Inputs", status = "warning", solidHeader = TRUE,
              "Select the range of columns in the uploaded expression file corresponding to each phenotype with the sliders.",
              sliderInput("colPheno1", "Phenotype 1:", min = 1, max = 30, value = c(1,15), ticks = FALSE),
              sliderInput("colPheno2", "Phenotype 2:", min = 1, max = 30, value = c(16,30), ticks = FALSE),
              "Input a p-value to identify differentially expressed transcripts.",
              textInput("text", label = NULL, value = 0.05)
            ),
            box(
              title = "ROC Curve Analysis", status = "primary", solidHeader = TRUE,
              collapsible = TRUE,
              plotOutput("plot3", height = 250)
            )
    ),
    
    tabItem(tabName = "validation",
            h2("Cross Validate Machine Learning Methods")
    )
  )
)


ui <- dashboardPage(
  dashboardHeader(title = "MVP Pipeline"), sidebar, body
)

server <- function(input, output) {
  
}

shinyApp(ui, server)