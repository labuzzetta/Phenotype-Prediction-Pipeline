##Phenotype Prediction Pipeline shinyApp.R##

library(shiny)
library(shinydashboard)

ui <- dashboardPage(
  dashboardHeader(title = "MVP Pipeline"),
  dashboardSidebar(
    menuItem("Predict Unknowns", tabName = "unknowns", icon = icon("magic"),
             badgeLabel = "new", badgeColor = "green"),
    menuItem("Cross Validation", tabName = "validation", icon = icon("line-chart"),
             badgeLabel = "new", badgeColor = "green")
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "unknowns",
              h2("Predict Phenotype of Unknown Samples")
      ),
      
      tabItem(tabName = "validation",
              h2("Cross Validate Machine Learning Methods")
      )
    )
  )
)

server <- function(input, output) {
  
}

shinyApp(ui, server)