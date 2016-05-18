library(shiny)
ui <- fluidPage(pageWithSidebar(
  headerPanel("Phenotype Prediction Pipeline"), 
  sidebarPanel(
    
  ), 
  mainPanel(
    
  )))
server <- function(input, output) {
  
}
shinyApp(ui, server)