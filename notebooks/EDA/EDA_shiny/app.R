#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)


features_amp_strain <- read_csv("data/features_amp_strain.csv")
# rename for consistency
features_amp_strain <- rename(features_amp_strain, "n.plasmids" = n_plasmids)
vars <- names(select(features_amp_strain, -c("strain", "AB", "resistance")))


# Define UI for application that draws a histogram
ui <- fluidPage(

  # Application title
  titlePanel("Pairwise dot plots"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      style = "position:fixed;width:inherit;",
      "Inputs",
      width = 3,
      selectInput("xcol",
                  "X variable:",
                  vars,
                  selected = "n.plasmids"),
      selectInput("ycol", 
                  "Y variable:",
                  vars,
                  selected="n.beta.lac"),
      sliderInput("size",
                  "Dot size",
                  value=2, 
                  min = 0.5, 
                  max = 10),
      sliderInput("alpha",
                  "Transparency:",
                  value=0.5, 
                  min = 0.1, 
                  max = 1)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("dot.plot"),
      plotOutput("bar.plot")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  library(ggplot2)
  library(dplyr)
  
  df <- readr::read_csv("data/features_amp_strain.csv") %>%
    rename("n.plasmids"=n_plasmids)
  
  # Data for the 1st plot
  selectedData <- reactive({
    select(df, input$xcol, input$ycol, resistance)
  })
  
  output$dot.plot <- renderPlot({
    x <- paste0("`",input$xcol,"`")
    y <- paste0("`",input$ycol,"`")
    
    ggplot(selectedData(), aes_string(x, y)) +
      geom_point(aes(color=resistance), size = input$size, alpha = input$alpha) +
      scale_color_brewer(palette="Set1", name="Resistance")+
      theme(legend.position = "bottom")
  })
  
  # Data for the 2nd plot
  selectedData <- reactive({
    select(df, input$xcol, input$ycol, resistance)
  })
  
  output$bar.plot <- renderPlot({
    x <- paste0("`",input$plot,"`")
    ggplot(selected_data(), aes_string(x))+
      geom_bar(aes(fill=resistance), position="dodge", alpha=0.8)+
      scale_fill_brewer(palette="Set1", name="Resistance")+
      theme(legend.position = "bottom") +
      xlab("BL")+
      ylab("count")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
