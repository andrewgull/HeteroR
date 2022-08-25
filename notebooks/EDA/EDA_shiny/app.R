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


features_amp_strain <- read_csv("../data/features_amp_strain.csv")
# rename for consistency
features_amp_strain <- rename(features_amp_strain, "n.plasmids" = n_plasmids)
vars <- names(select(features_amp_strain, -c("strain", "AB", "resistance")))
strains <- features_amp_strain$strain

# Define UI for application that draws a histogram
ui <- fluidPage(

  # Application title
  titlePanel("Explore the data!"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      style = "position:fixed;width:inherit;",
      "Inputs",
      width = 3,
      selectInput("xcol",
                  "X variable:",
                  vars,
                  selected = "n.rep.total"),
      selectInput("ycol", 
                  "Y variable:",
                  vars,
                  selected="ampC.n.rep.tot"),
      sliderInput("size",
                  "Dot size",
                  value=2, 
                  min = 0.5, 
                  max = 10),
      sliderInput("alpha",
                  "Transparency:",
                  value=0.5, 
                  min = 0.1, 
                  max = 1),
      
      selectInput("bar",
                  "Count data",
                  c("n.beta.lac", "n.plasmids", "n.genes.plus.strand", "n.genes.plasmids"),
                  selected="n.beta.lac"),
      
      selectInput("box",
                  "Median & distribution",
                  vars,
                  selected = "med.dist.oriC"),
      radioButtons("trans", 
                    "Y-axis transformation", 
                    choices = c("identity", "log", "sqrt"), 
                    selected = "identity"),
      checkboxInput("notch", 
                    "Notch",
                    value = FALSE),
      
      sliderInput("sample", 
                  "Sample size", 
                  value=10, 
                  min = 10, 
                  max=200)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("dot.plot"),
      plotOutput("bar.plot"),
      plotOutput("box.plot")
    )
  )
)

# Define server 
server <- function(input, output) {

  library(ggplot2)
  library(dplyr)
  
  df <- readr::read_csv("../data/features_amp_strain.csv") %>%
    rename("n.plasmids"=n_plasmids)
  
  # Data for the 1st plot
  selected_data1 <- reactive({
    select(df, input$xcol, input$ycol, resistance)
  })
  
  output$dot.plot <- renderPlot({
    x <- paste0("`",input$xcol,"`")
    y <- paste0("`",input$ycol,"`")
    
    ggplot(selected_data1(), aes_string(x, y)) +
      geom_point(aes(color=resistance), size = input$size, alpha = input$alpha) +
      geom_smooth(method = "lm") +
      scale_color_brewer(palette="Set1", name="Resistance")+
      theme(legend.position = "bottom")+
      ggtitle("Dot plot")
  })
  
  # Data for the 2nd plot
  selected_data2 <- reactive({
    select(df, input$bar, resistance)
  })
  
  output$bar.plot <- renderPlot({
    x <- paste0("`",input$bar,"`")
    ggplot(selected_data2(), aes_string(x))+
      geom_bar(aes(fill=resistance), position="dodge", alpha=0.8)+
      scale_fill_brewer(palette="Set1", name="Resistance")+
      theme(legend.position = "bottom") +
      xlab("N")+
      ylab("count")+
      ggtitle(paste0("Count data: ", input$bar))
  })
  
  # The 3rd plot
  output$box.plot <- renderPlot({
    y <- paste0("`",input$box,"`")
    ggplot(df, aes_string("resistance", y))+
      geom_violin(aes(fill=resistance), alpha=0.7)+
      geom_boxplot(aes(fill=resistance), alpha=0.2, notch = input$notch, varwidth = T)+
      geom_jitter(alpha=0.2, width=0.15, height = 0.1)+
      coord_trans(y = input$trans) +
      scale_fill_brewer(palette = "Set1") +
      xlab("")+ylab("") +
      ggtitle(input$box)+
      guides(fill="none")
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
