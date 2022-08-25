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
  # Set a theme
  # To preview themes use 
  # bslib::bs_theme_preview(bs_theme(bootswatch = "theme_name"))
  theme = bslib::bs_theme(bootswatch = "darkly"),
  
  # Application title
  titlePanel("Explore the data!"),
  
  fluidRow(
    column(3, 
           selectInput("xcol",
                       "X variable",
                       vars,
                       selected = "n.rep.total"),
           selectInput("ycol", 
                       "Y variable",
                       vars,
                       selected="ampC.n.rep.tot"),
           sliderInput("size",
                       "Dot size",
                       value=2, 
                       min = 0.5, 
                       max = 10),
           sliderInput("alpha",
                       "Opacity",
                       value=0.5, 
                       min = 0.1, 
                       max = 1)
           ),
    column(9, plotOutput("dot.plot"))
  ),
  
  fluidRow(
    column(3, 
            selectInput("bar",
                        "Count data",
                        c("n.beta.lac", "n.plasmids", "n.genes.plus.strand", "n.genes.plasmids"),
                        selected="n.beta.lac")
            ),
    column(9, plotOutput("bar.plot"))
  ),
  
  fluidRow(
    column(3, 
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
                         value = FALSE)
           ),
    column(9, plotOutput("box.plot"))
  ),
  
  fluidRow(
    column(12, sliderInput("sample", 
                           "Sample size", 
                           value=10, 
                           min = 10, 
                           max=length(strains))), 
  ),
  
  fluidRow(
    column(12, plotOutput("heatmap"))
  )
)

# Define server 
server <- function(input, output) {

  thematic::thematic_shiny()
  
  # read data
  df <- readr::read_csv("../data/features_amp_strain.csv") %>%
    rename("n.plasmids"=n_plasmids)
  
  df2 <- readr::read_csv("../data/amp_amr_types_strain.csv")
  
  # Dot plot:

  output$dot.plot <- renderPlot({
    x <- paste0("`",input$xcol,"`")
    y <- paste0("`",input$ycol,"`")
    
    ggplot(df, aes_string(x, y)) +
      geom_point(aes(color=resistance), size = input$size, alpha = input$alpha) +
      geom_smooth(method = "lm") +
      scale_color_brewer(palette="Set1", name="Resistance")+
      theme(legend.position = "bottom")+
      ggtitle("Dot plot")
  })
  
  # Bar plot
  
  output$bar.plot <- renderPlot({
    x <- paste0("`",input$bar,"`")
    ggplot(df, aes_string(x))+
      geom_bar(aes(fill=resistance), position="dodge", alpha=0.8)+
      scale_fill_brewer(palette="Set1", name="Resistance")+
      theme(legend.position = "bottom") +
      xlab("N")+
      ylab("count")+
      ggtitle(paste0("Count data: ", input$bar))
  })
  
  # Box plot
  
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
  
  # Heat map
  # select data and make tidy
  selected_data <- reactive({
    tidyr::gather(slice(df2, 1:input$sample), key="AMR.type", value="N", 2:22)
    })

  output$heatmap <- renderPlot({
    ggplot(selected_data(), aes(strain, AMR.type))+geom_tile(aes(fill=N))+
      theme(axis.text.x = element_text(angle = 45, vjust = 0.2, hjust=1))+
      xlab("")+
      scale_fill_distiller(palette = "Blues", direction = 1) +
      #scale_fill_viridis_c(direction = 1, alpha = 0.8)+
      guides(colour = "colorbar", size = "legend", shape = "legend") +
      ggtitle("Beta-lactamase types")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
