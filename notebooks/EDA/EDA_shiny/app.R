# Shiny web application for EDA

library(shiny)
library(tidyverse)

# Read the main data table with features
features_amp_strain <- read_csv("../data/features_amp_strain.csv")

# rename 'n_plasmids' for consistency of column names
features_amp_strain <- rename(features_amp_strain, "n.plasmids" = n_plasmids)

# get vars to use later in UI (selectInput)
vars <- names(select(features_amp_strain, -c("strain", "AB", "resistance")))

# get strain names to use in UI (heat map)
strains <- features_amp_strain$strain

################
### UI part ####
################
ui <- fluidPage(
  
  # bslib::bs_theme_preview(bs_theme(bootswatch = "theme_name"))
  theme = bslib::bs_theme(bootswatch = "darkly"),
  
  # Application title
  titlePanel("EDA"),
  
  # A separate panel for AMP data
  tabsetPanel(
    tabPanel("AMP",
             
             # 1st row with a heat map
             fluidRow(
               # Left column with controls
               column(2, sliderInput("sample", 
                                     "Sample size", 
                                     value = 10, 
                                     min = 10, 
                                     max = length(strains))),
               # Right column with the heat map itself
               column(10, plotOutput("heatmap"))
             ),
             
             # 2nd row with a bar pot for count data
             fluidRow(
               # first column with controls
               column(2, 
                      selectInput(inputId = "bar", label = "Count data",
                                  choices = c("n.beta.lac", "n.plasmids", "n.genes.plus.strand", "n.genes.plasmids"),
                                  selected = "n.beta.lac")
               ),
               # second column with the plot itself
               column(10, plotOutput("bar.plot"))
             ),
             
             # 3rd row with a box plot
             fluidRow(
               # Left column with controls for the box plot
               column(2, 
                      selectInput(inputId = "box", label = "Median & distribution",
                                  choices = vars, selected = "med.dist.oriC"),
                      radioButtons(inputId = "trans", label = "Y-axis transformation", 
                                   choices = c("identity", "log", "sqrt"), selected = "identity", 
                                   choiceNames = c("none", "log", "sqrt")),
                      checkboxInput(inputId = "notch", label = "Notch",
                                    value = FALSE)
               ),
               # Right column with the box plot itself
               column(10, plotOutput("box.plot"))
             ),
             
             # 4th row for a dot plot widget
             fluidRow(
               # Left side with widget's controls
               column(2, 
                      selectInput(inputId = "xcol",label = "X variable",
                                  choices = vars, selected = "n.rep.total"),
                      selectInput(inputId = "ycol", label = "Y variable",
                                  choices = vars, selected = "ampC.n.rep.tot"),
                      sliderInput(inputId = "size", label = "Dot size",
                                  value = 2, min = 0.5, max = 10),
                      sliderInput(inputId = "alpha", label = "Opacity",
                                  value = 0.5, min = 0.1, max = 1)
               ),
               # Right part with the dot plot itself
               column(10, plotOutput("dot.plot"))
             ),
             
    ),
    # Tab for CFX data
    tabPanel("CFX", "no data yet"),
    # Tab for MCN data
    tabPanel("MCN", "no data yet"),
    # Tab for GM data
    tabPanel("GM", "no data yet"),
    # Tab for NFT data
    tabPanel("NFT", "no data yet")
  )
  
  
)

###################
### Server part ###
###################
server <- function(input, output) {
  
  # for plots colored according to the chosen theme above
  thematic::thematic_shiny()
  
  # read data with main features
  df <- readr::read_csv("../data/features_amp_strain.csv") %>%
    rename("n.plasmids" = n_plasmids)
  
  # read data with BL types
  df2 <- readr::read_csv("../data/amp_amr_types_strain.csv")
  
  # Dot plot:
  output$dot.plot <- renderPlot({
    x <- paste0("`",input$xcol,"`")
    y <- paste0("`",input$ycol,"`")
    
    ggplot(df, aes_string(x, y)) +
      geom_point(aes(color = resistance), size = input$size, alpha = input$alpha) +
      geom_smooth(method = "lm") +
      scale_color_brewer(palette = "Set1", name = "Resistance") +
      theme(legend.position = "bottom") +
      ggtitle("Dot plot")
  })
  
  # Bar plot
  output$bar.plot <- renderPlot({
    x <- paste0("`",input$bar,"`")
    ggplot(df, aes_string(x)) +
      geom_bar(aes(fill = resistance), position = "dodge", alpha = 0.8) +
      scale_fill_brewer(palette = "Set1", name = "Resistance") +
      theme(legend.position = "bottom") +
      xlab("N") +
      ylab("count") +
      ggtitle(paste0("Count data: ", input$bar))
  })
  
  # Box plot
  output$box.plot <- renderPlot({
    y <- paste0("`",input$box,"`")
    ggplot(df, aes_string("resistance", y)) +
      geom_violin(aes(fill = resistance), alpha = 0.7) +
      geom_boxplot(aes(fill = resistance), alpha = 0.2, notch = input$notch, varwidth = T) +
      geom_jitter(alpha = 0.2, width = 0.15, height = 0.1) +
      coord_trans(y = input$trans) +
      scale_fill_brewer(palette = "Set1") +
      xlab("") + ylab("") +
      ggtitle(input$box) +
      guides(fill = "none")
  })
  
  # Heat map
  # the data should be made tidy first
  selected_data <- reactive({
    tidyr::gather(slice(df2, 1:input$sample), key = "AMR.type", value = "N", 2:22)
    })
  
  # plot
  output$heatmap <- renderPlot({
    ggplot(selected_data(), aes(strain, AMR.type))+geom_tile(aes(fill = N)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.2, hjust = 1), legend.position = "bottom" ) +
      xlab("") +
      scale_fill_distiller(palette = "Blues", direction = 1) +
      #scale_fill_viridis_c(direction = 1, alpha = 0.8)+
      guides(colour = "colorbar", size = "legend", shape = "legend") +
      ggtitle("Heatmap: beta-lactamase types")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
