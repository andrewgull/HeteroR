#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(ggplot2)
library(data.table)


get_table <- function(table_path, filt=FALSE){
  df <- fread(table_path)
  df <- select(df, -c(end_1, end_2, X, V1))
  if (filt){
    df <- filter(df, spans_center == "yes")
  } 
  return(df)  
}


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Genomic repeats statistics"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          style = "position:fixed;width:inherit;",
          "Inputs",
          width = 3,
            sliderInput("sample_size",
                        "Sample size:",
                        min = 1,
                        max =237,
                        value = 10),
            sliderInput("seed", 
                       "Set seed:",
                       min = 1,
                       max = 99,
                       value = 10
                       ),
            sliderInput("min_len",
                        "Minimum repeat length",
                        min=20,
                        max=1000,
                        value=20),
            checkboxInput("filter", 
                          "Exclude repeat pairs not spanning RG center",
                          value=FALSE)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput(outputId = "boxPlot"),
           plotOutput("barPlot"),
           plotOutput("boxPlot2"),
           plotOutput("boxPlot3"),
           dataTableOutput("sampleTable")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  sample_data <- reactive({
      # sample data
      set.seed(input$seed)
      # read in the main data table
      repeat_df <- get_table("data/repeats_summary_table.csv", filt = input$filter)
      strains <- unique(repeat_df$strain)
      strains_sample <- sample(strains, input$sample_size, replace = FALSE)
      # get sample table for rendering
      filter(bind_rows(lapply(strains_sample, function(x){filter(repeat_df, strain == x)})), length >= input$min_len)
    })
  
    # render outputs
    output$sampleTable <- renderDataTable({
        sample_data()
      }, 
    options = list(pageLength = 10)
    )
    
    output$boxPlot <- renderPlot({
        ggplot(sample_data(), aes(strain, length)) +
            geom_boxplot(outlier.size = 1.0, outlier.alpha = 0.2, size=0.2) +
            geom_violin(alpha=0.2, size=0.2) +
            theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
            xlab("") +
            ylab("length") +
            scale_y_continuous(breaks=c(20, 100, 500, 1000, 2000, 4000, 6000), trans="log") +
            ggtitle("Repeat lengths distribution")
    })
    
    output$barPlot <- renderPlot({
        # prepare counts
        repeat_counts <- sample_data() %>% group_by(strain) %>% summarise("n_repeats"=n())
        # make a plot
        ggplot(repeat_counts, aes(strain, n_repeats)) +
            geom_bar(stat="identity", fill="steelblue") +
            theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
            xlab("") +
            ylab("counts") +
            ggtitle("Repeat count per strain")
    })
    
    output$boxPlot2 <- renderPlot({
        # prepare counts
        repeat_counts <- sample_data() %>% group_by(record_id, strain) %>% summarise("n_repeats"=n())
        # make a plot
        ggplot(repeat_counts, aes(strain, n_repeats)) +
            geom_boxplot(outlier.size = 1.0, outlier.alpha = 0.5, size=0.2) +
            geom_violin(alpha=0.2, size=0.2) +
            theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
            xlab("") +
            ylab("counts") +
            coord_trans(y="sqrt") +
            scale_y_continuous(breaks=c(1, 10, 20, 50, 100, 500, 1000, 2000)) +
            ggtitle("Repeat counts per RG region")
        
    })
    
    output$boxPlot3 <- renderPlot({
        ggplot(sample_data(), aes(strain, AR_length)) +
            geom_boxplot(alpha=0.5, outlier.size = 1.0, outlier.alpha = 0.5, size=0.2) +
            geom_violin(alpha=0.2, size=0.2) +
            theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
            xlab("") +
            ylab("length") +
            ggtitle("Amplifiable region length")
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
