#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

df <- read.csv("data/heatmap.csv")
df_tidy <- gather(df, "strain", "n", 3:405)  # the last number equals number of strains plus 2
df_tidy$strain <- sub("_rgi_table", "", df_tidy$strain)
families <- unique(df$gene_family)

library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

# Define UI for application
ui <- fluidPage(
    # Application title
    titlePanel("Heatmap"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            style = "position:fixed;width:inherit;",
            "Inputs",
            width = 2,
            checkboxGroupInput("AMR_family", label = "AMR family", choices = families, selected = families)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("heatmap", width="100%")
        )
    )
)

# Define server logic 
server <- function(input, output) {
    
    sample_data <- reactive({
        #filter(df_tidy, gene_family == input$AMR_family)
        bind_rows(lapply(input$AMR_family, function(x){filter(df_tidy, gene_family == x)}))
    })
    
    output$heatmap <- renderPlot({
        ggplot(sample_data(), aes(strain, gene)) + 
            geom_tile(aes(fill=n)) +
            #scale_fill_gradient(low="white", high="steelblue", breaks=c(0, 1, 2)) +
            scale_fill_viridis() +
            theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, size=8)) +
            xlab("")
    }, height=1200, width=2800)
}

# Run the application 
shinyApp(ui = ui, server = server)
