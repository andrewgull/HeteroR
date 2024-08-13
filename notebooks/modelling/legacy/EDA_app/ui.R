# ui

ui <- fluidPage(
  
  # bslib::bs_theme_preview(bs_theme(bootswatch = "theme_name"))
  theme = bslib::bs_theme(bootswatch = "yeti"),
  
  # Application title
  titlePanel("PIP-TAZO Hetero-Resistance EDA"),
  
  # A separate panel for basic plots
  tabsetPanel(
    tabPanel("Basic plots",
             # row for a dot plot widget
             fluidRow(
               # Left side with widget's controls
               column(2, 
                      selectInput(inputId = "xcol",label = "X variable",
                                  choices = num_vars, selected = "coverage"),
                      selectInput(inputId = "ycol", label = "Y variable",
                                  choices = num_vars, selected = "read.max.len"),
                      radioButtons(inputId = "x.trans", label = "X-axis transformation", 
                                   choices = c("identity", "log", "sqrt"), selected = "identity", 
                                   choiceNames = c("none", "log", "sqrt")),
                      radioButtons(inputId = "y.trans", label = "Y-axis transformation", 
                                   choices = c("identity", "log", "sqrt"), selected = "identity", 
                                   choiceNames = c("none", "log", "sqrt"))
               ),
               # Right part with the dot plot itself
               column(10, plotOutput("dot.plot", brush = "plot_brush"))
             ),
             
             # 4th row with hover/brush output
             fluidRow(
               column(6),
               column(6, tableOutput("dot.plot.data"))
             ),
             
             # row with a bar pot for count data
             fluidRow(
               # first column with controls
               column(2, 
                      selectInput(inputId = "bar", label = "Count data",
                                  choices = char_vars,
                                  selected = "n.beta."),
                      selectInput(inputId = "bar.type", label = "bar type",
                                  choices = c("stack", "fill", "dodge"),
                                  selected = "dodge")
               ),
               # second column with the plot itself
               column(10, plotOutput("bar.plot"))
             ),
             
             # 2nd row with a box plot
             fluidRow(
               # Left column with controls for the box plot
               column(2, 
                      selectInput(inputId = "box", label = "Median & distribution",
                                  choices = num_vars, selected = "n.plasmids"),
                      radioButtons(inputId = "trans", label = "Y-axis transformation", 
                                   choices = c("identity", "log", "sqrt"), selected = "identity", 
                                   choiceNames = c("none", "log", "sqrt")),
                      checkboxInput(inputId = "notch", label = "Notch",
                                    value = FALSE)
               ),
               # Right column with the box plot itself
               column(10, plotOutput("box.plot"))
             )
             
    ),
    tabPanel("Normalization", "ORQ normalization\n", 
             fluidRow(
               column(width = 2, 
                      selectInput(inputId = "var", 
                                  label = "Variable",
                                  choices = vars_norm,
                                  selected = "n.plasmids")
               ),
               column(width = 10, 
                      plotOutput("histogram.orig"))
             ),
             fluidRow(
               column(width = 2),
               column(width = 10,
                      plotOutput("histogram.norm"))
             )
             
    ),
    # Tab for heatmap of genes
    tabPanel("Heatmap", 
             fluidRow(
               column(width = 12, plotlyOutput("heatmap"))
             )
    ),
    # Tab for PCA plots
    tabPanel("PCA",
             fluidRow(
               column(width = 2, 
                      sliderInput(inputId = "pca.comp",
                                  label = "Number of components",
                                  min = 2,
                                  max = 40,
                                  step = 1,
                                  value = 8)),
               column(width = 10,
                      plotOutput("pca")))
             ),
    # Tab for PCA 3D
    tabPanel("PCA 3D",
             fluidRow(
               column(width = 12,
                      plotlyOutput("pca.3d", height = 1000))
             )),
    # Tab for UMAP
    tabPanel("UMAP", 
             fluidRow(
               column(width = 2, 
                      sliderInput(inputId = "umap.neighb",
                                  label = "Number of neighbors",
                                  min = 4,
                                  max = 40,
                                  step = 1,
                                  value = 15),
                      sliderInput(inputId = "umap.comp",
                                  label = "Number of components",
                                  min = 2,
                                  max = 40,
                                  step = 1,
                                  value = 8),
                      sliderInput(inputId = "umap.dist",
                                  label = "Minimal distance",
                                  min = 0.01,
                                  max = 0.5,
                                  step = 0.01,
                                  value = 0.1)),
               column(width = 10,
                      plotOutput("umap"))
             )),
    tabPanel("UMAP 3D",
             fluidRow(
               column(width = 2, 
                      sliderInput(inputId = "umap.3d.neighb",
                                  label = "Number of neighbors",
                                  min = 4,
                                  max = 40,
                                  step = 1,
                                  value = 15),
                      sliderInput(inputId = "umap.3d.dist",
                                  label = "Minimal distance",
                                  min = 0.01,
                                  max = 0.5,
                                  step = 0.01,
                                  value = 0.1)),
               column(width = 10,
                      plotlyOutput("umap.3d", height = 1000)))
    ),
    tabPanel("UMAP 3D post EDA",
             fluidRow(
               column(width = 12,
                      plotlyOutput("umap.3d.post", height = 1000))
             ))
  )
)