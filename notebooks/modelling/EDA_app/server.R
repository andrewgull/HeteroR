# server

server <- function(input, output) {
  
  # for plots colored according to the chosen theme above
  thematic::thematic_shiny()
  
  # Histogram original
  output$histogram.orig <- renderPlot({
    features_strain %>% 
      ggplot(aes_string(x = input$var)) + 
      geom_histogram(bins = 30, color = "white", fill = "blue", alpha = 1/3) + 
      ggtitle("Original data")
  })
  
  # Histogram norm
  output$histogram.norm <- renderPlot({
    features_norm %>% 
      ggplot(aes_string(x = input$var)) + 
      geom_histogram(bins = 30, color = "white", fill = "red", alpha = 1/3) + 
      ggtitle("Normalized data")
  })
  
  # Dot plot:
  output$dot.plot <- renderPlot({
    x <- paste0("`", input$xcol,"`")
    y <- paste0("`", input$ycol,"`")
    
    ggplot(features_strain, aes_string(x, y)) +
      geom_point(aes(color = resistance), size = 2, alpha = 0.5) +
      scale_color_brewer(palette = "Set1", name = "Resistance") +
      coord_trans(y = input$y.trans, x = input$x.trans) +
      theme(legend.position = "bottom") +
      ggtitle("Dot plot")
  })
  
  # Dot plot hover table
  output$dot.plot.data <- renderTable({
    req(input$plot_brush)
    brushedPoints(select(features_strain, strain, input$xcol, input$ycol), input$plot_brush)
  })
  
  # Bar plot A
  output$bar.plot <- renderPlot({
    x <- paste0("`",input$bar,"`")
    ggplot(features_strain, aes_string(x)) +
      geom_bar(aes(fill = resistance), position = input$bar.type, alpha = 0.8) +
      scale_fill_brewer(palette = "Set1", name = "Resistance") +
      theme(legend.position = "bottom") +
      xlab("N") +
      ylab("count / proportion") +
      ggtitle(paste0("Count data: ", input$bar))
  })
  
  # Box plot
  output$box.plot <- renderPlot({
    y <- paste0("`",input$box,"`")
    ggplot(features_strain, aes_string("resistance", y)) +
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
  output$heatmap <- renderPlotly({
    vals <- scales::rescale(c(0:6))
    o <- order(vals, decreasing = FALSE)
    cols <- scales::col_numeric("Blues", domain = NULL)(vals)
    colz <- setNames(data.frame(vals[o], cols[o]), NULL)
    plot_ly(z = bl_amr_types_strain_td$N, 
                   x = bl_amr_types_strain_td$strain, 
                   y = bl_amr_types_strain_td$AMR.type, 
                   type = "heatmap", 
                   colorscale = colz)
  })
  
  output$umap <- renderPlot({
    umap_rec <- data_rec %>%
      step_umap(all_numeric_predictors(), 
                num_comp = input$umap.comp, 
                neighbors = input$umap.neighb,
                min_dist = input$umap.dist) 
    
    data_umap <- prep(umap_rec, retain = TRUE)
    
    data_umap$template %>% 
      plot_validation_results() +
      ggtitle("Supervised Uniform Manifold Approximation and Projection")
  }, height = 1000)
}