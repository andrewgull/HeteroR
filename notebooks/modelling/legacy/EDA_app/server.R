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
  
  # UMAP plot multi-panel
  output$umap <- renderPlot({
    umap_rec <- data_rec %>%
      step_umap(all_numeric_predictors(), 
                num_comp = input$umap.comp, 
                neighbors = input$umap.neighb,
                min_dist = input$umap.dist,
                outcome = "resistance") 
    
    data_umap <- prep(umap_rec, retain = TRUE)
    
    data_umap$template %>% 
      plot_validation_results() +
      ggtitle("Supervised Uniform Manifold Approximation and Projection")
  }, height = 1000)
  
  # UMAP 3D plotly
  output$umap.3d <- renderPlotly({
    umap_rec <- data_rec %>%
      step_umap(all_numeric_predictors(), 
                num_comp = 3, 
                neighbors = input$umap.3d.neighb,
                min_dist = input$umap.3d.dist,
                outcome = "resistance") 
    
    data_umap <- prep(umap_rec, retain = TRUE)
    
    umap3d <-
      plot_ly(
        data_umap$template,
        x = ~ UMAP1,
        y = ~ UMAP2,
        z = ~ UMAP3,
        color = ~ data_umap$template$resistance,
        colors = c('#cf280c', '#1b56f7')
      )
    
    umap3d <- umap3d %>% 
      add_markers(size = 2, text = ~ data_umap$template$strain)
    
    umap3d <- umap3d %>% 
      layout(scene = list(
        xaxis = list(title = 'UMAP1'),
        yaxis = list(title = 'UMAP2'),
        zaxis = list(title = 'UMAP3')
      ))
    
    umap3d
  })
  
  # UMAP 3D postEDA
  output$umap.3d.post <- renderPlotly({
    umap3d <- plot_ly(
      data_umap12_postEDA,
      x = ~UMAP1,
      y = ~UMAP2,
      z = ~UMAP3,
      color = ~ data_umap12_postEDA$resistance,
      colors = c('#cf280c', '#1b56f7'),
      symbol = ~ data_umap12_postEDA$prediction,
      symbols = c('circle', 'square', 'diamond')
    )
    
    umap3d <- umap3d %>% 
      add_markers(size = 2, text = ~ data_umap12_postEDA$strain)
    
    umap3d <- umap3d %>% 
      layout(scene = list(
        xaxis = list(title = 'UMAP1'),
        yaxis = list(title = 'UMAP2'),
        zaxis = list(title = 'UMAP3')
      ))
    
    umap3d
  })
  
  # PCA multi-panel
  output$pca <- renderPlot({
    pca_rec <- data_rec %>%
      step_pca(all_numeric_predictors(),
               num_comp = input$pca.comp)

    data_pca <- prep(pca_rec, retain = TRUE)
    
    data_pca$template %>% 
      plot_validation_results() + 
      ggtitle("Principal Component Analysis (HR12)")
  }, height = 1000)
  
  # PCA 3D plotly
  output$pca.3d <- renderPlotly({
    
    pca_rec <- data_rec %>%
      step_pca(all_numeric_predictors(),
               num_comp = 3)
    
    data_pca <- prep(pca_rec, retain = TRUE)
    
    pca3d <- plot_ly(data_pca$template, 
                     x = ~PC1, 
                     y = ~PC2, 
                     z = ~PC3, 
                     color = ~data_pca$template$resistance, 
                     colors = c('#cf280c', '#1b56f7')) 
    
    pca3d <- pca3d %>% 
      add_markers(size = 2,
                  text = ~data_pca$template$strain)
    
    pca3d <- pca3d %>% 
      layout(scene = list(xaxis = list(title = 'PC1'),
                          yaxis = list(title = 'PC2'),
                          zaxis = list(title = 'PC3'))) 
    
    pca3d
    
  })
  
}