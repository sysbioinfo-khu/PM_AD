# Set working directory

# Load the data
data <- readRDS('/home/somi/AD_spatial/sp_data/hpc/rds_OUT/Clustered_8_Visium_LogNorm.rds')
integrated_8_visium = data

# Adjust seurat_clusters for specific cases
integrated_visium_score_data$seurat_clusters[integrated_visium_score_data$seurat_clusters == 7] <- 6

# Function to convert p-values to significance labels
convert_p_value <- function(p) {
  if (is.na(p)) {
    return("NS")  # If p-value is NA, return "NS"
  }
  if (p < 0.0001) {
    return("***")  # p-value < 0.0001, return "***"
  } else if (p < 0.001) {
    return("**")  # p-value < 0.001, return "**"
  } else if (p < 0.05) {
    return("*")  # p-value < 0.05, return "*"
  } else {
    return("NS")  # p-value >= 0.05, return "NS"
  }
}

# Load necessary libraries
library(ggplot2)
library(gridExtra)

# Violin Plot creation function (updated)
plot_violin_for_clusters <- function(object, feature, clusters = 0:6, custom_colors = NULL) {
  
  # Default color settings (if no custom colors are provided)
  if (is.null(custom_colors)) {
    custom_colors <- c("HC" = "#26638C", "NL" = "#A19858", "LS" = "#80434D")
  }
  
  # Save the plot as a PNG file with a 7x1 grid
  png("ViolinPlot_Clusters_7x1.png", width = 3000, height = 18000, res = 300)  # DPI set to 300
  
  # Initialize a list to store the plots
  plot_list <- list()
  
  # Loop over the specified clusters to create the plots
  for (cluster_id in clusters) {
    
    # Filter the data for the current cluster
    cluster_data <- object[, object@meta.data$seurat_clusters == as.character(cluster_id)]
    
    # Perform DEG analysis (comparisons between groups)
    deg_results_LS_vs_HC <- FindMarkers(cluster_data, 
                                        ident.1 = "LS", ident.2 = "HC", 
                                        test.use = "wilcox", 
                                        group.by = 'disease', 
                                        logfc.threshold = 0.0001)
    deg_results_LS_vs_NL <- FindMarkers(cluster_data, 
                                        ident.1 = "LS", ident.2 = "NL", 
                                        test.use = "wilcox", 
                                        group.by = 'disease', 
                                        logfc.threshold = 0.0001)
    
    # Extract and convert p-values
    p_value_LS_vs_HC <- convert_p_value(deg_results_LS_vs_HC[feature, "p_val"])
    p_value_LS_vs_NL <- convert_p_value(deg_results_LS_vs_NL[feature, "p_val"])
    
    # Create the Violin Plot
    p <- VlnPlot(object = cluster_data, 
                 features = feature, 
                 group.by = 'disease', 
                 layer = 'data', alpha = 0.5) + 
      scale_x_discrete(limits = c("HC", "NL", "LS")) + 
      scale_fill_manual(values = custom_colors) + 
      theme_classic() +
      geom_signif(comparisons = list(c("HC", "LS"), c("NL", "LS")), 
                  annotations = c(p_value_LS_vs_HC, p_value_LS_vs_NL), 
                  step_increase = 0.1) + 
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # Automatic scaling with additional 10% space
      theme(
        axis.text.x = element_blank(),  # Remove x-axis text
        axis.title.x = element_blank(),  # Remove x-axis title
        axis.text.y = element_text(size = 12),  # y-axis text size
        axis.title.y = element_blank(),  # Remove y-axis title
        plot.title = element_text(size = 16)  # Title size
      ) +
      ggtitle(paste("Cluster", cluster_id))  # Dynamic title for each cluster
    
    # Add the plot to the list
    plot_list[[cluster_id + 1]] <- p
  }
  
  # Arrange the plots in a 7x1 grid
  grid.arrange(grobs = plot_list, ncol = 1)
  
  # Close the PNG device
  dev.off()
}

# Example usage: Create violin plots for the 'FCER1G' gene across clusters 0 to 6
plot_violin_for_clusters(integrated_8_visium, feature = 'FCER1G', clusters = 0:6)

# Another example to create violin plots for a different gene ('MS4A4A')
plot_violin_for_clusters(integrated_8_visium, feature = 'MS4A4A', clusters = 0:6)







