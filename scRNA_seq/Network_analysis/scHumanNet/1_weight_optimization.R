# Load necessary libraries
library(scHumanNet)
library(ACTIONet)
library(SCINET)
library(Seurat)
library(igraph)
library(SingleCellExperiment)
library(dplyr)
library(readxl)

### Process Placenta Term Data
### Other datasets can be converted into RDS files and this code will apply directly.

# Load the data
data <- readRDS("/data_hdd1/yang/scRNA-seq/placen_skin/RDS/Placenta_term_scHuman.rds")
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

### Specify the cell types for CSN calculation
data$annotation = data$annotation_merged
data <- subset(data, subset = annotation %in% c("vil.EVT", "vil.Ery", "vil.SCT", "vil.VCT", "vil.FB", "vil.Hofb"))
# data <- subset(data, subset = annotation %in% c("EVT", "SCT", "VCT", "fFB1", "fFB2", "HB")) # Early placenta 
data$celltype_condition <- data$annotation_merged
data2 <- SingleCellExperiment(assays = list(logcounts = data@assays$RNA$counts), colData = data@meta.data)

# Compute ACE
ace = reduce.ace(data2)
ace[['Labels']] <- data2$celltype_condition
ace <- compute.cluster.feature.specificity(ace, ace$Labels, "celltype_specificity_scores")



#### Estimate marker genes
Marker_gene_sc = data
Idents(Marker_gene_sc) = Marker_gene_sc$celltype_condition
marker <- FindAllMarkers(Marker_gene_sc, only.pos = TRUE, logfc.threshold = 0.25, min.diff.pct = 0.25, min.pct = 0.25)

gene_list <- split(marker$gene, marker$cluster)

# Jaccard Index Function
jaccard_index <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}


results_combined <- data.frame(
  Min_Edge_Weight = numeric(),
  Type = character(),
  Model1 = character(),
  Model2 = character(),
  Jaccard_Index = numeric(),
  stringsAsFactors = FALSE
)

# Set the range for min.edge.weight
min_edge_weights <- seq(2, 5, by = 0.1)

# Compute Jaccard Index for each min.edge.weight
valid_weights <- c()

for (weight in min_edge_weights) {
  # Run SCINET clustering
  Celltype.specific.networks <- run.SCINET.clusters(
    ace, 
    specificity.slot.name = "celltype_specificity_scores_feature_specificity", 
    thread_no = 40, 
    min.edge.weight = weight
  )
  
  # Extract gene sets for cell type-specific network
  genesets <- lapply(Celltype.specific.networks, function(x) V(x)$name)
  model_names <- names(Celltype.specific.networks)
  
  # Check if all gene sets have more than 100 genes
  if (all(sapply(genesets, length) >= 100)) {
    # Compute Jaccard index if all sets have at least 100 genes
    for (i in 1:length(model_names)) {
      for (j in i:length(model_names)) {
        model1 <- model_names[i]
        model2 <- model_names[j]
        
        common_genes <- length(intersect(genesets[[model1]], genesets[[model2]]))
        total_genes_union <- length(union(genesets[[model1]], genesets[[model2]]))
        
        # Compute Jaccard index
        jaccard_index_value <- common_genes / total_genes_union
        
        # Add to results dataframe
        results_combined <- rbind(results_combined, 
                                  data.frame(Min_Edge_Weight = weight, 
                                             Type = "Celltype Specific Network", 
                                             Model1 = model1, 
                                             Model2 = model2, 
                                             Jaccard_Index = jaccard_index_value))
      }
    }
    
    # Add valid weights
    valid_weights <- c(valid_weights, weight)
    
    # Compute Jaccard index between marker genes and cell type-specific networks
    for (model in model_names) {
      if (model %in% names(gene_list)) {
        jaccard_index_value_marker <- jaccard_index(genesets[[model]], gene_list[[model]])
        
        # Add to results dataframe
        results_combined <- rbind(results_combined,
                                  data.frame(Min_Edge_Weight = weight, 
                                             Type = "Marker Gene Network", 
                                             Model1 = model, 
                                             Model2 = "Marker Genes", 
                                             Jaccard_Index = jaccard_index_value_marker))
      }
    }
  }
}

# Repeating the process for weight range adjustment
start_weight <- 2  
valid_weights <- c()  

while (weight <= 5) {  
  # Run SCINET clustering
  Celltype.specific.networks <- run.SCINET.clusters(
    ace, 
    specificity.slot.name = "celltype_specificity_scores_feature_specificity", 
    thread_no = 40, 
    min.edge.weight = weight
  )
  
  # Extract gene sets for cell type-specific network
  genesets <- lapply(Celltype.specific.networks, function(x) V(x)$name)
  model_names <- names(Celltype.specific.networks)
  
  # Check if all gene sets have more than 100 genes
  if (all(sapply(genesets, length) >= 100)) {
    for (i in 1:length(model_names)) {
      for (j in i:length(model_names)) {
        model1 <- model_names[i]
        model2 <- model_names[j]
        
        common_genes <- length(intersect(genesets[[model1]], genesets[[model2]]))
        total_genes_union <- length(union(genesets[[model1]], genesets[[model2]]))
        
        # Compute Jaccard index
        jaccard_index_value <- common_genes / total_genes_union
        
        # Add to results dataframe
        results_combined <- rbind(results_combined, 
                                  data.frame(Min_Edge_Weight = weight, 
                                             Type = "Celltype Specific Network", 
                                             Model1 = model1, 
                                             Model2 = model2, 
                                             Jaccard_Index = jaccard_index_value))
      }
    }
    
    valid_weights <- c(valid_weights, weight)
    
    # Compute Jaccard index between marker genes and cell type-specific networks
    for (model in model_names) {
      if (model %in% names(gene_list)) {
        jaccard_index_value_marker <- jaccard_index(genesets[[model]], gene_list[[model]])
        
        # Add to results dataframe
        results_combined <- rbind(results_combined,
                                  data.frame(Min_Edge_Weight = weight, 
                                             Type = "Marker Gene Network", 
                                             Model1 = model, 
                                             Model2 = "Marker Genes", 
                                             Jaccard_Index = jaccard_index_value_marker))
      }
    }
  } else {
    break  # Exit loop if a set has fewer than 100 genes
  }
  
  weight <- weight + 0.1
}

valid_weights
results_combined


####### Ranking
Marker_net = results_combined[results_combined$Type == 'Marker Gene Network',]
CCN = results_combined[results_combined$Type == 'Celltype Specific Network',]

# Identify unique edge weights
unique(Marker_net$Min_Edge_Weight)
unique(CCN$Min_Edge_Weight)

##### Marker Gene Network Ranking
pivot_df <- Marker_net %>%
  tidyr::pivot_wider(names_from = Min_Edge_Weight, values_from = Jaccard_Index, 
                     names_prefix = "Weight_")

ranked_df <- pivot_df
# Rank by Weight columns
ranked_df[, 4:ncol(ranked_df)] <- t(apply(pivot_df[, 4:ncol(pivot_df)], 1, function(x) rank(-x, ties.method = "min")))
colnames(ranked_df)[4:ncol(ranked_df)] <- paste0("Rank_", colnames(pivot_df)[4:ncol(pivot_df)])
print(ranked_df)

# Normalize the rank for each cell type (0-1 range) - using Min-Max scaling 

###### Cell Type Specific Network Ranking
pivot_df2 <- CCN %>%
  tidyr::pivot_wider(names_from = Min_Edge_Weight, values_from = Jaccard_Index, 
                     names_prefix = "Weight_")

ranked_df2 <- pivot_df2
ranked_df2[, 4:ncol(ranked_df2)] <- t(apply(pivot_df2[, 4:ncol(pivot_df2)], 1, function(x) rank(x, ties.method = "min")))
colnames(ranked_df2)[4:ncol(ranked_df2)] <- paste0("Rank_", colnames(pivot_df2)[4:ncol(pivot_df2)])
ranked_df2

# Average ranking across models
average_ranked_df <- ranked_df2 %>%
  group_by(Type, Model1) %>% 
  summarise(across(starts_with("Rank_Weight"), mean, na.rm = TRUE), .groups = 'drop')

# Full join to combine Marker and CSN average ranked data
combined_df <- full_join(ranked_df, average_ranked_df, by = "Model1")

# Final dataset for analysis
Marker_ranked_df <- ranked_df %>%
  dplyr::select(-Type, -Model2)

CSN_average_ranked_df = average_ranked_df %>%
  dplyr::select(-Type)

CSN_average_ranked_df_sorted <- CSN_average_ranked_df %>%
  arrange(factor(Model1, levels = Marker_ranked_df$Model1))

result_df <- CSN_average_ranked_df_sorted[,-1] + Marker_ranked_df[,-1]
mean_values <- colMeans(result_df, na.rm = TRUE)
mean_df <- data.frame(
  Rank = names(mean_values),
  Mean_Value = as.numeric(mean_values)
)

#### Dot Plot for Extended Figure 6
dot_plot <- ggplot(mean_df, aes(x = Rank, y = Mean_Value)) +
  geom_point(size = 3) + 
  geom_line(group = 1) + 
  labs(title = "Dot Plot with Line of Mean Values",
       x = "Rank Weights",
       y = "Mean Value") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
dot_plot
dev.off()

# Combine two dataframes
combined_df <- Scaled_average_ranked_df %>%
  full_join(Scaled_average_ranked_df2, by = c("Model1"), suffix = c("_df1", "_df2"))

# Group by Model1 and calculate average
final_df <- combined_df %>%
  group_by(Model1) %>%
  summarise(across(starts_with("Rank_Weight"), mean, na.rm = TRUE), .groups = 'drop')

# Print final dataframe
final_df


##### Heatmap Generation
results_pivot <- CCN[CCN$Min_Edge_Weight == '3.7',] %>%
  dplyr::select(Model1, Model2, Jaccard_Index) %>%
  distinct() %>%
  tidyr::pivot_wider(names_from = Model2, values_from = Jaccard_Index, values_fill = list(Jaccard_Index = NA))

# Fill diagonal with 1
for (model in results_pivot$Model1) {
  results_pivot[results_pivot$Model1 == model, model] <- 1  
}

# Fill NaN values with Jaccard index values
for (model in results_pivot$Model1) {
  for (col in names(results_pivot)[-1]) {  
    if (is.na(results_pivot[results_pivot$Model1 == model, col])) {
      results_pivot[results_pivot$Model1 == model, col] <- results_pivot[results_pivot$Model1 == col, model]
    }
  }
}

# Print results
print(results_pivot)

# Heatmap Setup
library(ComplexHeatmap)
library(dplyr)
library(reshape2) 
library(RColorBrewer)  
library(circlize)

# Convert results_pivot to long format
results_long <- melt(results_pivot, id.vars = "Model1", variable.name = "Model2", value.name = "Jaccard_Index")

# Generate Jaccard index matrix
jaccard_matrix <- as.matrix(results_pivot[, -1])
rownames(jaccard_matrix) <- results_pivot$Model1

# Set quantile-based color scale
breaks <- quantile(jaccard_matrix, probs = seq(0, 1, length.out = 100), na.rm = TRUE)
reds_colors <- colorRamp2(c(0, 1), c("white", "red"))

# Generate Heatmap
heatmap <- Heatmap(jaccard_matrix,
                   name = "Jaccard Index",
                   col = reds_colors,
                   cluster_rows = TRUE,
                   cluster_columns = TRUE,
                   show_row_names = TRUE,
                   show_column_names = TRUE,
                   row_names_gp = gpar(fontsize = 10),
                   column_names_gp = gpar(fontsize = 10),
                   column_title = "CSN",
                   row_title = "CSN",
                   heatmap_legend_param = list(title = "Jaccard Index", at = seq(0, 1, 0.1)))
draw(heatmap)

dev.off()
