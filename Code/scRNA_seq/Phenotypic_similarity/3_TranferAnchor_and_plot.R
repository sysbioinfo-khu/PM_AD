# Set working directory and load required libraries
setwd("/home/yang/PM_AD/Submission/Figure_4")

# Load required libraries
library(Seurat)
library(dplyr)
library(reticulate)
library(Matrix)
library(stringr)
library(readr)
library(harmony)
library(reshape)
library(ggplot2)
library(useful)
library(cowplot)

# ------------------------------- Load and Preprocess Data -------------------------------

# Load and preprocess early placenta scRNA-seq data
Placenta_early <- readRDS("/data_hdd1/yang/scRNA-seq/placen_skin/RDS/Placenta_early.RDS")
Placenta_early <- NormalizeData(Placenta_early, normalization.method = "LogNormalize", scale.factor = 10000)
Placenta_early <- FindVariableFeatures(Placenta_early, selection.method = "vst", nfeatures = 2000)

# Load and preprocess healthy late placenta scRNA-seq data
Placenta_term <- readRDS("/home/yang/PM_AD/scRNA-seq/placen_skin/RDS/Placenta_term_all_healthy.RDS")
Placenta_term <- NormalizeData(Placenta_term, normalization.method = "LogNormalize", scale.factor = 10000)
Placenta_term <- FindVariableFeatures(Placenta_term, selection.method = "vst", nfeatures = 2000)

# Load and preprocess fetal skin scRNA-seq data
fetal_skin <- readRDS("/home/yang/PM_AD/scRNA-seq/placen_skin/RDS/Fetal_skin_all.RDS")
fetal_skin <- NormalizeData(fetal_skin, normalization.method = "LogNormalize", scale.factor = 10000)
fetal_skin <- FindVariableFeatures(fetal_skin, selection.method = "vst", nfeatures = 2000)

# Load and preprocess healthy adult skin scRNA-seq data
adult_healthy <- readRDS("/data_hdd1/yang/scRNA-seq/placen_skin/RDS/Adult_healthy_skin_all.RDS")
adult_healthy <- NormalizeData(adult_healthy, normalization.method = "LogNormalize", scale.factor = 10000)
adult_healthy <- FindVariableFeatures(adult_healthy, selection.method = "vst", nfeatures = 2000)

# ------------------------------- Find Phenotypic Similarity Between Datasets -------------------------------

# Query: placenta, Ref: fetal skin
anchors <- FindTransferAnchors(reference = fetal_skin, query = Placenta_early, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = fetal_skin$anno_final, dims = 1:30)
tmp1 <- AddMetaData(Placenta_early, metadata = predictions)
tmp2 <- data.frame(predicted.id = tmp1$predicted.id, final_clustering = tmp1$annotation)

# Query: Early Placenta, Ref: Adult healthy
anchors <- FindTransferAnchors(reference = adult_healthy, query = Placenta_early, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = adult_healthy$final_clustering, dims = 1:30)
tmp1 <- AddMetaData(Placenta_early, metadata = predictions)
tmp2 <- data.frame(predicted.id = tmp1$predicted.id, final_clustering = tmp1$annotation)

# Query: Placenta term, Ref: Adult healthy
anchors <- FindTransferAnchors(reference = adult_healthy, query = Placenta_term, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = adult_healthy$final_clustering, dims = 1:30)
tmp1 <- AddMetaData(Placenta_term, metadata = predictions)
tmp2 <- data.frame(predicted.id = tmp1$predicted.id, final_clustering = tmp1$annotation_merged)

# Query: Placenta term, Ref: Fetal skin
anchors <- FindTransferAnchors(reference = fetal_skin, query = Placenta_term, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = fetal_skin$anno_final, dims = 1:30)
tmp1 <- AddMetaData(Placenta_term, metadata = predictions)
tmp2 <- data.frame(predicted.id = tmp1$predicted.id, final_clustering = tmp1$annotation_merged)

# Query: Fetal skin, Ref: Placenta term
anchors <- FindTransferAnchors(reference = Placenta_term, query = fetal_skin, dims = 1:30)
predictions <- TransferData(anchorset = anchors, refdata = Placenta_term$annotation_merged, dims = 1:30)
tmp1 <- AddMetaData(fetal_skin, metadata = predictions)
tmp2 <- data.frame(predicted.id = tmp1$predicted.id, final_clustering = tmp1$anno_final)







# ------------------------------- Mean Data Calculation Function -------------------------------

#To create the plot, we need to summarize the data by calculating the mean of these values
#Function to calculate the mean of prediction scores by annotation
Mean_data <- function(Predictions, res_df) {
  tmp1 <- base::merge(Predictions, res_df, by = "row.names", all = TRUE)
  tmp2 <- tmp1 %>%
    dplyr::select(-Row.names, -predicted.id, -prediction.score.max) %>%
    dplyr::group_by(Annotation) %>%
    dplyr::summarise(dplyr::across(everything(), mean)) %>%
    tibble::column_to_rownames(var = "Annotation")
  tmp2 <- tmp2[-nrow(tmp2), ]  # Remove last row if necessary
  return(tmp2)
}

# ------------------------------- Process Data for Plotting -------------------------------

# Process data and transpose for plotting
new_df <- data.frame(final_clustering = tmp2[, c(2)])
rownames(new_df) <- rownames(tmp2)
colnames(new_df) <- c("Annotation")
res <- t(Mean_data(predictions, new_df))
res_df <- as.data.frame(res)
res_df$X <- rownames(res_df)

# ------------------------------- Plot Figure 4B -------------------------------
# You can create other plots in the paper by using this code with different results 
# This plot represents the query: early placenta, Ref: healthy adult skin cell.
# Sort data by 'HB' (Hofbauer cells) values 
data_sorted <- res_df[order(-res_df$HB), ]
data_sorted$X_clean <- gsub("^prediction.score\\.", "", data_sorted$X)

# Create barplot
ggplot(data_sorted, aes(x = reorder(X_clean, -HB), y = HB)) + 
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "X", y = "HB") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  )

# ------------------------------- Extended Data Figure 4 -------------------------------

# Melt the data for heatmap
data_melted <- reshape2::melt(res_df, id.vars = "X", variable.name = "CellType", value.name = "Score")
data_melted$X <- factor(data_melted$X, levels = unique(data_melted$X))
data_melted$CellType <- factor(data_melted$CellType, levels = unique(data_melted$CellType))

# Create heatmap plot
ggplot(data_melted, aes(x = CellType, y = X, fill = Score)) +
  geom_tile(color = "white", size = 0.5) + 
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(x = "Cell types of Placenta (early)", y = "Cell types of Adult skin") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank()
  ) +
  coord_fixed()  # Ensure aspect ratio is fixed

