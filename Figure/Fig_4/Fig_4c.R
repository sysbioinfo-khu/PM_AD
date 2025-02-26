# Set the working directory
setwd('/data_hdd1/yang/Source_data/Fig_4')

# Load necessary packages
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(reshape2)  
library(readxl)


# Load the data for the proportion plot (test_mean.txt file)
data = read_excel("cibersortx_early_plot.xlsx")
data = as.data.frame(data)
rownames(data) = data$Group
data = data[,c(-1,-2)]

data_subset <- data[rownames(data) %in% c("High PM-healthy", "High PM-AD"),]
data_subset = as.data.frame(t(data_subset))  
data_subset$CellType = rownames(data_subset)

# Convert data to long format (required for ggplot)
data_long <- melt(data_subset, id.vars = "CellType")

# Rename columns for easier reference
colnames(data_long) <- c("CellType", "Time", "Proportion")

# Reorder the "CellType" factor so that "HB" is at the bottom
data_long$CellType <- factor(data_long$CellType, 
                             levels = c('EVT', 'SCT', 'VCT', 'fFB1', 'fFB2', "HB"))

# Define color scheme: "HB" will be red, others will be grey
cell_colors <- ifelse(levels(data_long$CellType) == "HB", "red", "grey")

# Create the plot
ggplot(data_long,
       aes(x = Time, stratum = CellType, alluvium = CellType, y = Proportion, fill = CellType, label = CellType)) +
  geom_flow(stat = "alluvium", lode.guidance = "rightleft", color = "darkgray") +
  geom_stratum(width = 0.5) +
  geom_text(stat = "stratum", aes(label = CellType), size = 5, color = "black") +  # Add labels to the strata
  theme_classic() +
  scale_fill_manual(values = cell_colors) +  # Apply custom color scheme
  labs(x = "Time", y = "Proportion", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





###################----- Statistical analysis ----- ###################  
# Load data for statistical analysis 
data = read_excel("cibersortx_early.xlsx")
data = as.data.frame(data)
data = data[,c(-1)]
# Subset data for "G3" and "G4" groups
tet = data[which(data$Group == 'High PM-AD' | data$Group == 'High PM-healthy'),]

p_values <- c()
for (cell_type in colnames(tet)[-1]) {  
  model <- glm(formula = as.formula(paste(cell_type, "~ Group")), data = tet, family = quasibinomial)
  p_value <- anova(model, test = "LRT")$'Pr(>Chi)'[2]
  p_values <- c(p_values, p_value)
}

p_value_df <- data.frame(CellType = colnames(tet)[-1], PValue = p_values)

# Print the p-value results
print(p_value_df)

