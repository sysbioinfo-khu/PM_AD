setwd("/data_hdd1/yang/Source_data/Fig_2")

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
# Fig. 2a

data = read_excel("Fiugre2A.xlsx")
data = as.data.frame(data)
# triplet = data[data['PM_status'] =='high_PM',]
triplet = data[data['PM_status'] =='low_PM',]

triplet <- data.frame(
  index = triplet$TF_symbol,
  DNAm_low_RLM_target_vs_TF_estimate = triplet$DNAm_low_RLM_target_vs_TF_estimate,
  DNAm_low_RLM_target_vs_TF_pvalue = triplet$DNAm_low_RLM_target_vs_TF_pvalue,
  DNAm_high_RLM_target_vs_TF_estimate =triplet$DNAm_high_RLM_target_vs_TF_estimate,
  DNAm_high_RLM_target_vs_TF_pvalue = triplet$DNAm_high_RLM_target_vs_TF_pvalue
)


triplet_long <- data.frame(
  index = rep(c("Low DNAm", "High DNAm"), times = nrow(triplet)),
  group = rep(triplet$index, each = 2),  
  estimate = unlist(lapply(seq_len(nrow(triplet)), function(i) c(triplet$DNAm_low_RLM_target_vs_TF_estimate[i], triplet$DNAm_high_RLM_target_vs_TF_estimate[i]))),
  pvalue = unlist(lapply(seq_len(nrow(triplet)), function(i) c(triplet$DNAm_low_RLM_target_vs_TF_pvalue[i], triplet$DNAm_high_RLM_target_vs_TF_pvalue[i])))
)


get_color <- function(estimate, pvalue) {
  if (pvalue > 0.05) {
    return("gray")
  } else if (estimate > 0) {
    return("red")
  } else {
    return("blue")
  }
}

triplet_long$color <- mapply(get_color, triplet_long$estimate, triplet_long$pvalue)
triplet_long$group <- factor(triplet_long$group, levels = c("RELB", "NFE2L2", "RUNX1"))
triplet_long$index <- factor(triplet_long$index, levels = c("Low DNAm", "High DNAm"))

ggplot(triplet_long, aes(x = index, y = estimate, fill = color)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  scale_fill_identity() +
  facet_wrap(~group) +
  labs(x = "", y = "RLM Estimate", title = "TF-target gene regulation mediated by CpG methylation") +
  theme_classic()

dev.off()

# Fig. 2b

df = read_excel("Figure_2B.xlsx")
df = as.data.frame(df)
proximal = df[,c(2,3,4,5,13)]
distal = df[,c(7,8,9,10,11,12,13)]

df2 = proximal
# df2 = distal


# High PM
df_filtered <- df2 %>%
  filter(Group %in% c("High PM-healthy", "High PM-AD"))

# Low PM
# df_filtered <- df %>%
#   filter(Group %in% c("Low PM-healthy", "Low PM-AD"))


df_long <- df_filtered %>%
  pivot_longer(cols = -Group, names_to = "Column", values_to = "Value")

#  median, Q1 (1st quartile), Q3 (3rd quartile) 
df_stats <- df_long %>%
  group_by(Group, Column) %>%
  summarise(
    Median = median(Value, na.rm = TRUE),
    Q1 = quantile(Value, 0.25, na.rm = TRUE),
    Q3 = quantile(Value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# x-axis order
df_stats$Column <- factor(df_stats$Column, levels = colnames(df2)[-ncol(df2)])


ggplot(df_stats, aes(x = Column, y = Median, color = Group, group = Group)) +
  geom_line() +  
  geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.2) +  
  labs(x = "", y = "Median Values", title = "", color = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_classic() +
  scale_x_discrete(name = "Columns") +
  scale_y_continuous(name = "Methylation (Median)") +
  theme(legend.position = "top") +
  scale_color_manual(values = c("High PM-healthy" = "orange1", "High PM-AD" = "red"))

dev.off()


# Figure 2b right (gene expression)
gene_exp = df[,c('Group','FCER1G')]

# High PM
gene_exp2 <- gene_exp %>%
  filter(Group %in% c("High PM-healthy", "High PM-AD"))

ggplot(gene_exp2, aes(x = Group, y = FCER1G, fill = Group, color = Group)) +
  geom_boxplot(color = "black") + 
  geom_jitter(width = 0.2, size = 2) +  
  labs(title = "", x = "", y = "FCER1G Expression") + 
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("High PM-healthy" = "orange1", "High PM-AD" = "red")) +
  scale_color_manual(values = c("High PM-healthy" = "orange1", "High PM-AD" = "red"))

dev.off()


# Figure 2c
df = read_excel("Figure_2c.xlsx")
df = as.data.frame(df)
proximal = df[,c(2,3,4,5,11)]
distal = df[,c(7,8,9,10,11)]

df2 = proximal
# df2 = distal


df_long <- df2 %>%
  pivot_longer(cols = -Group, names_to = "Column", values_to = "Value")

#  median, Q1 (1st quartile), Q3 (3rd quartile) 
df_stats <- df_long %>%
  group_by(Group, Column) %>%
  summarise(
    Median = median(Value, na.rm = TRUE),
    Q1 = quantile(Value, 0.25, na.rm = TRUE),
    Q3 = quantile(Value, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# x-axis order
df_stats$Column <- factor(df_stats$Column, levels = colnames(df2)[-ncol(df2)])


ggplot(df_stats, aes(x = Column, y = Median, color = Group, group = Group)) +
  geom_line() +  
  geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.2) +  
  labs(x = "", y = "Median Values", title = "", color = "Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_classic() +
  scale_x_discrete(name = "Columns") +
  scale_y_continuous(name = "Methylation (Median)") +
  theme(legend.position = "top") +
  scale_color_manual(values = c("disease: Control" = "blue", "disease: Atopic dermatitis" = "red"))

dev.off()



# Figure 2d

gene_exp = read_excel("Figure_2d.xlsx")
gene_exp = as.data.frame(gene_exp)

# Adults
gene_exp2 <- gene_exp %>%
  filter(agegroup %in% c("Adults"))

gene_exp2$status <- factor(gene_exp2$status, levels = c("Healthy", "Non-lesional AD","Lesional AD"))

ggplot(gene_exp2, aes(x = status, y = FCER1G, fill = status, color = status)) +
  geom_boxplot(color = "black") +  
  geom_jitter(width = 0.2, size = 2) +  
  labs(title = "", x = "", y = "FCER1G Expression") +  
  theme_classic() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Healthy" ='blue',"Non-lesional AD" = "yellow2", "Lesional AD" = "red")) +
  scale_color_manual(values = c("Healthy" ='blue',"Non-lesional AD" = "yellow2", "Lesional AD" = "red"))

dev.off()

