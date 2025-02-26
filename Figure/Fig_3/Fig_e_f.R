setwd("/data_hdd1/yang/Source_data/Fig_3")

library(ggpubr)
library(reshape2)
library(readxl)
library(ggplot2)
library(ggrepel)

data = read_excel("Fig_3e_statistic.xlsx")
data = as.data.frame(data)
colnames(data)[1] ='Cell_Type'


Clinical_CpG_corr_P = data[,c('P_value_High PM_AD_vs_High PM_healthy','P_value_Low PM_AD_vs_Low PM_healthy')]
Clinical_CpG_corr_R = data[,c('Statistic_High PM_AD_vs_High PM_healthy','Statistic_Low PM_AD_vs_Low PM_healthy')]

Clinical_CpG_corr_P <- melt(Clinical_CpG_corr_P)
Clinical_CpG_corr_R <- melt(Clinical_CpG_corr_R)


pvals = c(Clinical_CpG_corr_P$value)
coefs = c(Clinical_CpG_corr_R$value)


variables = c(Clinical_CpG_corr_P$variable)
Celllist = c(data$Cell_Type)




df <- data.frame(
  Celltype = Celllist,
  Cells = variables,
  Coefficients = coefs,
  Pvalues = pvals
)


df$Celltype <- factor(df$Celltype, levels = rev(unique(data$Cell_Type)))

# pdf("xCell_gene_corr_figure.pdf", width = 5, height = 5)
p <- ggplot(df, aes(x = Cells, y = Celltype)) + theme_classic() +
  geom_point(aes(size = -log10(Pvalues), fill = Coefficients), shape = 21,
             colour = ifelse(df$Pvalues < 0.05,'black','white')) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = 'grey')+
  scale_size(range = c(2, 5)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x=element_blank())

p
dev.off()



### Fig. 3f

data = read_excel('Fig_3f_statistic.xlsx')
data = as.data.frame(data)
colnames(data)[1] ='Cell_Type'

# Add a column for significance (-log10(P_value))
data$Significance <- -log10(data$P_value_Spearman)

# Create the plot
ggplot(data, aes(x = Correlation_R_Spearman, y = Significance, size = Significance, color = Correlation_R_Spearman)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text_repel(aes(label = Cell_Type), 
                  size = 3, color = "black", max.overlaps = Inf, box.padding = 0.5) +
  labs(
    x = "Correlation of cell types with FCER1G",
    y = "Statistical significance (-log10(P))",
    color = "R value",
    size = "-log10(P)"
  ) +
  theme_classic()


dev.off()