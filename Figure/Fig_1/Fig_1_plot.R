setwd("/data_hdd1/yang/Source_data/Fig_1/")

library(readr)
fig_1a = read_excel('figure 1 data.xlsx', sheet = 'Fig. 1a')
colnames(fig_1a) = c("time_point",'OR','lower_CI','upper_CI','p value')
fig_1a$group = c(rep("pregnancy", 3), rep("first", 3), rep("second", 3), rep("third", 3))

# pdf('Odd_ratio.pdf')
ggplot(fig_1a, aes(x = time_point, y = OR)) +
  geom_point(size = 3, color = "black") +  
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.2) + 
  scale_y_continuous(limits = c(0, 3)) +  
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  
  facet_wrap(~group, scales = "free_x") +  
  theme_classic() +
  labs(
    title = "Odds Ratio with 95% CI",
    x = "Time point",
    y = "Odds ratio"
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1), 
    strip.text = element_text(size = 12)   
  )

dev.off()




# Fig 1c left plot
pathway_res = read_excel('figure 1 data.xlsx', sheet = 'Fig. 1c')
DMC = pathway_res[c(1, 2, 3)]
DMC = as.data.frame(DMC[complete.cases(DMC[, c(1, 2, 3)]), ])

DMC_generatio = melt(DMC[, c(-3)])
DMC_P = melt(DMC[, c(-2)])

variables = c(DMC_generatio$variable)
Pathway_list = c(DMC_generatio$DMC_pathways)

df <- data.frame(
  Pathways = Pathway_list,
  Cells = variables,
  gene_ratio = DMC_generatio,
  Pvalues = DMC_P
)

df$Pathways <- factor(df$Pathways, levels = rev(unique(df$Pathways)))

# pdf("DMC_pathways_plot.pdf", width = 9, height = 5)
ggplot(df, aes(x = gene_ratio.value, y = Pathways, color = 'red')) +
  geom_point(aes(size = 4, fill = 'red'), shape = 21,
             colour = ifelse(df$Pvalues.value < 0.05, 'black', 'white')) + 
  scale_x_continuous(limits = c(0.2, 1)) +
  labs(x = "Gene_ratio", y = "Pathways") + 
  theme_classic()
dev.off()

# Fig 1c right plot
DEG = pathway_res[c(-1,-2,-3,-4)]
DEG$`Wilcoxon-Rank sum P values` <- -log10(DEG$`Wilcoxon-Rank sum P values`)
DEG = as.data.frame(DEG[, c(-2, -3)])
library(reshape2)

DEG_2 = melt(DEG)

DEG_2$DEG_pathways <- factor(DEG_2$DEG_pathways, levels = rev(unique(DEG_2$DEG_pathways)))

# pdf("DEG_pathway_plot.pdf", width = 9, height = 5)
ggplot(DEG_2, aes(x = value, y = DEG_pathways, color = variable)) +
  geom_point(aes(size = 4, fill = variable), shape = 21,
             colour = ifelse(DEG_2$variable == 'High_PM_AD_mean', 'black', 'white')) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(-1.2, 1.2)) +
  labs(x = "Pathway's mean score") +
  scale_fill_manual(values = c("blue", "orange", "green", "red")) + 
  theme_classic()
dev.off()
