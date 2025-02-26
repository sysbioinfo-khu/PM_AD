# Set working directory
setwd("/data_hdd1/yang/submission_data/Triplet_input/")

# Load necessary libraries
library(devtools)
library(TFBSTools)
library(ReMapEnrich)
library(MethReg)
library(dplyr)


# Read data
Pla_PM1AD1_DMR <- read.csv("Placenta_PM1AD0_vs_PM1AD1_DMR_with_covariates.txt", row.names = 1)
Pla_PM1AD1_DMR <- Pla_PM1AD1_DMR[,-ncol(Pla_PM1AD1_DMR)] 
Pla_PM1AD1_DMR <- as.data.frame(t(Pla_PM1AD1_DMR)) 

# Create DNAm SummarizedExperiment object
Pla_PM1AD1_DMR.se <- make_dnam_se(
  dnam = Pla_PM1AD1_DMR,
  genome = "hg38",
  arrayType = "EPIC",
  betaToM = FALSE, 
  verbose = FALSE
)

# Load Remap data
remapCatalog2 <- bedToGranges("remap2022_nr_macs2_hg38_v1_0.bed")
gene_names <- sub(":.*", "", remapCatalog2$id)
remapCatalog2$id <- gene_names

create_triplet_distance_based()

# Create triplet data based on Remap2022 data
Proximal_triplet <- create_triplet_distance_based(
  region = Pla_PM1AD1_DMR.se,
  target.method = "genes.promoter.overlap",
  genome = "hg38",
  TF.peaks.gr = remapCatalog2,
  target.promoter.upstream.dist.tss = 2000,
  target.promoter.downstream.dist.tss = 2000,
  motif.search.window.size = 400,
  cores = 1
)



# Distal triplet
Distal_triplet <- create_triplet_distance_based(
  region = Pla_PM1AD1_DMR.se,
  genome = "hg38", 
  target.method = "window",
  target.window.size = 500 * 10^3,
  target.rm.promoter.regions.from.distal.linking = TRUE,
  motif.search.window.size = 400,
  TF.peaks.gr = remapCatalog2,
  cores = 10
)


regulons.dorothea <- dorothea::dorothea_hs

# Read RNA-seq data (DEG)
PM_pla_genes_wilcox <- read.csv("Placenta_PM1AD0_vs_PM1AD1_DEG_with_covariates.txt", check.names = FALSE, row.names = 1)
PM_pla_genes_wilcox <- PM_pla_genes_wilcox[,-25]
Pla_Ensemble2 <- PM_pla_genes_wilcox

# Read RNA-seq data (all genes) for TF activity inference
Pla_all_genes <- read.csv("Placenta_transriptome_cov_PM1AD0_PM1AD1.txt", check.names = FALSE, row.names = 1)
ENSG_cols <- data.frame(t(Pla_all_genes[, grepl("^ENSG", colnames(Pla_all_genes))]), check.names = FALSE)

# Create transcriptome data (SummarizedExperiment object)
Pla_gene.exp.chr21.se <- make_exp_se(
  exp = as.data.frame(Pla_Ensemble2),
  genome = "hg38",
  verbose = FALSE
)

Pla_gene.exp.chr21.se2 <- make_exp_se(
  exp = ENSG_cols,
  genome = "hg38",
  verbose = FALSE
)

# Calculate TF activity using dorothea.
Pla_rnaseq.tf.es <- get_tf_ES(
  exp = Pla_gene.exp.chr21.se2 %>% SummarizedExperiment::assay(),
  regulons = regulons.dorothea
)

# Interaction model
# Proximal_triplet
# Distal_triplet
results.interaction.model <- interaction_model(
  triplet = Distal_triplet, # Proximal_triplet or Distal_triplet
  dnam = Pla_PM1AD1_DMR.se,
  exp = Pla_gene.exp.chr21.se,
  dnam.group.threshold = 0.5,
  sig.threshold = 0.05,
  fdr = FALSE,
  stage.wise.analysis = FALSE,
  filter.correlated.tf.exp.dnam = FALSE,
  filter.triplet.by.sig.term = FALSE, 
  cores = 50,
  tf.activity.es = Pla_rnaseq.tf.es
)


# Stratified model => Assign the role of TF.
results.stratified.model <- stratified_model(
  triplet = results.interaction.model,
  dnam = Pla_PM1AD1_DMR.se,
  exp = Pla_gene.exp.chr21.se2,
  dnam.group.threshold = 0.5,
  tf.dnam.classifier.pval.thld = 0.05,
  tf.activity.es = Pla_rnaseq.tf.es,
  cores = 50
)


# Filter results
filtered_data <- results.stratified.model[!(results.stratified.model$DNAm.effect == 'ns' | is.na(results.stratified.model$DNAm.effect)), ]

# save the results
# write.csv(filtered_data, "/data_hdd1/yang/submission_data/Triplet_output/Promoter.txt")
write.csv(filtered_data, "/data_hdd1/yang/submission_data/Triplet_output/Distal.txt")
