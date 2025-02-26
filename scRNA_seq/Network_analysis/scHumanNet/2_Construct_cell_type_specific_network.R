library(scHumanNet)
library(ACTIONet)
library(SCINET)
library(Seurat)
library(igraph)
library(SingleCellExperiment)
library(dplyr)
library(readxl)


# Load the data
data <- readRDS("/data_hdd1/yang/scRNA-seq/placen_skin/RDS/Placenta_term_scHuman.rds")
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

### Specify the cell types for CSN calculation
data$annotation = data$annotation_merged
data <- subset(data, subset = annotation %in% c("vil.EVT", "vil.Ery", "vil.SCT", "vil.VCT", "vil.FB", "vil.Hofb")) # Term placenta
# data <- subset(data, subset = annotation %in% c("EVT", "SCT", "VCT", "fFB1", "fFB2", "HB")) # Early placenta 
data$celltype_condition <- data$annotation_merged
data2 <- SingleCellExperiment(assays = list(logcounts = data@assays$RNA$counts), colData = data@meta.data)

# Compute ACE
ace = reduce.ace(data2)
ace[['Labels']] <- data2$celltype_condition
ace <- compute.cluster.feature.specificity(ace, ace$Labels, "celltype_specificity_scores")


#### Construct cell type specific network.
Celltype.specific.networks = run.SCINET.clusters(ace, specificity.slot.name = "celltype_specificity_scores_feature_specificity", thread_no = 40, 
                                                 min.edge.weight = 3.7) # min.edge.weight 3.7 is determined by the 1_weight_optimization process.

sorted.net.list <- SortAddLLS(Celltype.specific.networks, reference.network = graph.hn3)



