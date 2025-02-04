# Set working directory
setwd("/home/yang/PM_AD/Submission/Figure_4")

# Set directory for the raw data
setwd('/home/yang/PM_AD/scRNA-seq/placen_skin/For_RDS/Placenta_term_all/')

# Load the necessary data
NAME_counts <- readMM("Placenta_term_all_counts.mtx")
NAME_genes <- read.csv("Placenta_term_all_genes.csv")
NAME_metadata <- read.csv('Placenta_term_all_metadata.csv', row.names = 1)
NAME_cells <- read.csv('Placenta_term_all_cell_names.csv', row.names = 1)

# Assign column names to counts data (cell names) and row names (gene names)
colnames(NAME_counts) <- NAME_cells$index
rownames(NAME_counts) <- NAME_genes$gene_short_name

# Create a Seurat object with the counts data and add metadata
NAME <- CreateSeuratObject(counts = NAME_counts, project = "skin", min.cells = 3, min.features = 200)
NAME <- AddMetaData(NAME, NAME_metadata)

# Save the Seurat object as an RDS file for later use
saveRDS(NAME, "/home/yang/PM_AD/scRNA-seq/placen_skin/RDS/Placenta_term_all_healthy.RDS")

# Once the RDS file is created, you can estimate phenotypic similarity using Seurat's TransferAnchor function
