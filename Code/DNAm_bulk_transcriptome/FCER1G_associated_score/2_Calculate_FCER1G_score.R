# Load necessary libraries


library(dplyr)
library(GSEABase)
library(GSVA)


# Read the transcriptome data
ad <- read.table("Placenta_transcriptome_data.txt", row.names = 1, header = TRUE, check.names = FALSE)
ad <- data.matrix(ad)

# Read FCER1G-directed genes from CSV file
FCER1G_directed_genes <- c('FCER1G', read.csv('FCER1G_score_input')$X)

# Create a GeneSet for FCER1G directed genes
gs1 <- GeneSet(setName = "FCER1G_directed_genes", FCER1G_directed_genes)



# Create a GeneSetCollection
gsc <- GeneSetCollection(gs1)

# Set up the GSVA parameter
gsvaPar <- gsvaParam(ad, gsc)

# Perform GSVA (Gene Set Variation Analysis)
gsva.es <- gsva(gsvaPar, verbose = FALSE)

# Save the GSVA results to a CSV file
write.csv(gsva.es, 'Figure3B_score.txt')