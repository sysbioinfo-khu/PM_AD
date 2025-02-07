# PM AD

This repository contains the code and methodology for the analysis of **epigenetic** and **transcriptomic data** to investigate placental involvement in **Atopic dermatitis (AD)**. The analysis involves multiple approaches, including bulk transcriptome and DNA methylation data integration, epigenetic data filtering, single-cell RNA sequencing (scRNA-seq), and spatial transcriptomics. 

The code is organized into different sections based on the methods used in the paper. 

## 1. DNAm_bulk_transcriptome

This section focuses on the integration of paired bulk transcriptome data and DNA methylation data to generate **triplet data**. The method for calculating the **FCER1G score**, as presented in Figure 3 of the paper, is outlined here. 

## 2. Epigenetic Data

In this section, epigenetic data are employed to further filter the triplet data.
- ChIP-seq (H3K27ac, H4K4me1, H3K4me3)
  - merge data (samtools)
  - peak calling
  
- Hi-C
  - FindTAD
    
- Filtering and refinement of triplet data using epigenetic signatures
  - bedintersect

## 3. scRNA-seq

In this section, 

### 3.1 Phenotypic Similarity
We identified the phenotypic similarities between different cell types across the developmental stage. 

### 3.2 Digital Cytometry (AutogeneS, CIBERSORTx)
**Digital cytometry** methods such as **AutogeneS** and **CIBERSORTx** are used to perform cell-type deconvolution from bulk tissue RNA-seq data. 


### 3.3 Conserved Gene Signature
This part of the analysis identifies **conserved gene signatures** across different conditions. These gene signatures are crucial for identifying the genes that are maintained throughout the developmental stages

### 3.4 RNA Velocity
**RNA velocity** analysis is employed to estimate the dynamics of gene expression and predict cellular transitions over time. 

### 3.5 Network Analysis
A **network-based approach** is used to identify key molecular interactions. 
We used scHumanNet to construct cell type-specific networks. Subsequently, we performed weight optimization to build an appropriately weighted CSN. The code for this analysis is included in this directory.

## 4. Spatial Transcriptome

**Spatial transcriptomics** is used to analyze the spatial distribution of gene expression within tissue sections. 

### 4.1 Cell2Location
The **Cell2Location** method is used to map the spatial distribution of cell types within tissue samples. This analysis provides insights into how different cell types are localized within the tissue and how this localization may be altered in the context of AD.

### 4.2 Cell-type Specific Expression 

