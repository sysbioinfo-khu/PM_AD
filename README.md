# PM AD

This repository contains the code and methodology for the analysis of **epigenetic** and **transcriptomic data** to investigate the effects of **prenatal particulate matter exposure** on **placental involvement** in **Atopic Dermatitis (AD)**. The analysis employs multiple approaches, including the integration of bulk transcriptome and DNA methylation data, epigenetic data, single-cell RNA sequencing (scRNA-seq), and spatial transcriptomics.

The code is organized into different sections based on the methods used in the paper. Preprocessed files used in the code can be provided upon request.

## 1. Figures
This section contains the code and data for generating figures presented in the paper. The code covers Figures 1 to 4, while Figures 5 and 6 are related to the scRNA-seq section, and Figure 7 is part of the Spatial Transcriptomics section.

## 2. DNAm_bulk_transcriptome
This section focuses on integrating paired bulk transcriptome and DNA methylation data to generate **triplet data**. The methodology for calculating the **FCER1G score**, as described in Figure 3 of the paper, is outlined here.

* **Triplet Data Generation**  
    1. `1_triplet_analysis`: Create CpG-TF-Target gene pairs based on DMC and DEG data (Remap2022).
    2. `TFBS_triplet`: Identify transcription factor binding sites (TFBS) in regions surrounding CpGs (JASPAR2022).

* **FCER1G Score Calculation**  
    1. `FCER1G_score_input`: Code to retrieve the gene list needed to calculate the FCER1G score.
    2. `Calculate_FCER1G_Score`: Code for computing the FCER1G score.
    3. `Plot`: Code to generate plots related to the FCER1G score.

## 3. Epigenetic Data
In this section, epigenetic data is used to further filter the triplet data.

1. `Convert_hic_to_cool`: Code to convert ENCODE Hi-C data into cool files for running hicFindTADs with HiCExplorer tools.
2. `Filtering_triplet`: Code for filtering triplets using epigenetic data (bedintersect).

## 4. scRNA-seq
This section contains the code used for the analysis of scRNA-seq data in the paper.

### 4.1 Phenotypic Similarity
We identified phenotypic similarities between different cell types across developmental stages.

1. `Conversion_anndata_to_seurat`: Code to convert anndata to a Seurat object for TransferAnchor.
2. `TransferAnchor_and_plot`: Code for running TransferAnchor and generating plots from the paper.

### 4.2 Digital Cytometry (AutogeneS, CIBERSORTx)
Digital cytometry methods, such as **AutogeneS** and **CIBERSORTx**, are used to perform cell-type deconvolution from bulk tissue RNA-seq data.

* **AutogeneS**  
    1. `AutogeneS`: Code for running AutogeneS.

* **CIBERSORTx**  
    1. `Skin_Mixture_file_creation`: Code for preparing input data to run CIBERSORTx.
    2. `Skin_Reference_file_creation`: Code to run CIBERSORTx (requires CIBERSORTx Docker image).

`Proportion_fig_and_test`: Code for generating proportion figures (Figures 4c and 5g) and performing statistical tests.

### 4.3 Conserved Gene Signature
This section identifies **conserved gene signatures** across different conditions. These signatures help identify genes that are maintained across developmental stages.

1. `Conserved_gene_code`: Code for identifying conserved genes shown in Figure 4.
2. `Conserved_gene_score`: Code for performing the analysis and generating plots for Figures 4d, 4f, and 4g.

### 4.4 RNA Velocity
**RNA velocity** analysis is used to estimate the dynamics of gene expression and predict cellular transitions over time.

1. `celldancer`: Code for calculating RNA velocity for each module.
2. `Dynamo_downstream`: Code for analyzing pseudotime, state transitions, and generating associated graphs using Dynamo.
3. `Gene_trend_plot`: Code for plotting gene trends along pseudotime, as shown in Figure 4c.

### 4.5 Network Analysis
A **network-based approach** is used to identify key molecular interactions. **scHumanNet** is utilized to construct cell-type-specific networks, followed by weight optimization to build an appropriately weighted Cell-Specific Network (CSN).

1. `1_weight_optimization`: Code for finding optimized weights for the SCINET algorithm in the process of building a cell-type-specific network.
2. `Construct_cell_type_specific_network`: Code to construct a cell-type-specific network based on the optimized weights.
3. `Conserved_network`: Code to identify conserved networks across different developmental stages.

**Network statistics** such as connectivity and differential hubness can be obtained using scHumanNet.

## 5. Spatial Transcriptome
**Spatial transcriptomics** is used to analyze the spatial distribution of gene expression within tissue sections.

### **Cell2Location**
The **Cell2Location** method is used to map the spatial distribution of cell types within tissue samples, providing insights into how cell localization is altered in AD.

1. `Annotation_cell2location`: Code for inferring cell-type abundance in spatial transcriptome data using scRNA-seq.
2. `Postprocessing_reduced`: Code for examining specific cell gene programs and gene expression in spatial data.
