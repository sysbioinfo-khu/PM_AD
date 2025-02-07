### PM AD: Epigenetic and Transcriptomic Analysis of Placental Data in Alzheimer's Disease
This repository contains the code and methodology for the analysis of epigenetic and transcriptomic data to investigate placental involvement in Alzheimer's Disease (AD). The analysis involves multiple approaches, including bulk transcriptome and DNA methylation data integration, epigenetic data filtering, single-cell RNA sequencing (scRNA-seq), and spatial transcriptomics. These methods are designed to generate valuable insights into the molecular mechanisms underlying AD and to identify potential biomarkers, such as the FCER1G gene, as explored in Figure 3 of our paper.

The code is organized into different sections based on the analysis methods used in the paper. Each section corresponds to a specific approach or analytical step taken in our research.

1. DNAm_bulk_transcriptome
This section focuses on the integration of paired bulk transcriptome data and DNA methylation data to generate triplet data. The method for calculating the FCER1G score, as presented in Figure 3 of the paper, is outlined here. This approach provides insights into the relationship between DNA methylation patterns and gene expression in the context of AD.

Key Steps:
Integration of transcriptome and DNA methylation data
Generation of triplet data (Transcriptome, DNA Methylation, and Gene Expression)
Computation of FCER1G score based on combined data
2. Epigenetic Data
In this section, epigenetic data are employed to further filter the triplet data. This additional filtering step enhances the robustness of the analysis and helps identify more precise molecular signatures relevant to AD pathogenesis. This analysis includes the use of DNA methylation patterns and their potential roles in modulating gene expression.

Key Steps:
Filtering and refinement of triplet data using epigenetic signatures
Analysis of DNA methylation's impact on gene expression
3. scRNA-seq
Single-cell RNA sequencing (scRNA-seq) provides a detailed view of gene expression at the single-cell level. The following sub-sections describe how scRNA-seq data were analyzed to explore various aspects of cellular behavior and gene expression in AD.

3.1 Phenotypic Similarity
This analysis explores the phenotypic similarities between different cellular populations in the context of AD. Using scRNA-seq data, the method identifies key cell types and their relationships, providing insights into the cellular heterogeneity of AD.

3.2 Digital Cytometry (AutogeneS, CIBERSORTx)
Digital cytometry methods such as AutogeneS and CIBERSORTx are used to perform cell-type deconvolution from bulk tissue RNA-seq data. This allows the identification of relative cell-type abundances, shedding light on the cellular composition in AD-affected tissues.

3.3 Conserved Gene Signature
This part of the analysis identifies conserved gene signatures across different conditions. These gene signatures are crucial for understanding the molecular pathways shared across AD samples, regardless of environmental or genetic differences.

3.4 RNA Velocity
RNA velocity analysis is employed to estimate the dynamics of gene expression and predict cellular transitions over time. This provides insights into how cells are transitioning in response to disease-related processes in AD.

3.5 Network Analysis
A network-based approach is used to identify key molecular interactions and pathways involved in AD. By constructing gene co-expression networks, we identify pivotal genes and regulatory relationships that could serve as potential biomarkers or therapeutic targets for AD.

4. Spatial Transcriptome
Spatial transcriptomics is used to analyze the spatial distribution of gene expression within tissue sections. This section focuses on two primary methods for understanding the spatial aspect of gene expression in AD.

4.1 Cell2Location
The Cell2Location method is used to map the spatial distribution of cell types within tissue samples. This analysis provides insights into how different cell types are localized within the tissue and how this localization may be altered in the context of AD.

4.2 Cell-type Specific Expression
In this analysis, cell-type-specific gene expression patterns are explored to better understand how specific cell types contribute to the pathophysiology of AD. The approach helps identify cellular markers or specific gene expression changes in response to AD progression.
   
