# Set working directory
setwd("/data_hdd1/yang/submission_data/Triplet_output/")

# Load necessary libraries
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(JASPAR2022)
library(TFBSTools)
library(rtracklayer)

# Load data (Triplet result files generated from the 1_Triplet_analysis.R)
# promoter_data <- read.csv('Promoter.txt')
distal_data <- read.csv('Distal.txt')

# data = promoter_data
data = distal_data


# Options for JASPAR
opts <- list()
opts[["species"]] <- 9606  # Human species
opts[["all_versions"]] <- TRUE
pfm_set <- getMatrixSet(JASPAR2022::JASPAR2022, opts)

# Initialize lists to store p-values and GFF3 information for each row
pvalues_list <- vector("list", nrow(data))
gff3_list <- vector("list", nrow(data))
strand_list <- vector("list", nrow(data))

# Loop through each row of the data
for (i in 1:nrow(data)) {
  # Extract TF name and regionID
  tf_name <- data$TF_symbol[i]
  regionID <- data$regionID[i]
  
  # Extract chromosome and position from regionID column(e.g., chr1:161217049-161217050)
  region_parts <- strsplit(regionID, ":")[[1]]
  chr <- region_parts[1]
  positions <- strsplit(region_parts[2], "-")[[1]]
  start_pos <- as.numeric(positions[1])
  end_pos <- as.numeric(positions[2])
  
  # Find matching TF matrix
  matched_matrices <- pfm_set[sapply(pfm_set, function(m) m@name == tf_name)]
  
  # If no matching TF is found, store NA for p-value and GFF3 information
  if(length(matched_matrices) == 0) {
    pvalues_list[[i]] <- NA
    gff3_list[[i]] <- NA
    next
  }
  
  # Retrieve the TF matrix and convert it to a PWM (Position Weight Matrix)
  matrix_id <- matched_matrices[[1]]@ID
  pfm <- TFBSTools::getMatrixSet(JASPAR2022::JASPAR2022, opts = list(ID = c(matrix_id)))
  pwm <- toPWM(pfm)
  
  # Extract the sequence for the specified region (+- 12 bp from regulatory CpG)
  subject <- getSeq(BSgenome.Hsapiens.UCSC.hg38, chr, start_pos - 12, end_pos + 12)
  
  # Search for TF binding sites in the region
  siteset <- searchSeq(pwm, subject, seqname="seq1", min.score="50%", strand="*")
  
  # Convert the binding sites to GFF3 format
  df <- writeGFF3(siteset)
  df$pval <- pvalues(siteset, type="sampling")[[1]]
  
  # Filter the GFF3 data to find regions that include regulatory CpG
  # 이 부분 + strand 일 경우는 13인데 - strand 일 경우는 어떨지 생각해야 함.
  filtered_df <- df[sapply(1:nrow(df), function(i) {
    start <- df$start[i]
    end <- df$end[i]
    
    # Check if the region includes regulatory CpG (cytosine methylation)
    (start <= 13 && end >= 13)
  }), ]
  
  # Calculate the minimum p-value if available
  if (length(filtered_df$pval) > 0) {
    min_pvalue <- min(unlist(filtered_df$pval), na.rm = TRUE)  # Extract the smallest p-value
    min_index <- which.min(unlist(filtered_df$pval))  # Find the index of the smallest p-value
  } else {
    min_pvalue <- NA
    min_index <- NA
  }
  
  # Store the p-value and GFF3 data in the lists
  pvalues_list[[i]] <- min_pvalue
  
  # Extract GFF3 information based on the index of the smallest p-value
  if (!is.na(min_index)) {
    gff3_data <- filtered_df[min_index, 9]
    gff3_list[[i]] <- gff3_data
  } else {
    gff3_list[[i]] <- NA
  }
  
  # Extract strand information based on the index of the smallest p-value
  if (!is.na(min_index)) {
    strand_data <- filtered_df[min_index, 7]
    strand_list[[i]] <- strand_data
  } else {
    strand_list[[i]] <- NA
  }
}


# Add the p-value and GFF3 information to the original data
data$TF_binding_pvalue <- unlist(pvalues_list)
data$GFF3_info <- unlist(gff3_list)

strand_list <- lapply(strand_list, function(x) if (is.null(x)) NA else x)
data$strand <- unlist(strand_list)


# Save the results to CSV files
# write.csv(data, 'Proximal_triplet_TFBS.txt')
write.csv(data, 'Distal_triplet_TFBS.txt')


