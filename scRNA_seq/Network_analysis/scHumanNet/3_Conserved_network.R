library(igraph)
library(readr)
Pla_early_data <- readRDS("/home/yang/PM_AD/Submission/Figure6/Placenta_early/HB_celltypeNet_list_LLS_V2.rds")
Pla_term_data <- readRDS("/home/yang/PM_AD/Submission/Figure6/Placenta_term/HB_celltypeNet_list_LLS_V2.rds")
Fetal_skin_data <- readRDS("/home/yang/PM_AD/Submission/Figure6/Fetal_skin/Fetal_CelltypeNet_list_LLS.rds")
Adult_data <- readRDS("/home/yang/PM_AD/Submission/Figure6/Adult_AD/celltypeNet_list_LLS.rds")

Module = read_excel("/home/yang/PM_AD/Submission/Table/Supplementary_Tables_S5_net.xlsx")


Pla_early_HB <- Pla_early_data[['HB']]
Pla_term_HB <- Pla_term_data[['vil.Hofb']]

Fetal_skin_Mac <- Fetal_skin_data[['M_fs_Macro']]

Adult_healthy_M2 <- Adult_data[['Healthy_Macro_2']]
Adult_AD_M2 <- Adult_data[['Eczema_Macro_2']]


### Change dataframe into network data
Pla_early_HB_net <- graph.data.frame(Pla_early_HB, directed = F)
Pla_term_HB_net <- graph.data.frame(Pla_term_HB, directed = F)

Fetal_skin_Mac_net <- graph.data.frame(Fetal_skin_Mac, directed = F)

Adult_healthy_M2_net <- graph.data.frame(Adult_healthy_M2, directed = F)
Adult_AD_M2_net <- graph.data.frame(Adult_AD_M2, directed = F)

Module[Module$modularity_class =='A',]$Label



### You can follow the code below to plot a network for only the genes in Module_A
# Define the function to filter the network based on Module_A genes
filter_subnetwork <- function(net, Module) {
  # Get the list of genes in the current network (net)
  net_genes <- unique(c(net[, 1], net[, 2]))
  
  # Get the genes from Module_A
  module_genes <- Module
  
  # Filter the genes in the network that are also part of the module
  module_genes_net <- as.character(module_genes[module_genes %in% net_genes])
  
  # Filter edges that only contain module genes
  filtered_edges <- net[net[, 1] %in% module_genes_net & net[, 2] %in% module_genes_net, ]
  
  # Return the filtered network as an igraph object
  net_sub <- graph.data.frame(filtered_edges, directed = FALSE)
  
  return(net_sub)
}

# Apply the function to filter networks for each dataset
# If you want to constructh network with  Module A => Module[Module$modularity_class =='A',]$Label 
Pla_early_HB_net2 = filter_subnetwork(Pla_early_HB, Module = Module$Label)
Pla_term_HB_net2 = filter_subnetwork(Pla_term_HB, Module = Module$Label)
Fetal_skin_Mac_net2 = filter_subnetwork(Fetal_skin_Mac, Module = Module$Label)
Adult_healthy_M2_net2 = filter_subnetwork(Adult_healthy_M2, Module = Module$Label)
Adult_AD_M2_net2 = filter_subnetwork(Adult_AD_M2, Module = Module$Label)

# Combine the filtered networks into a single network
combined_network <- igraph::union(Pla_early_HB_net2, Pla_term_HB_net2, Fetal_skin_Mac_net2, Adult_AD_M2_net2)
common_nodes <- intersect(V(Pla_early_HB_net2)$name,
                          intersect(V(Pla_term_HB_net2)$name,
                                    intersect(V(Fetal_skin_Mac_net2)$name, V(Adult_AD_M2_net2)$name)))


### Construct conversed network.
conserved_network <- induced_subgraph(combined_network, vids = common_nodes)



######
# Define genes in Module_A and Module_B
Module_A = Module[Module$modularity_class == 'A',]$Label
Module_B = Module[Module$modularity_class == 'B',]$Label

# Get the neighbors of FCER1G in the common network (using vertex names)
fcers1g_neighbors <- neighbors(conserved_network, v = which(V(conserved_network)$name == "FCER1G"))
fcers1g_neighbors_names <- V(conserved_network)$name[fcers1g_neighbors]  # Get the names of the neighbors

# Include FCER1G itself and its neighbors in the subgraph
fcers1g_subgraph <- induced_subgraph(conserved_network, vids = c("FCER1G", fcers1g_neighbors_names))

# Node size setting (degree-based)
node_degree <- igraph::degree(fcers1g_subgraph, mode = "all")
base_size <- 3
scale_factor <- 0.5
V(fcers1g_subgraph)$size <- base_size + (node_degree * scale_factor)

# Set node color: genes in Module_A are pink, genes in Module_B are green, others are sky blue
V(fcers1g_subgraph)$color <- ifelse(V(fcers1g_subgraph)$name %in% Module_A, "pink", 
                                    ifelse(V(fcers1g_subgraph)$name %in% Module_B, "green", "skyblue"))

# Color FCER1G red, even if it is in Module_A or Module_B
V(fcers1g_subgraph)$color[V(fcers1g_subgraph)$name == "FCER1G"] <- "red"

# Layout and plotting
layout <- layout_with_fr(fcers1g_subgraph, niter = 10000, repulserad = vcount(fcers1g_subgraph)^0.5, area = vcount(fcers1g_subgraph)^2)

plot(fcers1g_subgraph, layout = layout,
     vertex.color = V(fcers1g_subgraph)$color,
     vertex.size = V(fcers1g_subgraph)$size,
     vertex.label = V(fcers1g_subgraph)$name,
     vertex.label.cex = 0.4,
     vertex.label.color = "black",
     vertex.label.family = "Helvetica",
     vertex.frame.color = NA,
     main = "FCER1G and its Direct Neighbors in the Common Network")

# Add legend
legend("bottomleft", 
       legend = c("FCER1G", "Module A Genes", "Module B Genes", "Other Genes"), 
       col = c("red", "pink", "green", "skyblue"), 
       pch = 19, 
       bty = "n", 
       title = "Node Type")

dev.off()



center_node <- which(V(fcers1g_subgraph)$name == "FCER1G")

# Apply the star layout with FCER1G as the center node
layout <- layout_as_star(fcers1g_subgraph, center = center_node)

# Plot the graph with the updated color scheme
plot(fcers1g_subgraph, layout = layout,
     vertex.color = V(fcers1g_subgraph)$color,  # Apply the custom colors (pink, skyblue, red)
     vertex.size = V(fcers1g_subgraph)$size,
     vertex.label = V(fcers1g_subgraph)$name,
     vertex.label.cex = 0.4,
     vertex.label.color = "black",
     vertex.label.family = "Helvetica",
     vertex.frame.color = NA,
     main = "FCER1G and its Direct Neighbors in a Star Layout")

# Add legend with updated labels for the node types
legend("bottomleft", 
       legend = c("FCER1G", "Module Genes", "Other Genes"), 
       col = c("red", "pink", "skyblue"), 
       pch = 19, 
       bty = "n", 
       title = "Node Type")

dev.off()



