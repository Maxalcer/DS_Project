# Load packages
install.packages("ggraph")
install.packages("igraph")
library(ggraph)
library(igraph)
library(ggrepel)
library(ggplot2)
library(dplyr)

# Data filtering for genes from RNAseq and ATACseq

# Load gene list for RNAseq
rnagenes <- read.csv("../Neural_Network/important_genes/important_genes_RNA.csv")

# Load gene list for ATAC seq
atacgenes <- read.csv("../Neural_Network/important_genes/important_genes_ATAC.csv")

# Sort data frame by 'score' column with increasing order
sorted_rna <- rnagenes[order(-rnagenes$score), ]
sorted_atac <- atacgenes[order(-atacgenes$score), ]


# Group cell_type and sort by score in each group for scRNAseq
sorted_rna <- sorted_rna %>%
  arrange(class, desc(score)) %>%
  ungroup()

# Drop column gene_id from sorted_rna
sorted_rna <- sorted_rna %>% select(-gene_id)

# Group by cell_type and sort by score within each group for ATACseq
sorted_atac <- sorted_atac %>%
  arrange(class, desc(score)) %>%
  ungroup()

# Find common genes
common_genes <- intersect(sorted_atac, sorted_rna)
print(common_genes)

# There are no common genes

# Get top 5 genes for each cell type based on highest 5 scores
top_5_sorted_rna <- sorted_rna %>%
  group_by(class) %>%
  arrange(desc(score)) %>%
  slice_head(n = 5) %>%
  ungroup()

# Print top 5 highest scoring genes in ea class
print(top_5_sorted_rna)

# Get top 5 genes for each cell type based on top 5 scores
top_5_sorted_atac <- sorted_atac %>%
  group_by(class) %>%
  arrange(desc(score)) %>%
  slice_head(n = 5) %>%
  ungroup()

# Print top 5 highest scoring genes in ea class
print(top_5_sorted_atac)




# Biogrid interaction network for scRNAseq genes

# Make vector of file names for BioGRID interaction files
file_names <- c("/networks/BIOGRID-GENE-200216-4.4.235.tab3.txt",
"/networks/BIOGRID-GENE-200224-4.4.235.tab3.txt",
"/networks/BIOGRID-GENE-228288-4.4.235.tab3.txt",
"/networks/BIOGRID-GENE-201484-4.4.235.tab3.txt",
"/networks/BIOGRID-GENE-199546-4.4.235.tab3.txt",
"/networks/BIOGRID-GENE-229928-4.4.235.tab3.txt",
"/networks/BIOGRID-GENE-232372-4.4.235.tab3.txt",
"/networks/BIOGRID-GENE-228314-4.4.235.tab3.txt",
"/networks/BIOGRID-GENE-208227-4.4.235.tab3.txt",
"/networks/BIOGRID-GENE-219748-4.4.235.tab3.txt")

# Make vector of corresponding gene names
gene_names <- c("Hba-x",
                "Hbb-y",
                "Elmo1",
                "Mrc1",
                "Bcl11a",
                "Slc25a21",
                "Ppm1l", 
                "Bmp2k",
                "Vav3",
                "Hdac9")

# Make empty list to store graphs
graphs <- list()

# Loop through ea file name
for (i in seq_along(file_names)) {
  # Read the BioGRID interaction file
  interactions <- read.delim(file_names[i], header = TRUE, sep = "\t")
  
  # Make graph object
  g <- graph_from_data_frame(interactions, directed = FALSE, vertices = NULL)
  
  selected_columns <- interactions[, c("Official.Symbol.Interactor.A", "Official.Symbol.Interactor.B", "Experimental.System")]
  
  g <- graph_from_data_frame(selected_columns[, c("Official.Symbol.Interactor.A", "Official.Symbol.Interactor.B")], directed = FALSE)
  
  selected_columns$Experimental.System <- as.factor(selected_columns$Experimental.System)
  transparency_map <- as.numeric(selected_columns$Experimental.System) / max(as.numeric(selected_columns$Experimental.System))
  
  # Add transparency information to edges
  E(g)$transparency <- transparency_map
  
  # Store graph object in the list
  graphs[[i]] <- g
}

# Deal with graph objects

# Access first graph in the list
graph <- graphs[[8]]

# Visualize the graph using ggraph

# Plot graph using ggraph with 'fr' layout
ggraph(graph, layout = 'fr') + 
  geom_edge_link(color = 'gray', size = 0.5) + 
  geom_node_point(color = 'blue', size = 2) + 
  geom_node_text(aes(label = name), repel = TRUE, size = 4) +
  labs(title = "Bmp2k Entrez Gene ID Interactors") +
  theme_void() +
  theme(plot.title = element_text(size = 15))

# Get node labels from the graph
node_labels <- V(graph)$name

# Print node labels
print(node_labels)


# Biogrid interaction network for scATACseq genes

# Make vector of file names for your BioGRID interaction files
file_names_atac<- c("/networks/BIOGRID-GENE-229895-4.4.235.tab3.txt",
                    "/networks/BIOGRID-GENE-213774-4.4.235.tab3.txt",
                    "/networks/BIOGRID-GENE-198382-4.4.235.tab3.txt",
                    "/networks/BIOGRID-GENE-115175-4.4.235.tab3.txt",
                    "/networks/BIOGRID-GENE-205938-4.4.235.tab3.txt",
                    "/networks/BIOGRID-GENE-234763-4.4.235.tab3.txt",
                    "/networks/BIOGRID-GENE-200965-4.4.235.tab3.txt",
                    "/networks/BIOGRID-GENE-198101-4.4.235.tab3.txt",
                    "/networks/BIOGRID-GENE-208470-4.4.235.tab3.txt")

# Make vector of corresponding gene names
gene_names_atac <- c("Ube2o",
                     "Fars2",
                     "Bpgm",
                     "Zeb2",
                     "Iqgap1",
                     "Irf2bp2",
                     "Klf3",
                     "Ank1",
                     "Arhgap23")

# Make empty list to store graphs
graphs_atac <- list()

# Loop through file names
for (i in seq_along(file_names_atac)) {
  # Read the BioGRID interaction file
  interactions <- read.delim(file_names_atac[i], header = TRUE, sep = "\t")
  
  # Create graph object
  g <- graph_from_data_frame(interactions, directed = FALSE, vertices = NULL)
  
  selected_columns <- interactions[, c("Official.Symbol.Interactor.A", "Official.Symbol.Interactor.B", "Experimental.System")]
  
  g <- graph_from_data_frame(selected_columns[, c("Official.Symbol.Interactor.A", "Official.Symbol.Interactor.B")], directed = FALSE)
  
  selected_columns$Experimental.System <- as.factor(selected_columns$Experimental.System)
  transparency_map <- as.numeric(selected_columns$Experimental.System) / max(as.numeric(selected_columns$Experimental.System))
  
  # Add transparency information to edges
  E(g)$transparency <- transparency_map
  
  # Store the graph object in the list
  graphs_atac[[i]] <- g
}


# Access first graph in the list
graph <- graphs_atac[[8]]


# Plot the graph with uniform edge lengths
# Plot with adjusted size
ggraph(graph, layout = 'fr') + 
  geom_edge_link(color = 'gray', size = 0.5) + 
  geom_node_point(color = 'blue', size = 2) + 
  geom_node_text(aes(label = name), vjust = 1, hjust = 1, size = 3) + 
  theme_void()


# Get node labels from the graph
node_labels <- V(graph)$name

# Print the node labels
print(node_labels)
