library(Matrix)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(tidyverse)

path = "../data/unprocessed/atac_organ_fragments/"
files <- read.table(paste(path, "md5.txt", sep = ""))$V2

genes <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)

all_genes <- c()

for(file in files){
  fragments <- read_tsv(paste(path, file, sep = ""), 
                        col_names = FALSE,
                        num_threads = 4)
  
  atac_Granges <- GRanges(seqnames = fragments$X1,
                          ranges = IRanges(start = fragments$X2, 
                                           end = fragments$X3),
                          sample = fragments$X4)
  rm(fragments)
  
  nearest_genes <- nearest(atac_Granges, genes)
  atac_Granges$gene_id <- mcols(genes)$gene_id[nearest_genes]
  
  
  gene_symbols <- mapIds(org.Mm.eg.db,
                         keys = atac_Granges$gene_id,
                         column = "SYMBOL",
                         keytype = "ENTREZID",
                         multiVals = "first")
  
  atac_Granges$gene_symbol <- gene_symbols
  
  atac_df <- as_tibble(atac_Granges)
  
  rm(atac_Granges)
  
  atac_df <- atac_df %>%
    group_by(gene_symbol, sample) %>%
    summarise(count = n()) %>%
    ungroup()
  
  atac_df <- pivot_wider(atac_df, names_from = sample, values_from = count, values_fill = list(count = 0)) %>%
    drop_na(gene_symbol) %>%
    mutate(gene_symbol = as.character(gene_symbol)) %>%
    column_to_rownames(var = "gene_symbol")
  
  all_genes <- unique(c(all_genes, row.names(atac_df)))
  
  atac_matrix <- as(atac_df, "sparseMatrix")
  rm(atac_df)
  if(file == files[1]){
    total_matrix <- atac_matrix
  }
  else{
    temp1 <- Matrix(0, 
                    nrow = length(all_genes), 
                    ncol = ncol(atac_matrix),
                    dimnames = list(all_genes, colnames(atac_matrix)),
                    sparse = TRUE)
    
    temp2 <- Matrix(0, 
                    nrow = length(all_genes), 
                    ncol = ncol(total_matrix),
                    dimnames = list(all_genes, colnames(total_matrix)),
                    sparse = TRUE)
    
    temp1[rownames(atac_matrix),] <- atac_matrix
    rm(atac_matrix)
    temp2[rownames(total_matrix),] <- total_matrix
    
    total_matrix <- cbind(temp1, temp2)
    rm(temp1, temp2)
  }
}

saveRDS(total_matrix, "..data/unprocessed/ATAC_counts.rds")
