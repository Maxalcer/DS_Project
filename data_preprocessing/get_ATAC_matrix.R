library(Matrix)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)

path = "../data/unprocessed/atac_organ_fragments/"
files <- read.table(paste(path, "md5.txt", sep = ""))$V2

genes <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)

for(file in files[20]){
  
  fragments <- read_tsv(paste(path, file, sep = ""), 
                        col_names = FALSE,
                        num_threads = 16,
                        show_col_types = FALSE)
  
  print("file loaded")
  
  atac_Granges <- GRanges(seqnames = fragments$X1,
                          ranges = IRanges(start = fragments$X2, 
                                           end = fragments$X3),
                          sample = fragments$X4)
  rm(fragments)
  
  nearest_genes <- nearest(atac_Granges, genes)
  atac_Granges$gene_id <- mcols(genes)$gene_id[nearest_genes]
  
  print("nearest genes calculated")
  
  gene_symbols <- mapIds(org.Mm.eg.db,
                         keys = atac_Granges$gene_id,
                         column = "SYMBOL",
                         keytype = "ENTREZID",
                         multiVals = "first")
  
  atac_Granges$gene_symbol <- gene_symbols
  
  atac_df <- as_tibble(atac_Granges)
  
  print("got dataframe")
  
  rm(atac_Granges)
  
  atac_df <- atac_df %>%
    summarise(count = n(), .by = c("gene_symbol", "sample")) %>%
    arrange(gene_symbol, sample) %>%
    ungroup()
  
  print(atac_df)
  
  atac_df <- pivot_wider(atac_df, names_from = sample, values_from = count, values_fill = list(count = 0)) %>%
    drop_na(gene_symbol) %>%
    mutate(gene_symbol = as.character(gene_symbol)) %>%
    column_to_rownames(var = "gene_symbol")
  
  print("reshaped dataframe")
  
  atac_matrix <- as(atac_df, "sparseMatrix")
  
  print(paste("got matrix,", file))
  
  saveRDS(atac_matrix, paste("../data/unprocessed/atac_organ_counts/ATAC_counts",file,".rds", sep = ""), compress = TRUE)
  
  rm(atac_df, atac_matrix)
}


