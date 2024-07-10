library(Matrix)
library(rlist)

path = "../data/unprocessed/atac_organ_counts/"
files <- read.table(paste(path, "md5.txt", sep = ""))$V2

matrices <- list()

# create list of expression matrices from the organs (not included in git)
for(file in files){
  file <- paste(path,"/ATAC_counts",file,".rds", sep = "")
  matrices <- list.append(matrices, readRDS(file))
}

all_genes <- c()

# create list of all genes
for(mat in matrices){
  all_genes <- c(all_genes, rownames(mat))
}

all_genes <- sort(unique(all_genes))

exp_matrices <- list()

c <- 1

# extend matrices so that each matrix has the same rows (genes)
for(mat in matrices){
  temp <- Matrix(0, 
                 nrow = length(all_genes),
                 ncol = ncol(mat),
                 dimnames = list(all_genes, colnames(mat)))
  
  temp[rownames(mat),] <- mat
  
  exp_matrices <- list.append(exp_matrices, temp)
  print(paste("copied list", c, "of", length(matrices)))
  c <- c+1
}

rm(matrices)

# bind all matrices together
total_matrix <- do.call(cbind, exp_matrices)

saveRDS(total_matrix, "../data/unprocessed/ATAC_counts.rds")