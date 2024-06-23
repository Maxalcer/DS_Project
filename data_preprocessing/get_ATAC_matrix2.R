library(Matrix)
library(rlist)
path = "../data/unprocessed/atac_organ_counts/"
files <- read.table(paste(path, "md5.txt", sep = ""))$V2

matrices <- list()

for(file in files[1:3]){
  
  file <- paste(path,"/ATAC_counts",file,".rds", sep = "")
  matrices <- list.append(matrices, readRDS(file))
  
}

all_genes <- c()

for(mat in matrices){
  all_genes <- c(all_genes, rownames(mat))
}

all_genes <- sort(unique(all_genes))

exp_matrices <- list()

c <- 1

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

total_matrix <- cbind(exp_matrices)

saveRDS(total_matrix, "../data/unprocessed/ATAC_counts.rds")