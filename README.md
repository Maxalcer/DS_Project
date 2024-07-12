# Creating the Expression Matrices releated to Haematopoiesis: 

The necessary files are not included in the git because they are to large!

## Getting the scATAC-Seq Expression Matrix:

download OMIX928-21 and OMIX928-24 from https://ngdc.cncb.ac.cn/omix/release/OMIX928 (OMIX928-21 is ~ 30GB large)
put OMIX928-24 into "data/unprocessed" (create folder if necessary)

unpack OMIX928-21 to "data/unprocessed" (this should create a folder atac_organ_fragments with 22 tar.gz files)
Run the following files in this particular order:
1. "data_preprocessing/get_ATAC_matrix.rmd" (will require about 90GB of RAM and a couple of hours time)
2. "data_preprocessing/get_ATAC_matrix2.R" (will require about 50GB of RAM and a few hours time)
3. "data_preprocessing/ATAC_preprocessing.rmd" (will require about 30GB of RAM)

## Getting the scRNA-Seq Expression Matrix:

download gene_annotation.csv, cell_annotation.csv and gene_count_cleaned.RDS from https://oncoscape.v3.sttrcancer.org/atlas.gs.washington.edu.mouse.rna/downloads and put them into "data/unprocessed(gene_count_cleaned.RDS is ~ 10GB)

Run "data_preprocessing/RNA_preprocessing.rmd" (will require about 30GB of RAM)

# Differential Gene Expression Analysis:

Run "DESeq/DESeq_ATAC.rmd" and "DESeq/DESeq_RNA.rmd". The file paths are coherrent with the previous scripts

If explicit MA-Plots are required change the sub-trajectory-names in line 43-47

# UMAP: 

to create the UMAP plots run "UMAP/UMAP_RNA.rmd" and "UMAP/UMAP_ATAC.rmd" The file paths are coherrent with the previous scripts

to change the number of genes used for UMAP change the number in line 36

to change the UMAP parameters change them in the lines 44-49

if UMAP has been run with new parameters or a different number of genes, change the name of the output image

# Neural Network:

to get the traing data for the neural networks, run "data_preprocessing/get_training_data.rmd". To get the scRNA-Seq training data, change ATAC -> RNA in each file path. For a different number of input featureschange the number in line 27 (a change to the number of input features for the neural network is also required in this case)

to configure train and test the neural networks and perform the interpreatability run "Neural_Network/Neural_Network_RNA.ipynb" and "Neural_Network/Neural_Network_ATAC.ipynb" The file paths are coherrent withthe previous scripts

If changes are done to the neural network, change the names of the output files

To plot the learning curves run "Neural_Network/plot_learning_curves.ipynb". If new learning curves where created, change the file names accordingly
    