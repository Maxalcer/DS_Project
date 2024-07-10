Creating the Expression Matrices releated to Haematopoiesis: \n
The necessary files are not included in the git because they are to large! \n
\n
    Getting the scATAC-Seq Expression Matrix:\n
    download OMIX928-21 and OMIX928-24 from https://ngdc.cncb.ac.cn/omix/release/OMIX928 (OMIX928-21 is ~ 30GB large)\n
    put OMIX928-24 into "data/unprocessed" (create folder if necessary)\n
    unpack OMIX928-21 to "data/unprocessed" (this should create a folder atac_organ_fragments with 22 tar.gz files)\n
    Run the following files in this particular order:\n
    1. "data_preprocessing/get_ATAC_matrix.rmd" (will require about 90GB of RAM and a couple of hours time)\n
    2. "data_preprocessing/get_ATAC_matrix2.R" (will require about 50GB of RAM and a few hours time)\n
    3. "data_preprocessing/ATAC_preprocessing.rmd" (will require about 30GB of RAM)\n
\n
    Getting the scRNA-Seq Expression Matrix:\n
    download gene_annotation.csv, cell_annotation.csv and gene_count_cleaned.RDS from https://oncoscape.v3.sttrcancer.org/atlas.gs.washington.edu.mouse.rna/downloads and put them into "data/unprocessed" \n(gene_count_cleaned.RDS is ~ 10GB)\n
    Run "data_preprocessing/RNA_preprocessing.rmd" (will require about 30GB of RAM)\n
\n
Differential Gene Expression Analysis:\n
    Run "DESeq/DESeq_ATAC.rmd" and "DESeq/DESeq_RNA.rmd". The file paths are coherrent with the previous scripts\n
    If explicit MA-Plots are required change the sub-trajectory-names in line 43-47\n
\n
UMAP: \n
    to create the UMAP plots run "UMAP/UMAP_RNA.rmd" and "UMAP/UMAP_ATAC.rmd" The file paths are coherrent with the previous scripts\n
    to change the number of genes used for UMAP change the number in line 36\n
    to change the UMAP parameters change them in the lines 44-49\n
    if UMAP has been run with new parameters or a different number of genes, change the name of the output image\n
\n
Neural Network:\n
    to get the traing data for the neural networks, run "data_preprocessing/get_training_data.rmd". To get the scRNA-Seq training data, change ATAC -> RNA in each file path. For a different number of input features change the number in line 27 (a change to the number of input features for the neural network is also required in this case)\n
    to configure train and test the neural networks and perform the interpreatability run "Neural_Network/Neural_Network_RNA.ipynb" and "Neural_Network/Neural_Network_ATAC.ipynb" The file paths are coherrent with the previous scripts\n
    If changes are done to the neural network, change the names of the output files\n
    To plot the learning curves run "Neural_Network/plot_learning_curves.ipynb". If new learning curves where created, change the file names accordingly\n
    