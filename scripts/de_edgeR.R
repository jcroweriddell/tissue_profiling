library(edgeR)    
library(gplots)   
library(calibrate)

# Read data
sample <- read_csv(file = "sample_info_tissues.csv", col_names = TRUE)

# Create DGEList object to store input data 
dge <- FPKMresults_gene %>%
  as.matrix() %>%
  DGEList(counts = FPKMresults_gene, group = sample$experiment)

names(dge)
dge

# Ensure that both the number and order of genes in the gene_anno object be the same as
# that in the read_counts object
identical(row.names(FPKMresults), row.names(geneNumToName))
