library(edgeR)    
library(gplots)   
library(calibrate)

# Read data
dge <- read.table(file = "FPKMresults_gene.tsv.txt", 
                                     header = TRUE)
sample <- read_csv(file = "sample_info_tissues.csv", col_names = TRUE)

# Create DGEList object to store input data 
dge %<>% column_to_rownames("GeneName") %>%
  as.matrix() %>%
  edgeR::DGEList() %>%
  calcNormFactors()

names(dge)
dge

plot(hclust( d = dist(t(cpm(dge)),  method = "euclidean"),
             method = "ward.D"),
     hang = -1, xlab = "Samples", sub = "Distance=euclidean; Method=ward.D")
