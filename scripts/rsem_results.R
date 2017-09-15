###### Gene profiling using FPKM for sea snake tissues
## tutorials: http://biocluster.ucr.edu/~rkaundal/workshops/R_mar2016/RNAseq.html#sample-wise-correlation-analysis
# http://www.bioconductor.org/packages/release/bioc/vignettes/gage/inst/doc/RNA-seqWorkflow.pdf
## What is FPKM? https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/ 

# set working directory
getwd()
setwd("C:/Users/L033060262053/Documents/Research projects/Tail_photoreception/tissue_profiling/rsem_results/")
# load in data frame
FPKMresults <- read.delim("rsem_all_together.tsv")
head(FPKMresults)
# data exploration
# load dplyr package
install.packages("dplyr")
library(dplyr)

#select A.laevis tissues
#arrange in descending order of eye FPKM, 
ALA_FPKM <-select(FPKMresults, gene_id, ALA_eye_FPKM, ALAjuv_tailA2_FPKM, ALAjuv_tailB5_FPKM, ALAjuv_body_FPKM, ALA_vno_FPKM, ALA_tailA2_FPKM) %>% 
  arrange(desc(ALA_eye_FPKM))
head(ALA_FPKM)
View(ALA_FPKM)

#filter data by 'top gene expression' for eye tissue, cutoff FPKM = >1000
ALA_topEye <- ALA_FPKM %>% filter(ALA_eye_FPKM >1000)

# Option: select only one ALA column, arrange and filter out 'low expression' genes >10.0 FPKM
#ALA_eye <- select(FPKMresults, gene_id, ALA_eye_FPKM) %>%
  #filter(ALA_eye_FPKM >10) %>%
 # arrange(desc(ALA_eye_FPKM))
#head(ALA_eye)

######## frequency distribution of FPKM for eye tissues #############


######## Heat map across all tissues ##############
#install ggplot and RcolorBrewer
install.packages("ggplot2")
library(ggplot2)

install.packages("RColorBrewer")
library("RColorBrewer")

# transform data into matrix format
rnames <- ALA_topEye[,1] #assign Labels from column 1 to 'rnames'
mat_data <- data.matrix(ALA_topEye[,2:ncol(ALA_topEye)]) #tranform columns 2-5 in to matrix
rownames(mat_data) <- rnames # assign rownames

# customizing and plotting the heatmap

# create colour palette from red to green
my_palette <-colorRampPalette(c("red", "yellow", "green"))(n = 299)

# make heatmap
ALA_heatmap <-heatmap(mat_data, Rowv=NA, Colv = NA, col = cm.colors(256), 
                      scale = 'column', margins = c(5,10))

dev.off()

## 
######## MDS or PCA plot to check clustering of tissues ############
# useful paper explaining data visualisation in R https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4387895/ 
# install bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite()
# download edgeR (Bioconductor packags)
biocLite("limma")
biocLite("edgeR")
# guide to limma package http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf

library(limma)
library(edgeR)   

head(FPKMresults)
# fix counts from float to int
FPKMresultsInt <- read.table("all_counts_as_int.txt", header = TRUE)

read_counts <- FPKMresultsInt
sample_info <- read.csv("sample_info_tissues.csv")
#gene_anno <- 

#check counts
read_counts <- as.matrix(read_counts)
  
# Explore loaded R objects
class(read_counts) # print the data type of the object
dim(read_counts) # print the dimension of the object
head(read_counts) # print the first few (default is 6) rows of the object

sample_info # print the entire object as the data frame is small
#head(gene_anno) # print the first 6 rows of the object

#str(gene_anno) # an alternative way to view the structure of an object
# Create DGEList object to store input data 
dge <- new("DGEList",list(counts= read_counts, group = sample_info$experiment))

# print names of all components in the dge list object
names(dge)

# print a summary of the dge list object
dge

# Ensure that both the number and order of genes in the gene_anno object be the same as
# that in the read_counts object
#identical(row.names(read_counts), row.names(gene_anno))

## Data exploration
##########################################################################################
# Hierachical clustering of samples based on CPM values
plot(hclust( d = dist(t(cpm(dge)),  method = "euclidean"),
             method = "ward.D"),
     hang = -1, xlab = "Samples", sub = "Distance=euclidean; Method=ward.D")

# copy plot to pdf file
dev.copy2pdf(file="hierachical clustering of samples.pdf", height=4, width=6)

# MDS plot of samples
plotMDS(dge, top = 500,
         labels = row.names(sample_info), main = "MDS plot", cex=0.8,
         col = c( "red","blue")[ as.factor(sample_info$species)])

