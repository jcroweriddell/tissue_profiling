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

#select eye tissue only
#arrange in descending order FPKM, filter out 'low expression' genes >0.0
ALA_eye <-select(FPKMresults, gene_id, ALA_eye_FPKM) %>% 
  filter(ALA_eye_FPKM >10) %>%
  arrange(desc(ALA_eye_FPKM))
head(ALA_eye)

# plot frequency distribution of FPKM for eye tissue

# MDS or PCA plot to check clustering of tissues
# Heat map across all tissues

