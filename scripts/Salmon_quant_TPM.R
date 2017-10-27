library(magrittr)
library(tibble)
library(readr)
library(stringr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(pheatmap)
library(reshape2)
library(tibble)
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")

# Set working directory and read in data
setwd("/Users/L033060262053/Documents/Research projects/Tail_photoreception/tissue_profiling/salmon_quant/")
refGenome <- read.table(file = "tpm_quant/PMucros_quant_tpm_all.tsv", header = TRUE)
sampleinfo <- read_csv(file = "TPM_quant/sample_info_tissues.csv", col_names = TRUE)

# Change column names to tissue names
names(refGenome) 
refGenome <- select(refGenome, GeneName = Name, 
                  HMAJ_testis = TPM,
                  HMAJ_liver = TPM.1,
                  HMAJ_heart = TPM.2,
                  ALA_vno = TPM.3,
                  ALA_tailA2 = TPM.4,
                  ATEN_tailA2 = TPM.6,
                  ATEN_tailB5 = TPM.7,
                  ATEN_body = TPM.8,
                  ALA_eye = TPM.5,
                  BRH_vno = TPM.9,
                  ALAjuv_body = TPM.12,
                  ALAjuv_tailA2 = TPM.10,
                  ALAjuv_tailB5 = TPM.11,
                  NSC_vno = TPM.13)


# Convert TPM to FPKM
# Load TPMToFpkm function
TPMToFpkm <- function(TPM, geneLength){
  count <- apply(TPM, 2, function(e){ e/1000*geneLength })
  fpkm <- apply(count, 2, function(e){ e*1000*1000000/sum(e)/geneLength })
  return(fpkm)
}

# Subset refGenome data object
geneLength <- refGenome$geneLength ## Doesn't matter which geneLength column you choose. They're all the SAME!!!!!

## Getting TPM dataframe
TPM <- refGenome %>%
  column_to_rownames("GeneName") %>%
  select(-contains("geneLength"))

## Generating FPKM values
Fpkm <- TPMToFpkm(TPM = TPM, geneLength = geneLength) %>% 
  data.frame() %>% rownames_to_column("GeneName")

# Clustering analysis
dgerefGenome <- select(refGenome, GeneName, HMAJ_testis, HMAJ_liver, HMAJ_heart, ALA_tailA2, ATEN_tailB5, ATEN_tailA2,
                     ATEN_body, ALA_eye, ALAjuv_body, ALAjuv_tailA2, ALAjuv_tailB5)

dgeSkin <- select(refGenome, GeneName,ATEN_tailB5, ALA_tailA2, ATEN_tailA2,
                  ATEN_body, ALAjuv_body, ALAjuv_tailA2, ALAjuv_tailB5 )

## Turning counts object into DGElist and normalising
dgeSkin %<>% column_to_rownames("GeneName") %>%
  as.matrix() %>%
  edgeR::DGEList() %>%
  edgeR::calcNormFactors()

## Obtaining MDS data object 
mds <- limma::plotMDS(dgeSkin, gene.selection = "common", method = "logFC")

## Slightly ganky way of doing this, but extracting MDS coordinate object and joining with sample information.
plot <- mds@.Data[[3]] %>%
  as.data.frame() %>%
  rownames_to_column("group") %>%
  set_colnames(c("group","x","y")) %>%
  left_join(sampleinfo, by = c("group"="gene_id"))

## Plotting MDS using ggplot (nicer plot), includes overlap of experimental condition and species
#pdf(file = "plots/MDS_refGenome.pdf")
plot %>%
  ggplot(aes(x, y, colour = species, shape = species, label = tissue)) +
  geom_point(size = 4) +
  geom_text(hjust = 0, nudge_x = 0.2, check_overlap = TRUE) + 
  theme_bw(base_size = 16) +
  theme(legend.text = element_text(size = 12)) + 
  labs(x = "Dimension 1",
       y=  "Dimension 2")
#dev.off()

## Sample clustering
pdf(file = "plots/hclust_refGenome.pdf")
plot(hclust( d = dist(t(edgeR::cpm(dgerefGenome)),  method = "euclidean"),
             method = "ward.D"),
     hang = -1, xlab = "Samples", sub = "Distance=euclidean; Method=ward.D")
dev.off()




