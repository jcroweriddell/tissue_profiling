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
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")

# Set working directory and read in data
setwd("/Users/jennacrowe-riddell/Documents/tissue_profiling/salmon_quant/tpm_quant/")
PMucros <- read.table(file = "PMucros_quant_tpm_all.tsv", header = TRUE)
sampleinfo <- read_csv(file = "sample_info_tissues.csv", col_names = TRUE)

# Change column names to tissue names
names(PMucros) 
PMucros <- select(PMucros, GeneName = Name, 
                     geneLength = Length, 
                     HMAJ_testis = TPM, 
                     HMAJ_liver = TPM.1, 
                     HMAJ_heart = TPM.2, 
                     ALA_vno = TPM.3, 
                     ALA_tailA2 = TPM.4, 
                     ATEN_tailA2 = TPM.5,
                     ATEN_tailB5 = TPM.6, 
                     ATEN_body = TPM.7, 
                     ALA_eye = TPM.8, 
                     BRH_vno = TPM.9, 
                     ALAjuv_body = TPM.10,
                     ALAjuv_tailA2 = TPM.11, 
                     ALAjuv_tailB5 = TPM.12, 
                     NSC_vno = TPM.13)


# Convert TPM to FPKM
# Load tpmToFpkm function
tpmToFpkm <- function(tpm, geneLength){
  count <- apply(tpm, 2, function(e){ e/1000*geneLength })
  fpkm <- apply(count, 2, function(e){ e*1000*1000000/sum(e)/geneLength })
  return(fpkm)
}

# Subset PMucros data object
geneLength <- PMucros$geneLength ## Doesn't matter which geneLength column you choose. They're all the SAME!!!!!

## Getting TPM dataframe
tpm <- PMucros %>%
  column_to_rownames("GeneName") %>%
  select(-contains("geneLength"))

## Generating FPKM values
Fpkm <- tpmToFpkm(tpm = tpm, geneLength = geneLength) %>% data.frame()

# Clustering analysis
dgePMucros <- select(PMucros, GeneName, HMAJ_testis, HMAJ_liver, HMAJ_heart, ALA_vno, ALA_tailA2, ATEN_tailA2,
                     ATEN_body, ALA_eye, BRH_vno, ALAjuv_body, ALAjuv_tailA2, ALAjuv_tailB5, NSC_vno)

## Turning counts object into DGElist and normalising
dgePMucros %<>% column_to_rownames("GeneName") %>%
  as.matrix() %>%
  edgeR::DGEList() %>%
  edgeR::calcNormFactors()

## Obtaining MDS data object 
mds <- limma::plotMDS(dgePMucros, gene.selection = "common", method = "logFC")

## Slightly ganky way of doing this, but extracting MDS coordinate object and joining with sample information.
plot <- mds@.Data[[3]] %>%
  as.data.frame() %>%
  rownames_to_column("group") %>%
  set_colnames(c("group","x","y")) %>%
  left_join(sampleinfo, by = c("group"="gene_id"))

## Plotting MDS using ggplot (nicer plot), includes overlap of experimental condition and species
pdf(file = "plots/MDS_PMucros.pdf")
plot %>%
  ggplot(aes(x, y, colour = experiment, shape = species, label = tissue)) +
  geom_point(size = 4) +
  geom_text(hjust = 0, nudge_x = 0.2, check_overlap = TRUE) + 
  theme_bw(base_size = 16) +
  theme(legend.text = element_text(size = 12)) + 
  labs(x = "Dimension 1",
       y=  "Dimension 2")
dev.off()

## Sample clustering
pdf(file = "plots/hclust_PMucros.pdf")
plot(hclust( d = dist(t(edgeR::cpm(dgePMucros)),  method = "euclidean"),
             method = "ward.D"),
     hang = -1, xlab = "Samples", sub = "Distance=euclidean; Method=ward.D")
dev.off()

## Heatmap
# Gene annotation
TPM_Pmucros <- read.delim("PMucros_quant_tpm_all.tsv")
geneNumToName <- read_delim("gene_number_to_name.tsv", delim = "\t", col_names = c("GeneID", "GeneName"))

TPM_Pmucros <- select(TPM_Pmucros, GeneName = Name, HMAJ_testis = TPM, HMAJ_heart = TPM.1, ALA_vno = TPM.2, ALA_tailA2 = TPM.3, ATEN_tailA2 = TPM.4,
                     ATEN_tailB5 = TPM.5, ATEN_body = TPM.6, ALA_eye = TPM.7, BRH_vno = TPM.8, ALAjuv_body = TPM.9, ALAjuv_tailA2 = TPM.10, ALAjuv_tailB5 = TPM.11, NSC_vno = TPM.12)

## Heatmap of TPM
# Turn into data matrix and log (R function log1p = log(1+x). used here because it accounts for those cells = 0)
TPM_Pmucros <- TPM_Pmucros %>%
  column_to_rownames("GeneName") %>%
  select(HMAJ_heart, ALA_vno, HMAJ_testis, ATEN_tailA2, BRH_vno, NSC_vno, ALA_tailA2, ALAjuv_body, ALAjuv_tailA2, ALAjuv_tailB5, 
         ATEN_body, ATEN_tailB5, ALA_eye) %>% # set order of the tissues accroding to hclust
  as.matrix() %>%
  log1p()

Fpkm_heatmap <- Fpkm %>%
  as.matrix() %>%
  log1p()

colfunc <- brewer.pal(n = 9, name = "OrRd")
pheatmap(TPM_Pmucros, cluster_rows = FALSE, cluster_cols = FALSE, color = colfunc, cutree_rows = 2 #or 3)

pheatmap(Fpkm_heatmap,
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         color = colfunc, 
         cutree_rows = 2) #or 3))


