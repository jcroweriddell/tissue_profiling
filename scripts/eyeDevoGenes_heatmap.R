## Generate heatmap for genes involved in Eye development genes

library(magrittr)
library(tibble)
library(readr)
library(stringr)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(tidyselect)

# Set working directory and read in data
setwd("/Users/L033060262053/Documents/Research projects/Tail_photoreception/tissue_profiling/salmon_quant/tpm_quant/")
refGenome <- read.table(file = "PMucros_quant_tpm_all.tsv", header = TRUE)
sampleinfo <- read_csv(file = "sample_info_tissues.csv", col_names = TRUE)
eyeDevoGenes <- read_csv(file = "PMucros_genelist/PMucros_eyeDevo.txt", col_names = FALSE)
head(eyeDevoGenes)

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

# Change column names, remove unwanted characters
eyeDevoGenes$X1 <- gsub(">","", eyeDevoGenes$X1)

# Create new column for 'Predicted' & 'geneName'
eyeDevoGenes <- separate(eyeDevoGenes, X1, c("GeneName", "Predicted"), sep = " PREDICTED: ")
eyeDevoGenes$Predicted <- gsub("\\(", ",", eyeDevoGenes$Predicted)

eyeDevoGenes <- separate(eyeDevoGenes, Predicted, c("Predicted", "geneDescript"), sep = ",")
eyeDevoGenes$geneDescript <- gsub("\\)", "", eyeDevoGenes$geneDescript) 
eyeDevoGenes$geneDescript <- gsub("mRNA", "", eyeDevoGenes$geneDescript) 
eyeDevoGenes$geneDescript <- gsub("transcript variant*", "", eyeDevoGenes$geneDescript) 
View(eyeDevoGenes)

# Save VR gene list
write.table(eyeDevoGenes, "PMucros_R_analysis/eyeDevoGenes_XMdescription", sep = "\t")

# Create a vetor of XM values
geneList <- eyeDevoGenes$XM.geneid

# Filter refGenome by XM vOPN genes
eyeDevoGenesTPM <- filter(refGenome, GeneName %in% geneList)

# Then look at TPM count table across tissues!
View(eyeDevoGenesTPM)
write.table(eyeDevoGenesTPM, "PMucros_R_analysis/eyeDevoGenesTPM", sep = "\t")

# Heatmap
eyeDevoGenesTPMlog <- eyeDevoGenesTPM %>%
  column_to_rownames("GeneName") %>%
  select(ALA_tailA2, ALAjuv_tailB5, ALAjuv_tailA2, ALAjuv_body, ATEN_tailA2, ATEN_tailB5, ATEN_body,
         HMAJ_heart, HMAJ_liver, HMAJ_testis, ALA_eye, ALA_vno, BRH_vno, NSC_vno) %>%# set order of the tissues
  as.matrix() %>%
  log1p()

colfunc <- brewer.pal(n = 9, name = "OrRd")
pdf(file = "PMucros_R_analysis/eyeDevoGenes_heatmap.pdf")
pheatmap(eyeDevoGenesTPMlog, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         color = colfunc,
         cutree_rows = 4,
         breaks = c(0, 0.5, 1, 2, 4, 8, 3, 6, 12, 20))

dev.off()

