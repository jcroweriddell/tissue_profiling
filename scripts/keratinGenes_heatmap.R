## Generate heatmap for genes involved in Keratin

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
keratinGenes <- read_csv(file = "PMucros_genelist/PMucros_keratin.txt", col_names = FALSE)
head(keratinGenes)

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
keratinGenes$X1 <- gsub(">","", keratinGenes$X1)

# Create new column for 'Predicted' & 'geneName'
keratinGenes <- separate(keratinGenes, X1, c("XM.geneid", "Predicted"), sep = " PREDICTED: ")
keratinGenes$Predicted <- gsub("\\(", ",", keratinGenes$Predicted)

keratinGenes <- separate(keratinGenes, Predicted, c("Predicted", "geneName"), sep = ",")
keratinGenes$geneName <- gsub("\\)", "", keratinGenes$geneName) 
keratinGenes$geneName <- gsub("mRNA", "", keratinGenes$geneName) 
keratinGenes$geneName <- gsub("transcript variant*", "", keratinGenes$geneName) 
View(keratinGenes)

# Save VR gene list
write.table(keratinGenes, "PMucros_R_analysis/keratinGenes_XMdescription", sep = "\t")

# Create a vetor of XM values
geneList <- keratinGenes$XM.geneid

# Filter refGenome by XM vOPN genes
keratinGenesTPM <- filter(refGenome, GeneName %in% geneList)

# Then look at TPM count table across tissues!
View(keratinGenesTPM)
write.table(keratinGenesTPM, "PMucros_R_analysis/keratinGenesTPM", sep = "\t")

# Heatmap
keratinGenesTPMlog <- keratinGenesTPM %>%
  column_to_rownames("GeneName") %>%
  select(ALA_tailA2, ALAjuv_tailB5, ALAjuv_tailA2, ALAjuv_body, ATEN_tailA2, ATEN_tailB5, ATEN_body,
         HMAJ_heart, HMAJ_testis, ALA_eye, ALA_vno, BRH_vno, NSC_vno) %>%# set order of the tissues
  as.matrix() %>%
  log1p()

colfunc <- brewer.pal(n = 9, name = "OrRd")
pdf(file = "PMucros_R_analysis/keratinGenes_heatmap.pdf")
pheatmap(keratinGenesTPMlog, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         color = colfunc,
         breaks = c(0, 0.25, 0.5, 0.75, 1, 2, 3, 4.5, 6, 7),
         labels_row = keratinGenes$geneName,
         cutree_rows = 3)

dev.off()

