## Generate heatmap for genes involved in Retinol metabolism

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
refGenome <- read.table(file = "Chicken_quant_tpm_all.tsv", header = TRUE)
retGenes <- read_csv(file = "refGenomes/Gallus_genelist/Gallus_retGenes.txt", col_names = FALSE)

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
# Create new column for 'Predicted' & 'geneName'
retGenes$X1 <- gsub(">","", retGenes$X1)
retGenes <- separate(retGenes, X1, c("XM.geneid", "Predicted"), sep = " PREDICTED: ")
retGenes$Predicted <- gsub("\\(", ",", retGenes$Predicted)
retGenes <- separate(retGenes, Predicted, c("Predicted", "geneName"), sep = ",")
retGenes$geneName <- gsub("\\)", "", retGenes$geneName) 
retGenes$geneName <- gsub("mRNA", "", retGenes$geneName) 
retGenes$geneName <- gsub("transcript variant*", "", retGenes$geneName) 
View(retGenes)

# Save VR gene list
write.table(retGenes, "refGenomes/Gallus_heatmaps/retGenes_XMdescription", sep = "\t")

# Create a vetor of XM values
geneList <- retGenes$XM.geneid

# Filter refGenome by XM vOPN genes
retGenesTPM <- filter(refGenome, GeneName %in% geneList) %>% 
  left_join(y = retGenes, by = c("GeneName" = "XM.geneid")) %>% 
  select(-GeneName, -Predicted) %>% 
  arrange(geneName)

# Then look at TPM count table across tissues!
View(retGenesTPM)
write.table(retGenesTPM, "refGenomes/Gallus_heatmaps/retGenesTPM", sep = "\t")

# Heatmap
retGenesTPMlog <- retGenesTPM %>%
  column_to_rownames("geneName") %>%
  select(ALA_eye, ALA_tailA2,ALAjuv_tailB5, ALAjuv_tailA2, ALAjuv_body, ATEN_tailA2, ATEN_tailB5, ATEN_body,
         HMAJ_testis, HMAJ_heart, HMAJ_liver) %>%# set order of the tissues
  as.matrix() %>%
  log1p()

colfunc <- brewer.pal(n = 9, name = "OrRd")
#pdf(file = "PMucros_R_analysis/retGenes_heatmap.pdf")
quartz()
pheatmap(retGenesTPMlog, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         color = colfunc,
         border_color = NA,
         breaks = c(0, 0.25, 0.5, 0.75, 1, 2, 3, 4.5, 6, 7),
         main = "Gallus- retinal metabolism genes",
         filename = "Gallus_retGenes_heatmap.jpeg")
         #labels_row = retGenes$geneName)

#dev.off()

