## Generate heatmap for genes involved in Phototransduction
# fpkm

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
visGenes <- read_csv(file = "PMucros_genelist/PMucros_visGenes2.txt", col_names = FALSE)

# Change column names to tissue names
names(refGenome) 
refGenome <- select(refGenome, Length, GeneName = Name, 
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
# Load tpmToFpkm function
tpmToFpkm <- function(tpm, geneLength){
  count <- apply(tpm, 2, function(e){ e/1000*geneLength })
  fpkm <- apply(count, 2, function(e){ e*1000*1000000/sum(e)/geneLength })
  return(fpkm)
}

# Subset PMucros data object
geneLength <- refGenome$Length ## Doesn't matter which geneLength column you choose. They're all the SAME!!!!!

## Getting TPM dataframe
tpm <- refGenome %>%
  column_to_rownames("GeneName") %>%
  select(-contains("geneLength"))

## Generating FPKM values
FpkmRefGenome <- tpmToFpkm(tpm = tpm, geneLength = geneLength) %>% 
  data.frame() %>% 
  rownames_to_column("GeneName")

## Clean up visGenes table
# Change column names, remove unwanted characters
# Create new column for 'Predicted' & 'geneName'
visGenes$X1 <- gsub(">","", visGenes$X1)
visGenes <- separate(visGenes, X1, c("XM.geneid", "Predicted"), sep = " PREDICTED: ")
visGenes$Predicted <- gsub("\\(", ",", visGenes$Predicted)
visGenes <- separate(visGenes, Predicted, c("Predicted", "geneName"), sep = ",")
visGenes$geneName <- gsub("\\)", "", visGenes$geneName) 
visGenes$geneName <- gsub("mRNA", "", visGenes$geneName) 

View(visGenes)

# Save VR gene list
#write.table(visGenes, "PMucros_R_analysis/visGenes_XMdescription", sep = "\t")

# Create a vetor of XM values
geneList <- visGenes$XM.geneid

# Filter refGenome by XM vOPN genes
visGenesFPKM <- filter(FpkmRefGenome, GeneName %in% geneList) %>%
  left_join(y = visGenes, by = c("GeneName" = "XM.geneid")) %>%
  select(-GeneName, -Predicted) %>%
  arrange(geneName)

# Then look at TPM count table across tissues!
View(visGenesFPKM)
write.table(visGenesTPM, "PMucros_R_analysis/visGenesTPM", sep = "\t")

# Heatmap
visGenesFPKMlog <- visGenesFPKM %>%
  column_to_rownames("geneName") %>%
  select(ALA_eye, ALA_tailA2,ALAjuv_tailB5, ALAjuv_tailA2, ALAjuv_body, ATEN_tailA2, ATEN_tailB5, ATEN_body,
         HMAJ_testis, HMAJ_heart, HMAJ_liver) %>%# set order of the tissues
  as.matrix() %>%
  log1p()

colfunc <- brewer.pal(n = 9, name = "OrRd")

#pdf(file = "plots/visGenes_heatmap.pdf")
pheatmap(visGenesFPKMlog, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         color = colfunc,
         border_color = NA,
         breaks = c(0, 0.25, 0.5, 0.75, 1, 2, 4, 5.5, 7, 9),
         main = "Phototransduction genes")
         #filename = "/Users/L033060262053/Dropbox/Transcriptome_results/visual_heatmap.jpeg")
        # labels_row = visGenes$geneName)

#dev.off()

