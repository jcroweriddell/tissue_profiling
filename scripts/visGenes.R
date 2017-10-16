library(magrittr)
library(tibble)
library(readr)
library(stringr)
library(dplyr)
library(tidyverse)

# Set working directory and read in data
setwd("/Users/jennacrowe-riddell/Documents/tissue_profiling/salmon_quant/tpm_results")
refGenome <- read.table(file = "PMucros_quant_tpm_all.tsv", header = TRUE)
sampleinfo <- read_csv(file = "sample_info_tissues.csv", col_names = TRUE)
VRgenes <- read_csv(file = "VRgenes_PMucros.tsv", col_names = FALSE)

# Change column names to tissue names
names(refGenome) 
refGenome <- select(refGenome, GeneName = Name, HMAJ_testis = TPM, HMAJ_heart = TPM.1, ALA_vno = TPM.2, ALA_tailA2 = TPM.3, ATEN_tailA2 = TPM.4,
                    ATEN_tailB5 = TPM.5, ATEN_body = TPM.6, ALA_eye = TPM.7, BRH_vno = TPM.8, ALAjuv_body = TPM.9, ALAjuv_tailA2 = TPM.10, ALAjuv_tailB5 = TPM.11, NSC_vno = TPM.12)


# Change column names, remove unwanted characters
VRgenes <- select(VRgenes, XM.geneid = X1, mRNA = X2)
VRgenes$XM.geneid <- gsub(">","",VRgenes$XM.geneid)


# Create new column for 'Predicted'
VRgenes <- separate(VRgenes, XM.geneid, c("XM.geneid", "Predicted"), sep = " PREDICTED: ")

# This also kinda works
# VRgenes<- VRgenes %>% mutate(Predicted = unlist(strsplit(XM.geneid, split = ".1_")))

# Save VR gene list
write.table(VRgenes, "VRgenes_XMdescription", sep = "\t")

# Create a vetor of XM values
geneList <- VRgenes$XM.geneid

# Filter refGenome by vno tissues then XM VR values
refGenomeVno <- select(refGenome, GeneName, ALA_vno, NSC_vno, BRH_vno)
VRGenesTPM <- filter(refGenomeVno, GeneName %in% geneList)

# Then look at TPM count table across tissues!
View(VRGenesTPM)
write.table(VRGenesTPM, "VRgenesTPM_PMucros", sep = "\t")

# Heatmap
VRGenesTPMlog <- VRGenesTPM %>%
  column_to_rownames("GeneName") %>%
  # select(ALA_tailA2, ALAjuv_tailB5, ALAjuv_tail, ALAjuv_body, ATEN_tail, ATEN_tailB5, ATEN_body,
  # HMAJ_heart, HMAJ_testis, ALA_eye, ALA_vno, BRH_vno, NSC_vno) %>%# set order of the tissues
  as.matrix() %>%
  log1p()

colfunc <- brewer.pal(n = 9, name = "OrRd")
pdf(file = "plots/VRgenes_heatmap_PMucros.pdf")
pheatmap(VRGenesTPMlog, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         color = colfunc,
         border_colour = NA)
dev.off()

quartz()