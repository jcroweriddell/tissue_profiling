###############################################################################
## Required packages
###############################################################################

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
library(BiocInstaller)
library(edgeR)
library(tidyr)
# library(limma)

###############################################################################
## User defined functions
###############################################################################

# Function - TPM to FPKM
TPMToFpkm <- function(TPM, geneLength){
  
  count <- apply(TPM, 2, function(e){ e/1000*geneLength })
  fpkm <- apply(count, 2, function(e){ e*1000*1000000/sum(e)/geneLength })
  
  return(fpkm)
  
}

###############################################################################
## Data directories + Data import
###############################################################################

# Set working directory (what's going on here?)
dir <- "~/Documents/tissue_profiling/salmon_quant/tpm_quant/"
vis_dir <- paste0(dir,"refGenomes/PMucros_genelist/")

## Importing data
refGenome <- read.table(file = paste0(dir,"PMucros_quant_tpm_all.tsv"), header = TRUE)
sampleinfo <- read_csv(file = paste0(dir,"sample_info_tissues.csv"), col_names = TRUE)
visGenes <- read_csv(file = paste0(vis_dir,"PMucros_visGenes5.txt"), col_names = FALSE)

###############################################################################
## Count matrix clean-up
###############################################################################

## Gene length - used in TPM --> FPKM conversions below
gene_length <- refGenome[["Length"]]

## Changing column name to tissue name
refGenome <- select(refGenome, 
                    GeneName = Name, 
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

## Selecting columns (potentially just leave these out in the step above?)
dge <- select(refGenome, -c(NSC_vno,
                            BRH_vno,
                            ALA_vno))

###############################################################################
## Generating DGElist object + MDS plots
###############################################################################

## Generating DGElist object
dge %<>% 
  column_to_rownames("GeneName") %>%
  as.matrix() %>%
  edgeR::DGEList() %>%
  edgeR::calcNormFactors()

## Obtaining MDS data object 
mds <- plotMDS(dge, 
               gene.selection = "common", 
               method = "logFC")

## Generating dataframe to plot
plotObj <- mds@.Data[[3]] %>%
  as.data.frame() %>%
  rownames_to_column("group") %>%
  rename(x = V1, y = V2) %>%
  left_join(sampleinfo, by = c("group"="gene_id"))

## MDS plot
pdf(file = ..., width = , height = , ...) ## Prefer this to save files
plotgg <- plotObj %>%
  ggplot(aes(x, y, colour = Species, shape = Species, label = Tissue)) +
  geom_point(size = 4) +
  geom_text(hjust = 0, nudge_x = 0.1, check_overlap = TRUE) + 
  theme_bw(base_size = 15) +
  theme(legend.text = element_text(size = 12)) + 
  labs(x = "Dimension 1",
       y=  "Dimension 2") + 
  geom_vline(xintercept = 0, linetype = 3) + 
  geom_hline(yintercept = 0, linetype = 3) +
  theme_linedraw(base_size = 20)
dev.off()

###############################################################################
## Converting TPM to FPKM
###############################################################################

## Getting TPM dataframe - This hasn't been updated (geneLength isn't in any of the above regGenome objects...)
TPM <- refGenome %>%
  column_to_rownames("GeneName") #%>%
  #select(-contains("geneLength"))

## Generating FPKM values
Fpkm <- TPMToFpkm(TPM = TPM, geneLength = gene_length) %>% 
  data.frame() %>% 
  rownames_to_column("GeneName")

###############################################################################
## Generating visual genes table
###############################################################################

## Clean up visGenes table
visGenes %<>%
  mutate(X1=gsub(">","",X1)) %>%
  separate(X1, c("XM.geneid", "Predicted"), sep = " PREDICTED: ") %>%
  mutate(Predicted=gsub("\\(", ",", Predicted)) %>%
  separate(Predicted, c("Predicted", "geneName"), sep = ",") %>%
  mutate(geneName=gsub("\\)", "",geneName),
         geneName=gsub("mRNA","",geneName),
         geneName=gsub("transcript variant*","",geneName))

# Save VR gene list
write.table(visGenes, "refGenomes/PMucros_heatmaps/visGenes_XMdescription", sep = "\t")

# Create a vetor of XM values
geneList <- visGenes[["XM.geneid"]]

###############################################################################
## Filter refGenome by RefSeq (XM) gene IDs
###############################################################################

# Filter refGenome by RefSeq (XM) gene IDs
visGenesFpkm <- filter(refGenome, GeneName %in% geneList) %>%
  left_join(y = visGenes, by = c("GeneName" = "XM.geneid")) %>%
  select(-GeneName, -Predicted) %>%
  arrange(geneName)

# Then look at TPM count table across tissues!
View(visGenesFpkm)
write.table(visGenesFpkm, "refGenomes/PMucros_heatmaps/visGenesFpkm", sep = "\t")

###############################################################################
# Heatmap visualisation
###############################################################################

visGenesFpkmlog <- visGenesFpkm %>%
  column_to_rownames("geneName") %>%
  select(ALA_eye, ALA_tailA2,ALAjuv_tailB5, ALAjuv_tailA2, ALAjuv_body, 
         ATEN_tailA2, ATEN_tailB5, ATEN_body,HMAJ_testis, HMAJ_heart, HMAJ_liver) %>%# set order of the tissues
  as.matrix() %>%
  log1p()

colfunc <- brewer.pal(n = 9, name = "OrRd")

pheatmap(visGenesFpkmlog, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         color = colfunc,
         border_color = NA,
         breaks = c(0, 0.25, 0.5, 0.75, 1, 2, 4, 5.5, 7, 9),
         filename = "Pmucros_visGenes_Fpkmheatmap.jpeg")

dev.off()

sessionInfo()

