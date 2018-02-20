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

# Al - I set these paths for my computer. Changes these for your own

dir <- "../salmon_quant/tpm_quant/"
vis_dir <- "~/Projects_BioinformaticsHub/Kate_Sanders/tissue_profiling/pheatmap_jenna/"
write_dir <- "path/to/save/data/objects"

## NOTE 1: You shouldn't have to set `vis_dir` as it will be the working
## directory of this project. I set it here just so you could see 
## what I was doing.

## NOTE 2: The path to `write_dir` will be a relative path, which would look
## something like this `../path/to/output`. 

## Importing data
refGenome <- read.table(file = paste0(dir,"PMucros_quant_tpm_all.tsv"), header = TRUE)
sampleinfo <- read_csv(file = paste0(dir,"sample_info_tissues.csv"), col_names = TRUE)

## Importing visual genes of
fl <- list.files(path = vis_dir, pattern = "PMucros_", full.names = TRUE)
names(fl) <- c("opn","retinal","visGenes5")
vis_lst <- lapply(fl, function(x){
  read_csv(file = x, col_names = FALSE)
})

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

# ## MDS plot
# pdf(file = ..., width = , height = , ...) ## Prefer this to save files
plotObj %>%
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
# dev.off()

###############################################################################
## Converting TPM to FPKM
###############################################################################

## Getting TPM dataframe (this is literally just turning a column to rownames)
TPM <- refGenome %>%
  column_to_rownames("GeneName")

## Generating FPKM values
Fpkm <- TPMToFpkm(TPM = TPM, geneLength = gene_length) %>% 
  data.frame() %>% 
  rownames_to_column("GeneName")

###############################################################################
## Generating visual genes table
###############################################################################

visGenes <- lapply(names(vis_lst), function(x){
  df <- vis_lst[[x]] %>%
    mutate(X1=gsub(">","",X1)) %>%
    separate(X1, c("XM.geneid", "Predicted"), sep = " PREDICTED: ") %>%
    mutate(Predicted=gsub("\\(", ",", Predicted)) %>%
    separate(Predicted, c("Predicted", "geneName"), sep = ",") %>%
    mutate(geneName=gsub("\\)", "",geneName),
           geneName=gsub("mRNA","",geneName),
           geneName=gsub("transcript variant*","",geneName)) %>%
    mutate(group = x)
}) %>% bind_rows() %>%
  arrange(group)

## Setting the groups column as a factor variable - helps with orienting the data 
visGenes$group <- factor(visGenes$group)

# Create a vetor of XM values
geneList <- visGenes[["XM.geneid"]]

###############################################################################
## Filter refGenome by RefSeq (XM) gene IDs
###############################################################################

# Filter refGenome by RefSeq (XM) gene IDs
visGenesFpkm <- filter(refGenome, GeneName %in% geneList) %>%
  left_join(y = visGenes, by = c("GeneName" = "XM.geneid")) %>%
  select(-GeneName, -Predicted) %>%
  arrange(group)

###############################################################################
# Heatmap visualisation
###############################################################################

## Generating matrix object
visGenesFpkmlog <- visGenesFpkm %>%
  column_to_rownames("geneName") %>%
  select(-group) %>%
  select(ALA_eye, ALA_tailA2,ALAjuv_tailB5, ALAjuv_tailA2, ALAjuv_body,
         ATEN_tailA2, ATEN_tailB5, ATEN_body,HMAJ_testis, HMAJ_heart, HMAJ_liver) %>%# set order of the tissues
  as.matrix() %>%
  log1p()

## Setting the annotation dataframe (this is used to generate the category legend)
ann_df <- select(visGenesFpkm, geneName, Category = group) %>% 
  column_to_rownames("geneName")

## Colour for the heatmap
colfunc <- brewer.pal(n = 9, name = "OrRd")

## Colours for the categories
mycolors <- brewer.pal(8, "Accent")[c(1,4,5)]
names(mycolors) <- unique(ann_df$Category)
anno_colours <- list(Category = mycolors)

## plotting function
# pdf(filename = paste0(write_dir,"/pheatmap_visualGenes.png"),width = ...,height = ...)
pheatmap(visGenesFpkmlog, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         color = colfunc,
         border_color = NA,
         gaps_row = c(6,19),
         annotation_row = ann_df,
         annotation_colors = anno_colours,
         breaks = c(0, 0.25, 0.5, 0.75, 1, 2, 4, 5.5, 7, 9))
# dev.off()

###############################################################################
## Writing useful data objects
###############################################################################

saveRDS(object = dge, file = paste0(write_dir,"/dge.RDS")) 
saveRDS(object = Fpkm, file = paste0(write_dir,"/Fpkm.RDS")) 
saveRDS(object = visGenes, file = paste0(write_dir,"/visGenes.RDS")) 
saveRDS(object = visGenesFpkmlog, file = paste0(write_dir,"/visGenesFpkmlog.RDS"))  visGenesFpkmlog

###############################################################################

dev.off()

sessionInfo()

