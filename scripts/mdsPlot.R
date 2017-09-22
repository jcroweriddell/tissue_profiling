library(magrittr)
library(tibble)
library(readr)
library(dplyr)
library(ggplot2)
library(edgeR)

## Reading in counts object 
dge <- read.table(file = "~/Desktop/tissue_profiling-master/rsem_results/rsem_expectedCount_all_together_names.tsv", 
                  header = TRUE)
## Reading in sample object 
sample <- read_csv(file = "~/Desktop/tissue_profiling-master/rsem_results/sample_info_tissues.csv", col_names = TRUE)

## Turning counts object into DGElist and normalising
dge %<>% column_to_rownames("gene_id") %>%
  as.matrix() %>%
  edgeR::DGEList() %>%
  calcNormFactors()

## Obtaining MDS data object 
mds <- plotMDS(dge)

## Slightly ganky way of doing this, but extracting MDS coordinate object and joining with sample information.
plot <- mds@.Data[[3]] %>%
  as.data.frame() %>%
  rownames_to_column("group") %>%
  set_colnames(c("group","x","y")) %>%
  left_join(sample, by = c("group"="gene_id"))

## Plotting MDS using ggplot (nicer plot)
plot %>%
  ggplot(aes(x, y, colour = species, label = group)) +
  geom_point(size = 4) +
  geom_text(hjust = 0, nudge_x = 0.2, check_overlap = TRUE) + 
  theme_bw(base_size = 16) +
  theme(legend.text = element_text(size = 12)) + 
  labs(x = "Dimension 1",
       y=  "Dimension 2")
 



