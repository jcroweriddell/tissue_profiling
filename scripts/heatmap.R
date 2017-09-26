####### Heatmap and clustering to show FPKM/TPM across samples for a subset of genes #########
library(edgeR)
library(ape)

## Read in data
FPKMresultsNum <- read.table("rsem_all_together.tsv")

## Reading in counts object 
dge <- read.table(file = "rsem_expectedCount_all_together_names.tsv", 
                  header = TRUE)
## Reading in sample object 
sample <- read_csv(file = "sample_info_tissues.csv", col_names = TRUE)

## Turning counts object into DGElist and normalising
dge %<>% column_to_rownames("gene_id") %>%
  as.matrix() %>%
  edgeR::DGEList() %>%
  calcNormFactors()

## Clustering of samples to check sample reproducibility 
d <- cor(dge, method = "spearman")
  

d <- cor(rpkmDFeByg, method="spearman")
hc <- hclust(dist(1-d))
pdf("./results/Sample_tree_rpkm.pdf") # save the tree as PDF file
plot.phylo(as.phylo(hc), type="p", edge.col=4, edge.width=3, show.node.label=TRUE, no.margin=TRUE)
dev.off()
library()