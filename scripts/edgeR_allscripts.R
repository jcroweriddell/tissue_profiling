library(magrittr)
library(tibble)
library(readr)
library(dplyr)
library(ggplot2)
library(edgeR)

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

## Obtaining MDS data object 
mds <- plotMDS(dge, gene.selection = "common", method = "logFC")

## Slightly ganky way of doing this, but extracting MDS coordinate object and joining with sample information.
plot <- mds@.Data[[3]] %>%
  as.data.frame() %>%
  rownames_to_column("group") %>%
  set_colnames(c("group","x","y")) %>%
  left_join(sample, by = c("group"="gene_id"))

## Plotting MDS using ggplot (nicer plot), includes overlap of experimental condition and species
plot %>%
  ggplot(aes(x, y, colour = experiment, shape = species, label = tissue)) +
  geom_point(size = 4) +
  geom_text(hjust = 0, nudge_x = 0.2, check_overlap = TRUE) + 
  theme_bw(base_size = 16) +
  theme(legend.text = element_text(size = 12)) + 
  labs(x = "Dimension 1",
       y=  "Dimension 2")

## Same plot but with sequence data and species overlayed

plot %>%
  ggplot(aes(x, y, colour = seqYear, shape = species, label = tissue)) +
  geom_point(size = 4) +
  geom_text(hjust = 0, nudge_x = 0.2, check_overlap = TRUE) + 
  theme_bw(base_size = 16) +
  theme(legend.text = element_text(size = 12)) + 
  labs(x = "Dimension 1",
       y=  "Dimension 2")

## Sample clustering

plot(hclust( d = dist(t(cpm(dge)),  method = "euclidean"),
             method = "ward.D"),
     hang = -1, xlab = "Samples", sub = "Distance=euclidean; Method=ward.D")


## DE gene analyses
# Determine expression level filter
names(dge)
read_counts_sample1 <- dge$counts[ ,"ALA_eye_count" ]
cpm_sample1 <- cpm(read_counts_sample1)

for(count in c(5:10)){
  # find the CPM value of genes which have read counts between 5 and 10
  unique_cpm_value <- unique( cpm_sample1[ which( read_counts_sample1 == count ), ] )
  # print the CPM value rounded to 2 decimal places that corresponds to each count value between 5 and 10
  print( paste( count, "reads => CPM", round( unique_cpm_value, 2) ) )
}

# Filter out very lowly expressed genes, keeping genes that are expressed at a reasonable level (CPM>1) in at least 4 of the 8 samples
dge_expressed <- dge[ rowSums(cpm(dge) > 1) >= 4, ]
dim(dge_expressed$counts)

# Re-compute the library sizes after filtering
dge_expressed$samples$lib.size <- colSums(dge_expressed$counts)

# Visualise the library sizes for each sample
barplot(dge_expressed$samples$lib.size*1e-6, names=row.names( dge_expressed$samples ),
        xlab="Samples", ylab="Library size (millions)", cex.names = 0.6)
# Normalizing - Compute effective library sizes using the default TMM normalization:
dge_expressed <- calcNormFactors(dge_expressed, method="TMM")
dge_expressed$samples

### Defining the design matrix
# convert photoresponse information to factors
Photo <- as.factor(sample$experiment)
Photo

# create design matrix
design <- model.matrix(~ Photo)
design

# estimate a common negative binomial dispersion parameter for a DGE dataset with 
# a general experimental design
dispersions <- estimateGLMCommonDisp(dge_expressed, design, verbose=TRUE)

# estimate the abundance-dispersion trend by Cox-Reid approximate profile likelihood
dispersions <- estimateGLMTrendedDisp(dispersions, design)

# implement the empirical Bayes strategy proposed by McCarthy et al (2012) for estimating the tagwise negative binomial dispersions
dispersions <- estimateGLMTagwiseDisp(dispersions, design)

# plot biological coefficient of variation, indicates biological variation between replicate samples
plotBCV(dispersions, main="BCV plot")

# fits genewise negative binomial GLMs for design
fit <- glmFit(dispersions, design)
### BVC >0.6 too high

## Heatmap
log2cpm_expressed <- cpm(dge_expressed, log=T)
heatmap(log2cpm_expressed, col= heat.colors(299))
