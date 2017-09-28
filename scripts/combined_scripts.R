library(magrittr)
library(tibble)
library(readr)
library(dplyr)
library(ggplot2)
library(edgeR)
library(RColorBrewer)
library(tidyverse)

## Clustering data using MDS plot and dendrograms

## Set working directory
setwd("C:/Users/L033060262053/Documents/Research projects/Tail_photoreception/tissue_profiling/rsem_results/")

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

## Visualise RNA read data

x <- dge
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

## Removing lowly expressed genes 
## Check how many genes counts = 0 for all samples, set cut-off cpm >`1`in >3 samples

table(rowSums(x$counts==0)==9)
keep.exprs <- rowSums(cpm>1)>=3
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

## Compare pre/post filtering of lowly expressed gened
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", sample$gene_id, text.col=col, text.width=5, bty="n")

lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", sample$gene_id, text.col=col, text.width=5, bty="n")
# dev.off()

## Normalise by method of trimmed mean of M-values (TMM)
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

#####################################################################################

## DE gene analyses
## Determine expression level filter
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

#####################################################################################

## Gene annotation

## Heatmap
## Read data

FPKMresults <- read.delim("rsem_all_together.tsv")
TPMresults <- read.delim("rsem_TPM_all_together.tsv")
geneNumToName <- read_delim("gene_number_to_name.tsv", delim = "\t", col_names = c("GeneID", "GeneName"))

## Create table with FPKM and gene names

FPKMresults_gene <- FPKMresults %>%
  full_join(geneNumToName, by = c("gene_id" = "GeneID")) %>%
  select(GeneName, ALA_eye = ALA_eye_FPKM, ALA_vno = ALA_vno_FPKM , ALA_tailA2 = ALA_tailA2_FPKM, 
         ALAjuv_tail = ALAjuv_tailA2_FPKM, ALAjuv_tailB5 = ALAjuv_tailB5_FPKM,
         ALAjuv_body = ALAjuv_body_FPKM,ATEN_tail = ATEN_tailA2_FPKM, ATEN_tailB5= ATEN_tailB5_FPKM,
         ATEN_body =  ATEN_body_FPKM,NSC_vno = NSC_vno_FPKM, BRH_vno = BRH_vno_FPKM,
         HMAJ_heart =HMAJ_heart_FPKM, HMAJ_testis = HMAJ_testis_FPKM) 

## Create TPM

TPMresults_gene <- TPMresults %>%
  full_join(geneNumToName, by = c("gene_id" = "GeneID")) %>%
  select(GeneName, ALA_eye = ALA_eye_TPM, ALA_vno = ALA_vno_TPM , ALA_tailA2 = ALA_tailA2_TPM, 
         ALAjuv_tail = ALAjuv_tailA2_TPM, ALAjuv_tailB5 = ALAjuv_tailB5_TPM,
         ALAjuv_body = ALAjuv_body_TPM,ATEN_tail = ATEN_tailA2_TPM, ATEN_tailB5= ATEN_tailB5_TPM,
         ATEN_body =  ATEN_body_TPM,NSC_vno = NSC_vno_TPM, BRH_vno = BRH_vno_TPM,
         HMAJ_heart =HMAJ_heart_TPM, HMAJ_testis = HMAJ_testis_TPM)
head(TPMresults_gene)

## Filter by individual gene names

geneList <- c("GNAT1", "GNAT2", "GNAT3", "grk7-a", "GUCY2D",
              "GUCY2F",  "CUCA1A", "LRAT", "PDE6B", "PDE6C", "RDH8", "RPE65",
              "OPN3", "OPN4", "OPN5", "ROD1", "ABCA4", "ARRB2", "RHO", "L345_06817", "L345_18159", "L345_00724",
              "L345_18159", "DHRS9")

visGenesFPKM <- filter(FPKMresults_gene, GeneName %in% geneList)
View(visGenesFPKM)

visGenesTPM <- filter(TPMresults_gene, GeneName %in% geneList)
View(visGenesTPM)

## Save file

write.table(visGenesFPKM, "visGenesFPKM.txt", sep="\t")
write.table(visGenesTPM, "visGenesFPKM.txt", sep="\t")

#####################################################################################
## Heatmap of FPKM (or TPM) for visual genes
# Turn into data matrix and log (R function log1p = log(1+x). used here because it accounts for those cells = 0)

visGeneslogFPKM <- visGenesFPKM %>%
  column_to_rownames("GeneName") %>%
  as.matrix() %>%
  log1p()

colfunc <-colorRampPalette(bias = 1, c("royalblue", "springgreen","yellow","red"))
heatmap(visGeneslogFPKM, col= colfunc(50))

visGeneslogTPM <- visGenesTPM %>%
  column_to_rownames("GeneName") %>%
  as.matrix() %>%
  log1p()

heatmap(visGeneslogTPM, col= colfunc(50))

dev.off()



