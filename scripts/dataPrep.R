## Data pre-processing using CPM (counts per million) in edgeR
# tutorial http://www.bioconductor.org/help/workflows/RNAseq123/#normalising-gene-expression-distributions

library(edgeR)
library(RColorBrewer)

## Reading in counts object 
dge <- read.table(file = "rsem_expectedCount_all_together_names.tsv", 
                  header = TRUE)
## Reading in sample object 
sample <- read_csv(file = "sample_info_tissues.csv", col_names = TRUE)

dge %<>% column_to_rownames("gene_id") %>%
  as.matrix() %>%
  edgeR::DGEList() #%>%
 #calcNormFactors()
x <- dge
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

# Removing lowly expressed genes 
# Check how many genes counts = 0 for all samples, set cut-off cpm >`1`in >3 samples
table(rowSums(x$counts==0)==9)
keep.exprs <- rowSums(cpm>1)>=3
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

# Compare pre/post filtering of lowly expressed gened
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
#dev.off()

# normalise by method of trimmed mean of M-values (TMM)
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors



