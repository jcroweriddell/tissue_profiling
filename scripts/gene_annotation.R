
library(tidyverse)
# read data
FPKMresults <- read.delim("rsem_all_together.tsv")
TPMresults <- read.delim("rsem_TPM_all_together.tsv")
geneNumToName <- read_delim("gene_number_to_name.tsv", delim = "\t", col_names = c("GeneID", "GeneName"))

#create table with FPKM and gene names
FPKMresults_gene <- FPKMresults %>%
  full_join(geneNumToName, by = c("gene_id" = "GeneID")) %>%
  select(GeneName, ALA_eye = ALA_eye_FPKM, ALA_vno = ALA_vno_FPKM , ALA_tailA2 = ALA_tailA2_FPKM, 
         ALAjuv_tail = ALAjuv_tailA2_FPKM, ALAjuv_tailB5 = ALAjuv_tailB5_FPKM,
         ALAjuv_body = ALAjuv_body_FPKM,ATEN_tail = ATEN_tailA2_FPKM, ATEN_tailB5= ATEN_tailB5_FPKM,
         ATEN_body =  ATEN_body_FPKM,NSC_vno = NSC_vno_FPKM, BRH_vno = BRH_vno_FPKM,
         HMAJ_heart =HMAJ_heart_FPKM, HMAJ_testis = HMAJ_testis_FPKM) 

# create TPM
TPMresults_gene <- TPMresults %>%
  full_join(geneNumToName, by = c("gene_id" = "GeneID")) %>%
  select(GeneName, ALA_eye = ALA_eye_TPM, ALA_vno = ALA_vno_TPM , ALA_tailA2 = ALA_tailA2_TPM, 
         ALAjuv_tail = ALAjuv_tailA2_TPM, ALAjuv_tailB5 = ALAjuv_tailB5_TPM,
         ALAjuv_body = ALAjuv_body_TPM,ATEN_tail = ATEN_tailA2_TPM, ATEN_tailB5= ATEN_tailB5_TPM,
         ATEN_body =  ATEN_body_TPM,NSC_vno = NSC_vno_TPM, BRH_vno = BRH_vno_TPM,
         HMAJ_heart =HMAJ_heart_TPM, HMAJ_testis = HMAJ_testis_TPM)
head(TPMresults_gene)

# filter by individual gene names
geneList <- c("GNAT1", "GNAT2", "GNAT3", "grk7-a", "GUCY2D",
"GUCY2F",  "CUCA1A", "LRAT", "PDE6B", "PDE6C", "RDH8", "RPE65",
"OPN3", "OPN4", "OPN5", "ROD1", "ABCA4", "ARRB2", "RHO", "L345_06817", "L345_18159", "L345_00724",
"L345_18159", "DHRS9")

visGenesFPKM <- filter(FPKMresults_gene, GeneName %in% geneList)
View(visGenesFPKM)

visGenesTPM <- filter(TPMresults_gene, GeneName %in% geneList)
View(visGenesTPM)

# save file
write.table(visGenesFPKM, "visGenesFPKM.txt", sep="\t")
write.table(visGenesTPM, "visGenesFPKM.txt", sep="\t")
# Heatmap of FPKM (or TPM) for visual genes

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







