

library(tidyverse)
# read data
FPKMresults
TPMresults <- read.delim("rsem_TPM_all_together.tsv")
geneNumToName <- read_delim("gene_number_to_name.tsv", delim = "\t", col_names = c("GeneID", "GeneName"))

#create FPKM
FPKMresults_gene <- FPKMresults %>%
  full_join(geneNumToName, by = c("gene_id" = "GeneID")) %>%
  select(GeneName, ALA_eye_FPKM, ALA_vno_FPKM, ALA_tailA2_FPKM, ALAjuv_tailA2_FPKM, 
         ALAjuv_tailB5_FPKM, ALAjuv_body_FPKM, ATEN_tailA2_FPKM, ATEN_tailB5_FPKM, 
         ATEN_body_FPKM, NSC_vno_FPKM, BRH_vno_FPKM, HMAJ_heart_FPKM, 
         HMAJ_testis_FPKM) 
# create TPM
TPMresults_gene <- TPMresults %>%
  full_join(geneNumToName, by = c("gene_id" = "GeneID")) %>%
  select(GeneName, ALA_eye_TPM, ALA_vno_TPM, ALA_tailA2_TPM, ALAjuv_tailA2_TPM, 
         ALAjuv_tailB5_TPM, ALAjuv_body_TPM, ATEN_tailA2_TPM, ATEN_tailB5_TPM, 
         ATEN_body_TPM, NSC_vno_TPM, BRH_vno_TPM, HMAJ_heart_TPM, 
         HMAJ_testis_TPM) 

head(TPMresults_gene)
#ALA_A2tail<- FPKMresults_gene %>% arrange(desc(ALA_tailA2_FPKM))

View(ALA_A2tail)
head(FPKMResults_gene)

# I want to filter the counts by visual genes I'm interested in (create a vector then use match()?)


# filter by individual gene names
geneList <- c("GNAT1", "GNAT2", "GNAT3", "grk7-a", "GUCY2D",
"GUCY2F",  "CUCA1A", "LRAT", "PDE6B", "PDE6C", "RDH8", "RPE65",
"OPN3", "OPN4", "OPN5", "ROD1", "ABCA4", "ARRB2", "RHO", "L345_06817", "gar1", "L345_18159", "L345_00724",
"L345_18159", "RDH16", "DHRS9", "RDH8")

geneList

visGenesFPKM <- filter(FPKMresults_gene, GeneName %in% geneList)
View(visGenesFPKM)
visGenesTPM <- filter(TPMresults_gene, GeneName %in% geneList)
View(visGenesTPM)

# save file
write.table(visGenesFPKM, "visGenesFPKM.txt", sep="\t")

# Then create a heatmap of FPKM (or TPM) across tissues, ordered by visual genes of interest
