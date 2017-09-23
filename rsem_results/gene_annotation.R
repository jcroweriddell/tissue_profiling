

library(tidyverse)

geneNumToName <- read_delim("gene_number_to_name.tsv", delim = "\t", col_names = c("GeneID", "GeneName"))

FPKMresults_gene <- FPKMresults %>%
  full_join(geneNumToName, by = c("gene_id" = "GeneID")) %>%
  select(GeneName, ALA_eye_FPKM, ALA_vno_FPKM, ALA_tailA2_FPKM, ALAjuv_tailA2_FPKM, 
         ALAjuv_tailB5_FPKM, ALAjuv_body_FPKM, ATEN_tailA2_FPKM, ATEN_tailB5_FPKM, 
         ATEN_body_FPKM, NSC_vno_FPKM, BRH_vno_FPKM, HMAJ_heart_FPKM, 
         HMAJ_testis_FPKM) 

ALA_A2tail<- FPKMresults_gene %>% arrange(desc(ALA_tailA2_FPKM))

View(ALA_A2tail)
head(FPKMResults_gene)

# I want to filter the counts by visual genes I'm interested in (create a vector then use match()?)

# Then create a heatmap of FPKM across tissues, ordered by visual genes of interest
