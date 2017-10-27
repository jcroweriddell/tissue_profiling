
library(readr)
library(dplyr)

df <- read_tsv(file = "/Users/L033060262053/Documents/Research projects/Tail_photoreception/tissue_profiling/salmon_quant/tpm_quant/Chicken_quant_tpm_all.tsv", col_names = TRUE)

fd <- select(df, GeneName = Name, 
                    HMAJ_testis = TPM,
                    HMAJ_liver = TPM_1,
                    HMAJ_heart = TPM_2,
                    ALA_vno = TPM_3,
                    ALA_tailA2 = TPM_4,
                    ATEN_tailA2 = TPM_6,
                    ATEN_tailB5 = TPM_7,
                    ATEN_body = TPM_8,
                    ALA_eye = TPM_5,
                    BRH_vno = TPM_9,
                    ALAjuv_body = TPM_12,
                    ALAjuv_tailA2 = TPM_10,
                    ALAjuv_tailB5 = TPM_11,
                    NSC_vno = TPM_13)
testis <- select(fd, GeneName, HMAJ_testis) %>% arrange(desc(HMAJ_testis)) %>%
  
write(x = testis$GeneName,file = "/Users/L033060262053/Desktop/test.txt")
