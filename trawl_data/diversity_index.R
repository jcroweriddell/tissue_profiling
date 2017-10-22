library(dplyr)
library(magrittr)
library(tibble)
library(tidyverse)
library(ggplot2)
library(tidyr)

setwd("/Users/L033060262053/Documents/working_dir/activityPatterns/")

Dof <- read.csv("Dof_combined.csv", header = TRUE, sep = ",")# Diversity index between sites within EG

# Brillouin index (Hb) function
Hb <- function(ns) { 
  N <- sum(ns) 
  (lfactorial(N) - sum(lfactorial(ns)))/N 
} 

# Create data.frame with # of snakes caught per species per Site
ns <- Dof %>% 
  filter(Fishery == "EG") %>%
  group_by(Site) %>% 
  count(Species)

head(ns)

Hb(ns)


