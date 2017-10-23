library(dplyr)
library(magrittr)
library(tibble)
library(tidyverse)
library(ggplot2)
library(tidyr)

setwd("/Users/L033060262053/Documents/working_dir/activityPatterns/")

# site <- c(1,1,1,1,1,10,11,11,11,12,13,13,13,14,14)
# species <- LETTERS[1:15]
# n <- c(4,5,11,5,1,4,12,4,3,4,7,2,3,3,2)
# df <- data.frame(site,species,n)

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

## Al's solution
df %>% 
  split(f = .$site) %>%
  lapply(function(x){
    vec <- sum(x$n)
    (lfactorial(vec) - sum(lfactorial(x$n)))/vec
  }) %>%
  bind_rows()
  


