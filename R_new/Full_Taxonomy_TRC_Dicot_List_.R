## Script to create TRC Dicot list for Beth

rm(list = ls())

# devtools::install_github('levisc8/Fun_Phylo_Package')
library(FunPhylo)
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(taxize)

data(tyson)

SppList <- tyson$spp.list

# Create names list
SppNames <- gsub('_', ' ', unique(SppList$Species))

# Get full classification if it does exist
FullTaxonomy <- classification(SppNames, db = 'itis', rows = 1)

# remove NAs given by someGenus_spp
FullTaxonomy <- FullTaxonomy[!is.na(FullTaxonomy)]

# Reshape the data bunch to put it into a single data frame
 
FullTaxaTable <- ldply(FullTaxonomy,
                       .fun = dcast,
                       formula = rank ~ .,
                       value.var = 'name',
                       .progress = progress_text()) %>%
  setNames(c('Species', 'rank', 'value')) %>%
  spread(key = rank, value = value) %>%
  select(-species) %>%
  mutate(Species = gsub(' ', '_', .$Species)) %>%
  right_join(., SppList, by = 'Species') %>%
  .[!duplicated(.$Species), ]


# save output!

write.csv(FullTaxaTable, file = '../Data/TRC_Dicot_List.csv', 
          row.names = FALSE, quote = FALSE)
  


