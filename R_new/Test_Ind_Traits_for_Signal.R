# Compute phylogenetic signal for individual and whole suite of traits
# I'll use Blomberg's K for continuous traits, D-test for binary traits 
# (Fritz + Purvis 2010), and a 
# phylogenetic Mantel Test for the whole suite (Hardy + Pavoine 2012). 

rm(list = ls(all = T))

library(FunPhylo)
library(phytools)
library(caper)
library(vegan)
library(ade4)
library(dplyr)

TraitData <- tyson$traits %>% filter(Species.Name != 'Cerastium_spp.')
phylo <- tyson$phylo

TraitPhy <- drop.tip(phylo, setdiff(phylo$tip.label, TraitData$Species.Name))

GrowthForm <- c(names(TraitData)[10:16])
Legume <- 'N_Fixer'

DispersalTraitNames<-c(names(TraitData)[17:25])
# continuous traits for testing: SLA, Height, Toughness, Wood Density

SLA <- TraitData$SLA
names(SLA) <- TraitData$Species.Name
Height <- TraitData$Height
names(Height) <- TraitData$Species.Name
Tough <- TraitData$Tough
names(Tough) <- TraitData$Species.Name
WoodDens <- TraitData$WoodDens
names(WoodDens) <- TraitData$Species.Name

outData <- data.frame(Trait = c('SLA','Height','Leaf toughness','Stem density',
                                'Flower Period',
                                GrowthForm, Legume, DispersalTraitNames, 'All combined'),
                      DataFormat = c(rep('Continuous',4), 'Circular',
                                     rep('Binary', 
                                         length(DispersalTraitNames) +
                                           length(GrowthForm) + 1),
                                     'Gower Distance Matrix'),
                      Test = c(rep("Blomberg's K", 4),
                               'Mantel Test',
                               rep('Phylogenetic D',
                                   length(DispersalTraitNames) +
                                     length(GrowthForm) + 1),
                               'Mantel Test'),
                      N = NA,
                      Value = NA,
                      P_valB = NA,
                      P_valR = NA,
                      stringsAsFactors = F)

# TraitData: data frame with binary/continuous values
# TraitTreeDist: data frame of patristic trait distances
# TraitPhy: phylo that TraitTreeDist is derived from
SLATest <- phylosig(TraitPhy, SLA, test = T)
outData[outData$Trait == 'SLA', 'Value'] <- SLATest$K
outData[outData$Trait == 'SLA', 'P_valB'] <- SLATest$P
outData[outData$Trait == 'SLA', 'N'] <- length(SLA[!is.na(SLA)])

Ht <- phylosig(TraitPhy, Height, test = T)
outData[outData$Trait == 'Height', 'Value'] <- Ht$K
outData[outData$Trait == 'Height', 'P_valB'] <- Ht$P
outData[outData$Trait == 'Height', 'N'] <- length(Height[!is.na(Height)])

ToughTest <- phylosig(TraitPhy, Tough, test = T)
outData[outData$Trait == 'Leaf toughness', 'Value'] <- ToughTest$K
outData[outData$Trait == 'Leaf toughness', 'P_valB'] <- ToughTest$P
outData[outData$Trait == 'Leaf toughness', 'N'] <- length(Tough[!is.na(Tough)])

WoodDensTest <- phylosig(TraitPhy, WoodDens, test = T)
outData[outData$Trait == 'Stem density', 'Value'] <- WoodDensTest$K
outData[outData$Trait == 'Stem density', 'P_valB'] <- WoodDensTest$P
outData[outData$Trait == 'Stem density', 'N'] <- length(WoodDens[!is.na(WoodDens)])

# Binary traits for testing: Growth Form levels, N_fixer, Woody, dispersal
binTraits <- TraitData[ ,c(1, 7,9:25)]

TraitPhy$node.label <- 1:length(TraitPhy$node.label)

source('R_new/SE_Phylo_D.R')

for(i in names(binTraits)[-c(1)]){
 PhyD <- SE_phylo.d(binTraits, TraitPhy, "Species.Name", i)
 outData[outData$Trait == i, 'Value'] <- PhyD$DEstimate
 outData[outData$Trait == i, 'P_valB'] <- PhyD$Pval0
 outData[outData$Trait == i, 'P_valR'] <- PhyD$Pval1
 outData[outData$Trait == i, 'N'] <- length(TraitData[ ,i][!is.na(TraitData[ ,i])])
}



# Next, i'll do a mantel test with the trait and phylogenetic distances
# using square root of phylo distances per Hardy+Pavoine 2012
# Additionally, I'm going to create a distance matrix for flower period
# and try a Mantel test on that as well.

Circ_to_DistMat <- function(x, y, M) {
  out <- sqrt(1 - abs(1 - 2 * abs(x/M - y/M)))
  return(out)
}

FlowerPeriod <- TraitData$Flower.Period
names(FlowerPeriod) <- TraitData$Species.Name

FlowerPeriodDist <- outer(FlowerPeriod, 
                          FlowerPeriod,
                          FUN = Circ_to_DistMat,
                          M = 12)

names(FlowerPeriodDist) <- rownames(FlowerPeriodDist)

regTraitDist <- make_regional_trait_dist(TraitData,
                                         names(TraitData)[-c(1,5)]) %>%
                as.matrix()
traitPhyDists <- cophenetic(TraitPhy) %>% sqrt

regTraitDist <- regTraitDist[rownames(regTraitDist) %in% rownames(traitPhyDists),
                             colnames(regTraitDist) %in% colnames(traitPhyDists)] %>%
                .[rownames(traitPhyDists), colnames(traitPhyDists)]

FlowerPeriodDist <- FlowerPeriodDist[rownames(FlowerPeriodDist) %in% rownames(traitPhyDists),
                             colnames(FlowerPeriodDist) %in% colnames(traitPhyDists)] %>%
  .[rownames(traitPhyDists), colnames(traitPhyDists)]

if(identical(rownames(regTraitDist), rownames(traitPhyDists)) &
   identical(colnames(regTraitDist), colnames(traitPhyDists))){
 
  allResult <- mantel(traitPhyDists, regTraitDist)
  outData[outData$Trait == 'All combined', 'Value'] <- allResult$statistic
  outData[outData$Trait == 'All combined', 'P_valB'] <- allResult$signif
  outData[outData$Trait == 'All combined', 'N'] <- dim(traitPhyDists)[1]
}

if(identical(rownames(FlowerPeriodDist), rownames(traitPhyDists)) &
   identical(colnames(FlowerPeriodDist), colnames(traitPhyDists))){
  
  allResult <- mantel(traitPhyDists, FlowerPeriodDist)
  outData[outData$Trait == 'Flower Period', 'Value'] <- allResult$statistic
  outData[outData$Trait == 'Flower Period', 'P_valB'] <- allResult$signif
  outData[outData$Trait == 'Flower Period', 'N'] <- dim(FlowerPeriodDist)[1]
}


outData$Value <- round(outData$Value, 4)

add_stars <- function(x){
  if(is.na(x)){
    x <- NA
  } else if(x <= .001){
    x <- paste0(x, "***")
  } else if(x <= .01 & x > .001){
    x <- paste0(x, "**")
  } else if(x <= .05 & x > .01){
    x <- paste0(x, '*')
  } else if(x <= .1 & x > .05){
    x <- paste0(x, '+')
  } else {
    x <- paste(x)
  }
}

outData$P_valB <- sapply(outData$P_valB, FUN = function(x) add_stars(x))
outData$P_valR <- sapply(outData$P_valR, FUN = function(x) add_stars(x))

# write.csv(outData, '../Figures/Phylo_Signal_in_Functional_Traits.csv', na = "",
#            row.names = F)
