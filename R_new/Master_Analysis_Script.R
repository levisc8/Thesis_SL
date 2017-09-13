# master script w/version control---------------------
# so I don't keep renaming things
# depends on FunPhylo, pez, dplyr, and some others

rm(list = ls(all = T))
library(FunPhylo)
library(pez)
library(dplyr)

data(tyson)

spp.list <- tyson$spp.list
for(i in 2:4){
  spp.list[ ,i] <- as.numeric(spp.list[ ,i])
}
phylo <- tyson$phylo
communities <- tyson$communities
demo.data <- tyson$demo.data


# regional scale analysis-----------
# create Tyson scale phylo distance matrix
regionalTysonDists <- cophenetic(phylo) %>% data.frame()

diag(regionalTysonDists) <- NA
regionalTysonDists <- regionalTysonDists/max(regionalTysonDists, na.rm = T)
# next, create extract data for all exotic species. Remove cerastium_spp. It is
# included in the species list so that it matches with the traits data, but it 
# is also probably covered by one of the three Cerastium species on the species list.
exotics <- filter(spp.list, Exotic == 1 & Species != 'Cerastium_spp.')

exotics$MPD <- NA
exotics$NND <- NA
exotics$Focal <- NA

# extract mpd and nnd for all exotics
for(x in unique(exotics$Species)){
  exotics[exotics$Species == x, 'MPD'] <- mean(regionalTysonDists[ ,x], na.rm = T)
  exotics[exotics$Species == x, 'NND'] <- min(regionalTysonDists[ ,x], na.rm = T)
  exotics[exotics$Species == x, 'Focal'] <- ifelse(x %in% demo.data$Species,
                                                   1, 0)
}

# models
nndGLM <- glm(Invasive ~ NND, data = exotics,
              family = binomial())
mpdGLM <- glm(Invasive ~ MPD, data = exotics,
              family = binomial())

# results
summary(nndGLM)
summary(mpdGLM)

# regional scale figures----------------------
# plot results
xx <- seq(0, 1, .1)

plot(exotics$MPD, exotics$Invasive, 
     ylab = 'Pr(Invasive)', xlab = 'Mean Pairwise Distance',
     type = 'n')
points(exotics$MPD[exotics$Focal == 0],
       exotics$Invasive[exotics$Focal == 0],
       col = 'blue', pch = 19)
points(exotics$MPD[exotics$Focal == 1],
       exotics$Invasive[exotics$Focal == 1],
       col = 'red', pch = 19, cex = 1.2)
legend('bottomright', c('Focal Species', 'Other Exotics @ Tyson'),
       col = c('red', 'blue'), pch = 19)
lines(xx, predict(mpdGLM, data.frame(MPD = xx), type = 'response'),
      lty = 2, lwd = 1, col = 'red')

plot(exotics$NND, exotics$Invasive, 
     ylab = 'Pr(Invasive)', xlab = 'Nearest Neighbor Distance',
     type = 'n')
points(exotics$NND[exotics$Focal == 0],
       exotics$Invasive[exotics$Focal == 0],
       col = 'blue', pch = 19)
points(exotics$NND[exotics$Focal == 1],
       exotics$Invasive[exotics$Focal == 1],
       col = 'red', pch = 19, cex = 1.2)
legend('bottomright', c('Focal Species', 'Other Exotics @ Tyson'),
       col = c('red', 'blue'), pch = 19)
lines(xx, predict(nndGLM, data.frame(NND = xx), type = 'response'),
      lty = 2, lwd = 1, col = 'red')

# local scale effect sizes, phylo only---------------
# time for effect size of competition! phylogenetic only to start
# we'll get to the traits later

