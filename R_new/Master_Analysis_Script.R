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

demo.data$MEPPInv <- NA
demo.data$MPD <- NA
demo.data$NND <- NA
demo.data$AWMPD <- NA
demo.data$AWNND <- NA

for(x in unique(demo.data$Species)){
  # make local phylo and functional distance matrices. The functional
  # distance matrices are really just place holders, they won't be used at
  # all because a = 1. I have to make them though because the function I wrote 
  # to do this requires a functional distance matrix too because I'm an idiot.
  phyloMat <- make_local_phylo_dist(x, communities, phylo)
  funMat <- make_local_trait_dist(x, communities, 
                                  trait.data = tyson$traits,
                                  traits = names(tyson$traits)[-1],
                                  scale = 'scaledBYrange')
  
  # Abundance weighted data
  AWDat <- rarefy_FPD(x, phylo.mat = phyloMat,
                      fun.mat = funMat,
                      n.rare = 11, a = 1,
                      p = 2, abundance.weighted = T,
                      community.data = communities)
  # Unweighted data
  UWDat <- rarefy_FPD(x, phylo.mat = phyloMat,
                      fun.mat = funMat,
                      n.rare = 11, a = 1,
                      p = 2, abundance.weighted = F,
                      community.data = NULL)  
  
  # extract all data, as well as MEPP invasive rating for focal species
  demo.data[demo.data$Species == x, 'MPD'] <- UWDat$rare.mpd
  demo.data[demo.data$Species == x, 'NND'] <- UWDat$rare.nnd
  demo.data[demo.data$Species == x, 'AWMPD'] <- AWDat$rare.mpd
  demo.data[demo.data$Species == x, 'AWNND'] <- AWDat$rare.nnd
  demo.data[demo.data$Species == x, 
            'MEPPInv'] <- communities[communities$exotic_species == x &
                                      communities$community == x,
                                                              'Invasive']
  
}

# log transform abundance weighted data, but retain original values too
demo.data$logAWMPD <- log(demo.data$AWMPD)
demo.data$logAWNND <- log(demo.data$AWNND)

# regressions!
mpdLM <- lm(ESCR2 ~ MPD + CRBM, data = demo.data)
nndLM <- lm(ESCR2 ~ NND + CRBM, data = demo.data)
awmpdLM <- lm(ESCR2 ~ AWMPD + CRBM, data = demo.data)
awnndLM <- lm(ESCR2 ~ AWNND + CRBM, data = demo.data)
logAWmpdLM <- lm(ESCR2 ~ logAWMPD + CRBM, data = demo.data)
logAWnndLM <- lm(ESCR2 ~ logAWNND + CRBM, data = demo.data)

# prepare model outputs for table format
mpdCT <- summary(mpdLM)$coefficients %>% data.frame() 
nndCT <- summary(nndLM)$coefficients %>% data.frame()
awmpdCT <- summary(awmpdLM)$coefficients %>% data.frame()
awnndCT <- summary(awnndLM)$coefficients %>% data.frame()
logAWmpdCT <- summary(logAWmpdLM)$coefficients %>% data.frame()
logAWnndCT <- summary(logAWnndLM)$coefficients %>% data.frame()

# more table formatting BS
modelTable <- rbind(mpdCT, nndCT, awmpdCT,awnndCT,
                    logAWmpdCT, logAWnndCT)
Parameter <- rownames(modelTable)
Scale <- rep(c("Local", NA, NA), 6)
Type <- rep(c("Linear", NA, NA), 6)
ResponseVariable <- rep(c('Log-Response Ratio', NA, NA), 6)
R2adj <- c(summary(mpdLM)$adj.r.squared, NA, NA,
           summary(nndLM)$adj.r.squared, NA, NA,
           summary(awmpdLM)$adj.r.squared, NA, NA,
           summary(awnndLM)$adj.r.squared, NA, NA,
           summary(logAWmpdLM)$adj.r.squared, NA, NA,
           summary(logAWnndLM)$adj.r.squared, NA, NA)

# combine all data, clean up parameter column
names(modelTable)[2:4] <- c('Std.Error', 'Test.Statistic', 'p.value')
modelTable <- cbind(Scale, Type, ResponseVariable, Parameter,
                    modelTable, R2adj,
                    stringsAsFactors = F)


modelTable$Parameter <- str_replace_all(modelTable$Parameter,
                                        '[:punct:]|[:digit:]', "")

modelTable[ ,5:9] <- round(modelTable[ ,5:9], 4)

add_stars <- function(x){
  if(x < .001){
    x <- paste0(x, "***")
  } else if(x < .01 & x > .001){
    x <- paste0(x, "**")
  } else if(x < .05 & x > .01){
    x <- paste0(x, '*')
  } else if(x < .1){
    x <- paste0(x, '+')
  } else {
    x <- paste(x)
  }
}
  
# add significance levels
modelTable$p.value <- lapply(modelTable$p.value,
                             FUN = function (x) add_stars(x)) %>% unlist()

# write table to csv so we can include in paper
# write.csv(modelTable, '../Figures/Phylo_models_output.csv', na = "",
          # row.names = F)

# plots. Need to choose which ones we actually want to use
par(mfrow = c(3,2))
plot(ESCR2~MPD, data = demo.data,
     type = 'n', xlab = 'MPD', ylab = 'Effect size of competition')
points(demo.data$MPD[demo.data$MEPPInv == 0],
       demo.data$ESCR2[demo.data$MEPPInv == 0],
       col = 'blue', pch = 19)
points(demo.data$MPD[demo.data$MEPPInv == 1],
       demo.data$ESCR2[demo.data$MEPPInv == 1],
       col = 'red', pch = 19)
abline(mpdLM, col = 'orange', lty = 2)

plot(ESCR2~NND, data = demo.data, 
     type = 'n', xlab = 'NND', ylab = 'Effect size of competition')
points(demo.data$NND[demo.data$MEPPInv == 0],
       demo.data$ESCR2[demo.data$MEPPInv == 0],
       col = 'blue', pch = 19)
points(demo.data$NND[demo.data$MEPPInv == 1],
       demo.data$ESCR2[demo.data$MEPPInv == 1],
       col = 'red', pch = 19)
abline(nndLM, col = 'orange', lty = 2)

plot(ESCR2~AWMPD, data = demo.data,
     type = 'n', xlab = 'AWMPD', ylab = 'Effect size of competition')
points(demo.data$AWMPD[demo.data$MEPPInv == 0],
       demo.data$ESCR2[demo.data$MEPPInv == 0],
       col = 'blue', pch = 19)
points(demo.data$AWMPD[demo.data$MEPPInv == 1],
       demo.data$ESCR2[demo.data$MEPPInv == 1],
       col = 'red', pch = 19)
abline(awmpdLM, col = 'orange', lty = 2)

plot(ESCR2~AWNND, data = demo.data,
     type = 'n', xlab = 'AWNND', ylab = 'Effect size of competition')
points(demo.data$AWNND[demo.data$MEPPInv == 0],
       demo.data$ESCR2[demo.data$MEPPInv == 0],
       col = 'blue', pch = 19)
points(demo.data$AWNND[demo.data$MEPPInv == 1],
       demo.data$ESCR2[demo.data$MEPPInv == 1],
       col = 'red', pch = 19)
abline(awnndLM, col = 'orange', lty = 2)

plot(ESCR2~logAWMPD, data = demo.data,
     type = 'n', xlab = 'logAWmpd', ylab = 'Effect size of competition')
points(demo.data$logAWMPD[demo.data$MEPPInv == 0],
       demo.data$ESCR2[demo.data$MEPPInv == 0],
       col = 'blue', pch = 19)
points(demo.data$logAWMPD[demo.data$MEPPInv == 1],
       demo.data$ESCR2[demo.data$MEPPInv == 1],
       col = 'red', pch = 19)
abline(logAWmpdLM, col = 'orange', lty = 2)

plot(ESCR2~logAWNND, data = demo.data,
     type = 'n', xlab = 'logAWnnd', ylab = 'Effect size of competition')
points(demo.data$logAWNND[demo.data$MEPPInv == 0],
       demo.data$ESCR2[demo.data$MEPPInv == 0],
       col = 'blue', pch = 19)
points(demo.data$logAWNND[demo.data$MEPPInv == 1],
       demo.data$ESCR2[demo.data$MEPPInv == 1],
       col = 'red', pch = 19)
abline(logAWnndLM, col = 'orange', lty = 2)


# functional phylogenetic regressions-------------
# next, determine best models using latest version of invasives FPD
# since abundances don't seem to matter much, I'll leave those out (for
# sake of computing time) and use the shiny app to find R2~a curve with highest
# R2

# ----- leaving RStudio Beep Boop, report back w/ best models-----



