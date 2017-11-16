# master script w/version control---------------------
# so I don't keep renaming things
# depends on FunPhylo, pez, dplyr, and some others

rm(list = ls(all = T))
library(FunPhylo)
library(pez)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(tidyr)
library(vegan)
 
data(tyson)

# do a little data munging and variable conversion
spp.list <- tyson$spp.list
for(i in 2:4){
  spp.list[ ,i] <- as.numeric(spp.list[ ,i])
}
# remove allium vineale since it's a monocot. However, I wanted to keep
# it in the data set to be distributed because someone else might find it
# useful. We also have trait data for it
phylo <- tyson$phylo
communities <- tyson$communities %>% filter(community != 'Allium_vineale')
demo.data <- tyson$demo.data

# regional scale analysis-----------
# create Tyson scale phylo distance matrix
regionalTysonDists <- cophenetic(phylo) %>% data.frame()

diag(regionalTysonDists) <- NA

# next, create extract data for all exotic species. Remove cerastium_spp. It is
# included in the species list so that it matches with the traits data, but it 
# is also probably covered by one of the three Cerastium species on the species list.
exotics <- filter(spp.list, Exotic == 1 & Species != 'Cerastium_spp.')

exotics$MPD <- NA
exotics$NND <- NA
exotics$Focal <- NA

# extract mpd and nnd for all exotics
for(x in unique(exotics$Species)){
  exotics[exotics$Species == x, 'MPD'] <- mean(regionalTysonDists[ ,x],
                                               na.rm = T)
  exotics[exotics$Species == x, 'NND'] <- min(regionalTysonDists[ ,x],
                                              na.rm = T)
  exotics[exotics$Species == x, 'Focal'] <- ifelse(x %in% demo.data$Species,
                                                   'Y', 
                                                   'N')
}

# models
nndGLM <- glm(Invasive ~ NND, data = exotics,
              family = binomial())
mpdGLM <- glm(Invasive ~ MPD, data = exotics,
              family = binomial())

# results
summary(nndGLM)
summary(mpdGLM)

plt.blank <- theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   axis.line = element_line(colour = "black"),
                   axis.title = (element_text(size = 14)))

regional.nnd.plt <- ggplot(data = exotics, aes(x = NND, y = Invasive)) +
  geom_jitter(aes(color = Focal), alpha = 0.4,
              width = 0, height = 0.075,
              show.legend = FALSE) + 
  stat_smooth(method = 'glm', method.args = list(family = 'binomial'),
              se = TRUE, alpha = 0.1, color = 'grey') + 
  scale_color_manual(values = c('blue', 'red')) +
  xlab('') +
  scale_y_continuous('',
                     breaks = seq(0, 1, 1),
                     limits = c(-0.15, 1.15)) +
  scale_x_continuous('', 
                     breaks = seq(0, max(exotics$NND, na.rm = T),
                                  40),
                     limits = c(0, max(exotics$NND, na.rm = T) + 5)) +
  plt.blank +
  annotate('text', x = 177, y = 1.15, label = 'A',
           size = 5)

regional.mpd.plt <- ggplot(data = exotics, aes(x = MPD, y = Invasive)) +
  geom_jitter(aes(color = Focal), alpha = 0.4,
              width = 0, height = 0.075,
              show.legend = FALSE) + 
  stat_smooth(method = 'glm', method.args = list(family = 'binomial'),
              se = TRUE, alpha = 0.1, color = 'grey') + 
  scale_color_manual(values = c('blue', 'red')) +
  xlab('') +
  scale_y_continuous('',
                     breaks = seq(0, 1, 1),
                     limits = c(-0.15, 1.15)) +
  scale_x_continuous('', 
                     breaks = seq(200, 280, 20),
                     limits = c(200,
                                max(exotics$MPD, na.rm = T) + 5)) +
  plt.blank +
  annotate('text', x = 295, y = 1.1, label = 'B',
           size = 5)

# regional scale figures----------------------
# plot results
# xx <- seq(0, 1, .1)
# 
# plot(exotics$MPD, exotics$Invasive, 
#      ylab = 'Pr(Invasive)', xlab = 'Mean Pairwise Distance',
#      type = 'n')
# points(exotics$MPD[exotics$Focal == 0],
#        exotics$Invasive[exotics$Focal == 0],
#        col = 'blue', pch = 19)
# points(exotics$MPD[exotics$Focal == 1],
#        exotics$Invasive[exotics$Focal == 1],
#        col = 'red', pch = 19, cex = 1.2)
# legend('bottomright', c('Focal Species', 'Other Exotics @ Tyson'),
#        col = c('red', 'blue'), pch = 19)
# lines(xx, predict(mpdGLM, data.frame(MPD = xx), type = 'response'),
#       lty = 2, lwd = 1, col = 'red')
# 
# plot(exotics$NND, exotics$Invasive, 
#      ylab = 'Pr(Invasive)', xlab = 'Nearest Neighbor Distance',
#      type = 'n')
# points(exotics$NND[exotics$Focal == 0],
#        exotics$Invasive[exotics$Focal == 0],
#        col = 'blue', pch = 19)
# points(exotics$NND[exotics$Focal == 1],
#        exotics$Invasive[exotics$Focal == 1],
#        col = 'red', pch = 19, cex = 1.2)
# legend('bottomright', c('Focal Species', 'Other Exotics @ Tyson'),
#        col = c('red', 'blue'), pch = 19)
# lines(xx, predict(nndGLM, data.frame(NND = xx), type = 'response'),
#       lty = 2, lwd = 1, col = 'red')

# local scale effect sizes, phylo only---------------
# time for effect size of competition! phylogenetic only to start
# we'll get to the traits later

communities$MPD.Local <- NA
communities$NND.Local <- NA

All.Local.exotics <- subset(communities, exotic_species == '')
All.Local.exotics$Focal <- NA

demo.data$MEPPInv <- NA
demo.data$MPD <- NA
demo.data$NND <- NA
demo.data$AWMPD <- NA
demo.data$AWNND <- NA
demo.data$Regional_MPD <- NA
demo.data$Regional_NND <- NA

for(x in unique(demo.data$Species)){
  cat('Crunching data for species: ', x, '\n')
  # make local phylo and functional distance matrices. The functional
  # distance matrices are really just place holders, they won't be used at
  # all because a = 1. I have to make them though because the function I wrote 
  # to do this requires a functional distance matrix too because I'm an idiot.
  demo.data[demo.data$Species == x, 'Regional_MPD'] <- exotics[exotics$Species == x,
                                                               'MPD']
  demo.data[demo.data$Species == x, 'Regional_NND'] <- exotics[exotics$Species == x,
                                                               'NND']
  # Per Masha's suggestion, try logistic regression a la regional analysis,
  # but for local communities. Is presence/absence related to novelty?
  local.exotics <- filter(communities, exotic_species == x &
                            alien == 1)
  
  phyloMat <- make_local_phylo_dist(x, communities, phylo)
  funMat <- make_local_trait_dist(x, communities, 
                                  trait.data = tyson$traits,
                                  traits = names(tyson$traits)[-1],
                                  scale = 'scaledBYrange')
  # This is not the most efficient code...
  # Loc.Mat <- phyloMat/max(phyloMat)
  Loc.Mat <- phyloMat
  diag(Loc.Mat) <- NA  
  
  for(y in unique(local.exotics$community)){
    local.exotics[local.exotics$community == y,
                  'MPD.Local'] <- mean(Loc.Mat[ ,y], na.rm = TRUE)
    local.exotics[local.exotics$community == y,
                  'NND.Local'] <- min(Loc.Mat[ ,y], na.rm = TRUE)
    local.exotics[local.exotics$community == y,
                  'Focal'] <- exotics[exotics$Species == y, 'Focal']
  }
  All.Local.exotics <- rbind(All.Local.exotics, local.exotics)
  
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

# Repeat invasive glm analysis, but at local spatial scale-----------
Local.Ex.For.Reg <- All.Local.exotics %>% group_by(community) %>%
  summarise(Invasive = mean(Invasive),
            Focal = first(Focal),
            MPD.Local = mean(MPD.Local),
            NND.Local = mean(NND.Local))

summary(Local.Ex.For.Reg$NND.Local)

LocalMPDInvGLM <- glm(Invasive ~ MPD.Local, family = binomial(),
                      data = Local.Ex.For.Reg)
LocalNNDInvGLM <- glm(Invasive ~ NND.Local, family = binomial(),
                      data = Local.Ex.For.Reg)

message('Local Scale Regional Style Analysis')
summary(LocalMPDInvGLM)
summary(LocalNNDInvGLM)

local.nnd.plt <- ggplot(data = Local.Ex.For.Reg,
                        aes(x = NND.Local, y = Invasive)) +
  geom_jitter(aes(color = Focal), alpha = 0.4,
              width = 0, height = 0.075,
              show.legend = FALSE) + 
  # stat_smooth(method = 'glm', method.args = list(family = 'binomial'),
  #             se = TRUE, alpha = 0.2) + 
  scale_color_manual(values = c('blue', 'red')) +
  xlab('Nearest Neighbor Distance') +
  scale_y_continuous('',
                     breaks = seq(0, 1, 1),
                     limits = c(-0.15, 1.15)) +
  scale_x_continuous('Nearest Neighbor Distance', 
                     breaks = seq(0, max(Local.Ex.For.Reg$NND.Local),
                                  60),
                     limits = c(0, max(Local.Ex.For.Reg$NND.Local))) +
  plt.blank +
  annotate('text', x = 240, y = 1.15, label = 'C',
           size = 5)

local.mpd.plt <- ggplot(data = Local.Ex.For.Reg,
                        aes(x = MPD.Local, y = Invasive)) +
  geom_jitter(aes(color = Focal), alpha = 0.4,
              width = 0, height = 0.075,
              show.legend = FALSE) + 
  # stat_smooth(method = 'glm', method.args = list(family = 'binomial'),
  #             se = TRUE, alpha = 0.2) + 
  scale_color_manual(values = c('blue', 'red')) +
  xlab('Nearest Neighbor Distance') +
  scale_y_continuous('',
                     breaks = seq(0, 1, 1),
                     limits = c(-0.15, 1.15)) +
  scale_x_continuous('Mean Pairwise Distance', 
                     breaks = seq(140,
                                  280,
                                  40),
                     limits = c(min(Local.Ex.For.Reg$MPD.Local),
                                max(Local.Ex.For.Reg$MPD.Local))) +
  plt.blank +
  annotate('text', x = 275, y = 1.15, label = 'D',
           size = 5)

ggdraw() +
  draw_plot(regional.nnd.plt, x = .1, y = .5, 
            height = .45, width = .45) +
  draw_plot(regional.mpd.plt, x = .55, y = .5,
            height = .45, width = .45) + 
  draw_plot(local.nnd.plt, x = .1, y = 0.05,
            height = .45, width = .45) +
  draw_plot(local.mpd.plt,
            x = .55, y = 0.05,
            height = .45, width = .45) + 
  annotate('text', x = .1, y = .05, label = 'Local Scale',
           size = 5) +
  annotate('text', x = .1, .981, label = 'Regional Scale',
           size = 5) + 
  annotate('text', x = .08, y = .55, label = 'Focal Species') +
  annotate('text', x = .055, y = .522, label = 'Yes') + 
  annotate('point', x = .025, y = .522, color = 'red', alpha = .2) + 
  annotate('text', x = .055, y = .492, label = 'No') + 
  annotate('point', x = .025, y = .492, color = 'blue', alpha = .2) +
  annotate('text', x = .105, y = .885, label = 'Invasive') +
  annotate('text', x = .11, y = .637, label = 'Exotic') +
  annotate('text', x = .105, y = .435, label = 'Invasive') +
  annotate('text', x = .11, y = .192, label = 'Exotic')

# local ESCR regressions-------------
mpdLM <- lm(ESCR2 ~ MPD + CRBM, data = demo.data)
nndLM <- lm(ESCR2 ~ NND + CRBM, data = demo.data)
logAWmpdLM <- lm(ESCR2 ~ logAWMPD + CRBM, data = demo.data)
logAWnndLM <- lm(ESCR2 ~ logAWNND + CRBM, data = demo.data)
RegMPDLM <- lm(ESCR2 ~ Regional_MPD + CRBM, data = demo.data)
RegNNDLM <- lm(ESCR2 ~ Regional_NND + CRBM, data = demo.data)


message('standardized biomass')
summary(mpdLM)
summary(nndLM)
summary(logAWmpdLM)
summary(logAWnndLM)
summary(RegMPDLM)
summary(RegNNDLM)

# prepare model outputs for table format
mpdCT <- summary(mpdLM)$coefficients %>% data.frame() 
nndCT <- summary(nndLM)$coefficients %>% data.frame()
logAWmpdCT <- summary(logAWmpdLM)$coefficients %>% data.frame()
logAWnndCT <- summary(logAWnndLM)$coefficients %>% data.frame()
regMPDCT <- summary(RegMPDLM)$coefficients %>% data.frame()
regNNDCT <- summary(RegNNDLM)$coefficients %>% data.frame()

# more table formatting BS
modelTable <- rbind(mpdCT, nndCT,
                    logAWmpdCT, logAWnndCT,
                    regMPDCT, regNNDCT)
Parameter <- rownames(modelTable)
Scale <- c(rep(c("Local", NA, NA), 4),
           rep(c('Regional', NA, NA), 2))
Type <- rep(c("Linear", NA, NA), 6)
ResponseVariable <- rep(c('Log-Response Ratio', NA, NA), 6)
R2adj <- c(summary(mpdLM)$adj.r.squared, NA, NA,
           summary(nndLM)$adj.r.squared, NA, NA,
           summary(logAWmpdLM)$adj.r.squared, NA, NA,
           summary(logAWnndLM)$adj.r.squared, NA, NA,
           summary(RegMPDLM)$adj.r.squared, NA, NA,
           summary(RegNNDLM)$adj.r.squared, NA, NA)

# combine all data, clean up parameter column
names(modelTable)[2:4] <- c('Std.Error', 'Test.Statistic', 'p.value')
modelTable <- cbind(Scale, Type, ResponseVariable, Parameter,
                    modelTable, R2adj,
                    stringsAsFactors = F)


modelTable$Parameter <- stringr::str_replace_all(modelTable$Parameter,
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

#  write table to csv so we can include in paper
# write.csv(modelTable, '../Figures/Phylo_models_output.csv', na = "",
#             row.names = F)
# 
# # make model table for logistic regressions
# regmpdCT <- summary(mpdGLM)$coefficients %>% data.frame() 
# regnndCT <- summary(nndGLM)$coefficients %>% data.frame()
# locmpdCT <- summary(LocalMPDInvGLM)$coefficients %>% data.frame()
# locnndCT <- summary(LocalNNDInvGLM)$coefficients %>% data.frame()
# 
# logisticModelTable <- rbind(regmpdCT,
#                             regnndCT,
#                             locmpdCT,
#                             locnndCT)
# Parameter <- rownames(logisticModelTable)
# Scale <- c(rep(c('Regional', NA), 2),
#            rep(c('Local', NA), 2))
# Type <- rep(c('Logistic',NA), 4)
# ResponseVariable <- rep(c('Invasive Classification', NA), 4)
# 
# names(logisticModelTable)[2:4] <- c('Std.Error', 'Test.Statistic', 'p.value')
# 
# logisticModelTable <- cbind(Scale,
#                             Type, 
#                             ResponseVariable, 
#                             Parameter,
#                             logisticModelTable,
#                             stringsAsFactors = F)
# 
# logisticModelTable$Parameter <- stringr::str_replace_all(logisticModelTable$Parameter,
#                                                          '[:punct:]|[:digit:]', "")
# 
# logisticModelTable[ ,5:8] <- round(logisticModelTable[ ,5:8], 4)
# logisticModelTable$p.value <- lapply(logisticModelTable$p.value,
#                                      FUN = function (x) add_stars(x)) %>% unlist()
# write.csv(logisticModelTable, file = '../Figures/Logistic_regression_outputs.csv',
#           na = "", row.names = FALSE)

# plots! build first, draw later
bad.mets <- c('AWNND','AWMPD')
forPlot <- gather(demo.data, Metric, Magnitude, MPD:logAWNND) %>% 
  filter(!Metric %in% bad.mets )
forPlot$MEPPInv <- ifelse(forPlot$MEPPInv == 1, 'Y', 'N')
 
# create column of predictions for each regression that needs them
# ggplot will not draw trend lines with covariates, so we have to make
# them ourselves :(
x <- seq(min(demo.data$MPD), max(demo.data$MPD), length.out = 14)
y <- seq(min(demo.data$CRBM), max(demo.data$CRBM), length.out = 14)
MPD.Preds <- predict(mpdLM, data.frame(MPD = x,
                                       CRBM = y),
                     type = 'response',
                     interval = 'confidence',
                     se.fit = TRUE)$fit %>%
  data.frame %>% 
  cbind(., x)

x <- seq(min(demo.data$NND), max(demo.data$NN), length.out = 14)
NND.Preds <- predict(nndLM, data.frame(NND = x,
                                      CRBM = y),
                     type = 'response',
                     interval = 'confidence',
                     se.fit = TRUE)$fit %>%
  data.frame %>% 
  cbind(., x)

x <- seq(min(demo.data$logAWMPD), max(demo.data$logAWMPD), length.out = 14)
AWMPD.Pred <- predict(logAWmpdLM, data.frame(logAWMPD = x,
                                             CRBM = y),
                      type = 'response',
                      interval = 'confidence',
                      se.fit = TRUE)$fit %>%
  data.frame %>% 
  cbind(., x)

x <- seq(min(demo.data$logAWNND), max(demo.data$logAWNND), length.out = 14)
AWNND.Pred <- predict(logAWnndLM, data.frame(logAWNND = x,
                                             CRBM = y),
                      type = 'response',
                      interval = 'confidence',
                      se.fit = TRUE)$fit %>%
  data.frame %>% 
  cbind(., x)

loc.lrr.mpd.plt <- ggplot(data = filter(forPlot, Metric == 'MPD'),
                          aes(x = Magnitude,
                              y = ESCR2)) +
  scale_x_continuous('MPD', 
                     breaks = seq(140,
                                  260,
                                  40),
                     limits = c(135, 260)) +
  scale_y_continuous('',
                     breaks = seq(0, 3, 1),
                     limits = c(-.5, 3)) +
  geom_point(alpha = .4, aes(color = MEPPInv),
             show.legend = FALSE,
             size = 2) + 
  scale_color_manual(values = c('red','blue')) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  geom_line(data = MPD.Preds, 
            aes(x = x, 
                y = fit),
            color = 'black',
            alpha = .8,
            size = .8) 

loc.lrr.nnd.plt <- ggplot(data = filter(forPlot, Metric == 'NND'),
                          aes(x = Magnitude, y = ESCR2)) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE,
             size = 2) +
  scale_x_continuous('NND', 
                     breaks = seq(0,
                                  240,
                                  60),
                     limits = c(0, 240)) +
  scale_y_continuous('',
                     breaks = seq(0, 3, 1),
                     limits = c(-.5, 3)) +
  scale_color_manual(values = c('red','blue')) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  geom_line(data = NND.Preds, 
            aes(x = x, 
                y = fit),
            color = 'black',
            alpha = .8,
            size = .8) 


loc.lrr.aw.mpd.plt <- ggplot(data = filter(forPlot, Metric == 'logAWMPD'),
                             aes(x = Magnitude, y = ESCR2)) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE,
             size = 2) +
  scale_x_continuous('log(abundance-weighted MPD)', 
                     breaks = seq(4.5,
                                  6,
                                  .3),
                     limits = c(4.5, 6)) +
  scale_y_continuous('',
                     breaks = seq(0, 3, 1),
                     limits = c(-.5, 3)) +
  scale_color_manual(values = c('red','blue')) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  geom_line(data = AWMPD.Pred, 
            aes(x = x, 
                y = fit),
            color = 'black',
            alpha = .8,
            size = .8) 


loc.lrr.aw.nnd.plt <- ggplot(data = filter(forPlot, Metric == 'logAWNND'),
                             aes(x = Magnitude, y = ESCR2)) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE,
             size = 2) +
  scale_x_continuous('log(abundance-weighted NND)', 
                     breaks = seq(-3,
                                  2,
                                  .8),
                     limits = c(-3.1, 2)) +
  scale_y_continuous('',
                     breaks = seq(0, 3, 1),
                     limits = c(-.5, 3)) +
  scale_color_manual(values = c('red','blue')) +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  geom_line(data = AWNND.Pred, 
            aes(x = x, 
                y = fit),
            color = 'black',
            alpha = .8,
            size = .8) 


reg.lrr.mpd.plt <- ggplot(data = filter(forPlot, Metric == 'Regional_MPD'),
                             aes(x = Magnitude, y = ESCR2)) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE,
             size = 2) +
  scale_x_continuous('MPD', 
                     breaks = seq(200,
                                  250,
                                  10),
                     limits = c(200, 250)) +
  scale_y_continuous('',
                     breaks = seq(0, 3, 1),
                     limits = c(-.5, 3)) +
  scale_color_manual(values = c('red','blue')) +
  geom_hline(yintercept = 0, linetype = 'dotted')


reg.lrr.nnd.plt <- ggplot(data = filter(forPlot, Metric == 'Regional_NND'),
                          aes(x = Magnitude, y = ESCR2)) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE,
             size = 2) +
  scale_x_continuous('NND', 
                     breaks = seq(0,
                                  100,
                                  20),
                     limits = c(0, 110)) +
  scale_y_continuous('',
                     breaks = seq(0, 3, 1),
                     limits = c(-.5, 3)) +
  scale_color_manual(values = c('red','blue')) +
  geom_hline(yintercept = 0, 
             linetype = 'dotted')

# create blank canvas and get to drawing!

ggdraw() +
  draw_plot(loc.lrr.mpd.plt, x = .1, y = .67,
            width = .45, height = .33) + 
  draw_plot(loc.lrr.nnd.plt, x = .55, y = .67,
            width = .45, height = .33) +
  draw_plot(loc.lrr.aw.mpd.plt, x = .1, y = .34,
            width = .45, height = .33) +
  draw_plot(loc.lrr.aw.nnd.plt, x = .55, y = .34,
            width = .45, height = .33) +
  draw_plot(reg.lrr.mpd.plt, x = .1, y = 0,
            width = .45, height = .33) +
  draw_plot(reg.lrr.nnd.plt, x = .55, y = 0,
            width = .45, height = .33) +
  annotate('text', x = .03, y = .55,
           label = 'Effect size of competition', size = 5,
           angle = 90) +
  annotate('text', x = .52, y = .97, label = 'A') +
  annotate('text', x = .97, y = .97, label = 'B') +
  annotate('text', x = .52, y = .65, label = 'C') +
  annotate('text', x = .97, y = .65, label = 'D') +
  annotate('text', x = .52, y = .3, label = 'E') +
  annotate('text', x = .97, y = .3, label = 'F') +
  annotate('text', x = .10, y = .785, label = 'Local',
           size = 5) +
  annotate('text', x = .10, y = .46, label = 'Local',
           size = 5) + 
  annotate('text', x = .09, y = .115, label = 'Regional',
           size = 5) + 
  annotate('text', x = .075, y = .97, label = 'MEPP Invasive',
           size = 3.5) +
  annotate('point', x = .013, y = .97, color = 'blue', alpha = 0.4) +
  annotate('text', x = .071, y = .95, label = 'MEPP Exotic',
           size = 3.5) +
  annotate('point', x = .013, y = .95, color = 'red', alpha = 0.4)
  
# functional phylogenetic regressions-------------
# next, determine best models using latest version of invasives FPD
# since abundances don't seem to matter much, I'll leave those out (for
# sake of computing time) and use the shiny app to find R2~a curve with highest
# R2

# ----- leaving RStudio Beep Boop, report back w/ best models-----

# ----- Beep Boop and we're back. Best model was Phylogeny + SLA, Ht, Flower period,
# and leaf toughness with an a-value of .3-.35 (results vary due to rarefying, 
# but the max R^2 is always ~.85). Creating NMDS plots for all species w/ those
# traits and arranging them in order of descending effect size of competition

# NMDS Plots for traits + phylogeny from best model ---------------
demo.data <- arrange(demo.data, desc(ESCR2))

traits <- c('Height', 'WoodDens', 'GrowthForm')
par(mfrow = c(2,2))
trait.data <- tyson$traits

# for(x in unique(demo.data$Species)){
#   
#   cat('Plotting species: ', x, '\n')
#   # make functional and phylogenetic distance matrices
#   phyloMat <- make_local_phylo_dist(x, communities, phylo)
#   funMat <- make_local_trait_dist(x, communities, 
#                                   trait.data = trait.data,
#                                   traits = traits,
#                                   scale = 'scaledBYrange') %>% as.matrix %>%
#             data.frame()
#   
#   # subset phylomat so it matches functional mat
#   phyloMat <- phyloMat[rownames(funMat),
#                        names(funMat)] %>% data.frame()
#   
#   # a little bit of defensive programming to ensure we don't actually
#   # mash up mismatched matrices
#   if(!identical(rownames(funMat), rownames(phyloMat)) |
#      !identical(names(funMat), names(phyloMat))){
#     stop('check code for subsampling phylo matrix')
#   }
#   # since the best model seems to depend on rarefied outputs, I'm
#   # using both values that resulted in best models. There shouldn't
#   # be any major differences though
#   UWDat35 <- func_phy_dist(FDist = funMat, PDist = phyloMat,
#                          phyloWeight = .35, p = 2)  
#   
#   UWDat3 <- func_phy_dist(FDist = funMat, PDist = phyloMat,
#                            phyloWeight = .3, p = 2) 
#   
#   # run NMDS and plot results for a = .35 and a = .3
#   # TraitNMDS <- metaMDS(UWDat35, trymax = 200, trace = 0)
#   # plot(TraitNMDS, type = 'n', main = paste0(x, ' a = .35'))
#   # orditorp(TraitNMDS, display = 'sites', air = .01,
#   #          select = !rownames(UWDat35) %in% x,
#   #          col='red')
#   #   text(TraitNMDS$points[rownames(TraitNMDS$points) == x,1],
#   #      TraitNMDS$points[rownames(TraitNMDS$points) == x,2],
#   #      label = x, col = 'blue', cex = .7)
#   #   
#   #   TraitNMDS <- metaMDS(UWDat3, trymax = 200, trace = 0)
#   #   plot(TraitNMDS, type = 'n', main = paste0(x, ' a = .3'))
#   #   orditorp(TraitNMDS, display = 'sites', air = .01,
#   #            select = !rownames(UWDat3) %in% x,
#   #            col='red')
#   #   text(TraitNMDS$points[rownames(TraitNMDS$points) == x,1],
#   #        TraitNMDS$points[rownames(TraitNMDS$points) == x,2],
#   #        label = x, col = 'blue', cex = .7)
# }
# 

# R^2~a value for best models from Shiny Phylo Fun App

a_seq <- seq(0, 1, .025)
R2dat <- data.frame(A = a_seq,
                    NND = rep(NA, length(a_seq)),
                    MPD = rep(NA, length(a_seq)))
mod.data <- demo.data

mod.data[ , paste0('nna_', a_seq)] <- NA
mod.data[ , paste0('mpa_', a_seq)] <- NA  

for(x in unique(demo.data$Species)){
  cat('Calculating FPD for species: ', x, '\n')
  
  
  # phylo and functional distance matrices
  phylo.mat <- make_local_phylo_dist(x, communities, phylo)
  fun.mat <- make_local_trait_dist(x, communities, trait.data,
                                   traits = traits,
                                   scale = 'scaledBYrange')
  
  # run for each level of a
  for(a in a_seq){
    # FPD matrix for rarefied communities
    FPD <- rarefy_FPD(x, phylo.mat = phylo.mat,
                      fun.mat = fun.mat,
                      n.rare = 11, a = a, p = 2,
                      abundance.weighted = FALSE,
                      community.data = NULL)
    
    # store results for each species
    mod.data[mod.data$Species == x, paste0('nna_', a)] <- FPD$rare.nnd
    mod.data[mod.data$Species == x, paste0('mpa_', a)] <- FPD$rare.mpd
    
  }
}

# if(any(is.na(mod.data))) { # alias for is.nan. Apparently doesn't work on data.frames
#   mod.data[is.na(mod.data)] <- NA
# }
# 
# # log transform data
# mod.data[, 14:95] <- log(mod.data[ ,14:95])

# run models for each level of a, extract R^2
for(a in a_seq){
  i <- which(a_seq == a)
  nnd.form <- as.formula(paste('ESCR2 ~ ', paste0('nna_', a), " + CRBM"))
  mpd.form <- as.formula(paste('ESCR2 ~ ', paste0('mpa_', a), " + CRBM"))
  
  R2dat$NND[i] <- r2_calc(mod.data, nnd.form)
  R2dat$MPD[i] <- r2_calc(mod.data, mpd.form) 
}

# find peak in curve
maxr2 <- max(R2dat[ ,2:3])
if(maxr2 %in% R2dat$MPD) {
  maxr2met <- 'MPD'
  maxr2A <- R2dat[which(R2dat$MPD == maxr2), 'A']
} else if(maxr2 %in% R2dat$NND) {
  maxr2met <- 'NND'
  maxr2A <- R2dat[which(R2dat$NND == maxr2), 'A']
  
}

# plot results
Fig <- ggplot(data = R2dat, aes(x = A)) +
  geom_point(aes(y = NND, color = 'NND'),
             alpha = 0.4,
             show.legend = FALSE) +
  geom_line(aes(y = NND, color = 'NND'), 
            alpha = 0.4,
            show.legend = FALSE) + 
  geom_point(aes(y = MPD, color = 'MPD'),
             alpha = 0.4,
             show.legend = FALSE) +
  geom_line(aes(y = MPD, color = 'MPD'), 
            alpha = 0.4,
            show.legend = FALSE) +
  scale_x_continuous('', limits = c(0,1)) +
  scale_y_continuous('', limits = c(0,1)) + 
  scale_color_manual('',
                     values = c("red", "blue")) + 
  annotate('text',
           label = paste0('Maximum R^2: ', round(maxr2, 3)),
           x = .85, y = .95) +
  annotate('text', 
           label = paste0('Best a-value: ', maxr2A),
           x = .85, y = .9) +
  annotate('text', 
           label = paste0('Best Performing Metric: ', maxr2met),
           x = .85, y = .85)

ggdraw() +
  draw_plot(Fig,
            x = .1, y = .1,
            height = 0.88, width = 0.88) +
  annotate('text', x = 0.075, y = 0.55,
           label = 'Adjusted R-Squared',
           angle = 90, size = 5) +
  annotate('text', x = 0.225, y = 0.09,
           label = 'Functional\n information \nonly') +
  annotate('text', x = 0.93, y = 0.09,
           label = 'Phylogenetic \ninformation \nonly') +
  annotate('text', x = 0.575, y = .11,
           label = 'Phylogenetic scaling parameter',
           size = 5) +
  annotate('text', x = 0.07, y = 0.95, 
           label = 'Metric', size = 4.5) +
  annotate('text', x = 0.077, y = 0.92,
           label = 'NND', size = 3.75) + 
  annotate('text', x = 0.078, y = 0.89,
           label = 'MPD', size = 3.75) +
  annotate('point', x = 0.04, y = 0.92,
           color = 'blue', alpha = 0.4) +
  annotate('point', x = 0.04, y = 0.89, 
          color = 'red', alpha = 0.4)


