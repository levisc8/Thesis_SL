# Re-run analyses with Smith Trees to see if any difference is made

# Script to produce all of the analyses and figures in the main text of the paper.

library(FunPhylo)
library(pez)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(tidyr)
library(stringr)

data(tyson)

# do a little data munging and variable conversion

spp.list <- tyson$spp.list
for(i in 2:4){
  spp.list[ ,i] <- as.numeric(spp.list[ ,i])
}

# remove allium vineale since it's a monocot. However, I wanted to keep
# it in the data set to be distributed because someone else might find it
# useful.

phylo <- tyson$phylo_all
communities <- tyson$communities %>% filter(community != 'Allium_vineale')
demo.data <- tyson$demo.data

# regional scale analysis-----------
# create Tyson scale phylo distance matrix
regionalTysonDists <- cophenetic(phylo) %>% data.frame()

diag(regionalTysonDists) <- NA

# next, extract data for all exotic species. Remove cerastium_spp. It is
# included in the species list so that it matches with the traits data, but it 
# is also probably covered by one of the three Cerastium species on the species 
# list.

exotics <- filter(spp.list, Exotic == 1 & Species != 'Cerastium_spp.')

exotics$MPD <- NA
exotics$NND <- NA
exotics$Focal <- NA

# extract mpd and nnd for all exotics. THis also creates a "Focal" dummy variable
# so that focal species can be colored for plotting

for(x in unique(exotics$Species)){
  exotics[exotics$Species == x, 'MPD'] <- mean(regionalTysonDists[ ,x],
                                               na.rm = TRUE)
  exotics[exotics$Species == x, 'NND'] <- min(regionalTysonDists[ ,x],
                                              na.rm = TRUE)
  exotics[exotics$Species == x, 'Focal'] <- ifelse(x %in% demo.data$Species,
                                                   'Y', 
                                                   'N')
}


# local scale effect sizes, phylo only---------------
# phylogenetic only to start
# trait * phylogeny comes later

communities$MPD.Local <- NA
communities$NND.Local <- NA

All.Local.exotics <- subset(communities, exotic_species == '') %>%
  mutate(Focal = NA)

demo.data$MEPPInv <- NA
demo.data$MPD <- NA
demo.data$NND <- NA
demo.data$AWMPD <- NA
demo.data$AWNND <- NA
demo.data$Regional_MPD <- NA
demo.data$Regional_NND <- NA

for(x in unique(demo.data$Species)) {
  cat('Crunching data for species: ', x, '\n')
  
  # make local phylo and functional distance matrices. The functional
  # distance matrices are really just place holders, they won't be used at
  # all because a = 1. I have to make them though because the function I wrote 
  # to do this requires a functional distance matrix (the reason for this will
  # become apparent further down when we start varying "a").
  
  demo.data[demo.data$Species == x, 'Regional_MPD'] <- exotics[exotics$Species == x,
                                                               'MPD']
  demo.data[demo.data$Species == x, 'Regional_NND'] <- exotics[exotics$Species == x,
                                                               'NND']
  
  local.exotics <- filter(communities, exotic_species == x &
                            alien == 1)
  
  phyloMat <- make_local_phylo_dist(x, communities, phylo)
  funMat <- make_local_trait_dist(x, communities, 
                                  trait.data = tyson$traits,
                                  traits = names(tyson$traits)[-1],
                                  scale = 'scaledBYrange')
  
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
  AWDat <- rarefy_FPD(x, 
                      phylo.mat = phyloMat,
                      fun.mat = funMat,
                      n.rare = 11,
                      a = 1,
                      p = 2, 
                      abundance.weighted = TRUE,
                      community.data = communities)
  # Unweighted data
  UWDat <- rarefy_FPD(x, 
                      phylo.mat = phyloMat,
                      fun.mat = funMat,
                      n.rare = 11,
                      a = 1,
                      p = 2, 
                      abundance.weighted = FALSE,
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

# more table formatting
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


modelTable$Parameter <- str_replace_all(modelTable$Parameter,
                                        '[:punct:]|[:digit:]', "")

modelTable[ ,5:9] <- round(modelTable[ ,5:9], 4)

add_stars <- function(x) {
  
  if(is.na(x)) return(NA_character_)
  
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
write.csv(modelTable, '../Eco_Letters_Manuscript/Figures/Phylo_models_output_all_OTB.csv', 
          na = "",
          row.names = FALSE)


# Figure 1 plots. THese are assembled on a single canvas using cowplot to add 
# annotations

bad.mets <- c('AWNND','AWMPD')
forPlot <- gather(demo.data, Metric, Magnitude, MPD:logAWNND) %>% 
  filter(!Metric %in% bad.mets )
forPlot$MEPPInv <- ifelse(forPlot$MEPPInv == 1, 'Y', 'N')

# generate predictions from each significant model so we can plot them as
# trendlines. Tragically, the stat_smooth model formula interface won't work with
# additional covariates (or, rather, I can't figure out a way to make it work.
# It may actually be possible).

x <- seq(min(demo.data$MPD), max(demo.data$MPD), length.out = 14)
y <- seq(min(demo.data$CRBM), max(demo.data$CRBM), length.out = 14)
MPD.Preds <- predict(mpdLM, data.frame(MPD = x,
                                       CRBM = y),
                     type = 'response',
                     interval = 'confidence',
                     se.fit = TRUE)$fit %>%
  data.frame %>% 
  cbind(., x)

x <- seq(min(demo.data$NND), max(demo.data$NND), length.out = 14)
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

plt.blank <- theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(fill = NA,
                                                   size = 1.25,
                                                   color = 'black'),
                   axis.title = element_text(size = 18),
                   axis.line = element_line(size = 1.3),
                   axis.text = element_text(size = 16))

loc.lrr.mpd.plt <- ggplot(data = filter(forPlot, Metric == 'MPD'),
                          aes(x = Magnitude,
                              y = ESCR2)) +
  scale_x_continuous('MPD', 
                     breaks = seq(160,
                                  260,
                                  20),
                     limits = c(145, 260)) +
  scale_y_continuous('',
                     breaks = seq(0, 3.5, 1),
                     limits = c(-.5, 3.5)) +
  geom_point(alpha = .4, 
             aes(color = MEPPInv),
             show.legend = FALSE,
             size = 4) + 
  scale_color_manual(values = c('red','blue')) +
  geom_hline(yintercept = 0, linetype = 'dotted',
             alpha = 0.5,
             size = 2) +
  geom_line(data = MPD.Preds, 
            aes(x = x, 
                y = fit),
            color = 'black',
            alpha = .8,
            size = .8)  + 
  plt.blank

loc.lrr.nnd.plt <- ggplot(data = filter(forPlot, Metric == 'NND'),
                          aes(x = Magnitude, y = ESCR2)) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE,
             size = 4) +
  scale_x_continuous('NND', 
                     breaks = seq(0,
                                  240,
                                  60),
                     limits = c(0, 240)) +
  scale_y_continuous('',
                     breaks = seq(0, 3.5, 1),
                     limits = c(-0.5, 3.5)) +
  scale_color_manual(values = c('red','blue')) +
  geom_hline(yintercept = 0, linetype = 'dotted',
             alpha = 0.5,
             size = 2) +
  geom_line(data = NND.Preds, 
            aes(x = x, 
                y = fit),
            color = 'black',
            alpha = .8,
            size = .8) +
  plt.blank


loc.lrr.aw.mpd.plt <- ggplot(data = filter(forPlot, Metric == 'logAWMPD'),
                             aes(x = Magnitude, y = ESCR2)) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE,
             size = 4) +
  scale_x_continuous('log(abundance-weighted MPD)', 
                     breaks = seq(4.8,
                                  5.7,
                                  .20),
                     limits = c(4.7, 5.7)) +
  scale_y_continuous('',
                     breaks = seq(0, 3.5, 1),
                     limits = c(-.5, 3.5)) +
  scale_color_manual(values = c('red','blue')) +
  geom_hline(yintercept = 0, linetype = 'dotted',
             alpha = 0.5,
             size = 2) +
  geom_line(data = AWMPD.Pred, 
            aes(x = x, 
                y = fit),
            color = 'black',
            alpha = .8,
            size = .8) + 
  plt.blank


loc.lrr.aw.nnd.plt <- ggplot(data = filter(forPlot, Metric == 'logAWNND'),
                             aes(x = Magnitude, y = ESCR2)) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE,
             size = 4) +
  scale_x_continuous('log(abundance-weighted NND)', 
                     breaks = seq(-3,
                                  2,
                                  .8),
                     limits = c(-3.15, 2)) +
  scale_y_continuous('',
                     breaks = seq(0, 3.5, 1),
                     limits = c(-0.75, 3.5)) +
  scale_color_manual(values = c('red','blue')) +
  geom_hline(yintercept = 0, linetype = 'dotted',
             alpha = 0.5,
             size = 2) +
  geom_line(data = AWNND.Pred, 
            aes(x = x, 
                y = fit),
            color = 'black',
            alpha = .8,
            size = .8) + 
  plt.blank


reg.lrr.mpd.plt <- ggplot(data = filter(forPlot, Metric == 'Regional_MPD'),
                          aes(x = Magnitude, y = ESCR2)) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE,
             size = 4) +
  scale_x_continuous('MPD', 
                     breaks = seq(200,
                                  250,
                                  10),
                     limits = c(200, 250)) +
  scale_y_continuous('',
                     breaks = seq(0, 3.5, 1),
                     limits = c(-.5, 3.5)) +
  scale_color_manual(values = c('red','blue')) +
  geom_hline(yintercept = 0, linetype = 'dotted',
             alpha = 0.5,
             size = 2) + 
  plt.blank


reg.lrr.nnd.plt <- ggplot(data = filter(forPlot, Metric == 'Regional_NND'),
                          aes(x = Magnitude, y = ESCR2)) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE,
             size = 4) +
  scale_x_continuous('NND', 
                     breaks = seq(0,
                                  140,
                                  20),
                     limits = c(0, 150)) +
  scale_y_continuous('',
                     breaks = seq(0, 3.5, 1),
                     limits = c(-.5, 3.5)) +
  scale_color_manual(values = c('red','blue')) +
  geom_hline(yintercept = 0, 
             linetype = 'dotted',
             alpha = 0.5,
             size = 2) +
  plt.blank

# Draw the plots and add axis labels where needed

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
           label = 'Effect size of competition', size = 6,
           angle = 90) +
  annotate('text', x = .52, y = .97,
           label = 'A', size = 5) +
  annotate('text', x = .97, y = .97,
           label = 'B', size = 5) +
  annotate('text', x = .52, y = .65, 
           label = 'C', size = 5) +
  annotate('text', x = .97, y = .65, 
           label = 'D', size = 5) +
  annotate('text', x = .52, y = .3, 
           label = 'E', size = 5) +
  annotate('text', x = .97, y = .3,
           label = 'F', size = 5) +
  annotate('text', x = .10, y = .785, 
           label = 'Local',
           size = 5) +
  annotate('text', x = .10, y = .46, 
           label = 'Local',
           size = 5) + 
  annotate('text', x = .09, y = .115, 
           label = 'Regional',
           size = 5) + 
  annotate('text', x = .0675, y = .97, 
           label = 'Invasive',
           size = 5) +
  annotate('point', x = .013, y = .97, color = 'blue', 
           alpha = 0.4, size = 4) +
  annotate('text', x = .06, y = .95, 
           label = 'Exotic',
           size = 5) +
  annotate('point', x = .013, y = .95,
           color = 'red', alpha = 0.4,
           size = 4)

ggsave(filename = 'LRR_Regressions_Phylo_Only_all_OTB_For_Manuscript.png',
       path = '../Eco_Letters_Manuscript/Figures',
       height = 8.5,
       width = 12.5,
       units = 'in',
       dpi = 600)

ggsave(filename = "Figure_1_all_OTB.pdf",
       path = '../Eco_Letters_Manuscript/Figures',
       height = 8.5,
       width = 12.5,
       units = 'in',
       dpi = 600)

# functional phylogenetic regressions-------------
# next, determine best models using latest version of invasives FPD
# https://sam-levin.shinyapps.io/Invasives_FPD/

# ----- leaving RStudio Beep Boop, report back w/ best models-----

# ----- Beep Boop and we're back. Best model was Phylogeny + SLA, Ht, Flower period,
# and leaf toughness with an a-value of 0.3-0.4 (best a's vary due to rarefying, 
# but the max R^2 is always ~0.8 and it's always in this range). 

demo.data <- arrange(demo.data, desc(ESCR2))

traits <- c('Height', 'SLA', 'Tough', 'Flower.Period')
trait.data <- tyson$traits

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
  phylo.mat <- make_local_phylo_dist(x, 
                                     communities, 
                                     phylo)
  fun.mat <- make_local_trait_dist(x, 
                                   communities, 
                                   trait.data,
                                   traits = traits,
                                   scale = 'scaledBYrange')
  
  # run for each level of a
  for(a in a_seq){
    # FPD matrix for rarefied communities
    FPD <- rarefy_FPD(x, phylo.mat = phylo.mat,
                      fun.mat = fun.mat,
                      n.rare = 11,
                      a = a, 
                      p = 2,
                      abundance.weighted = FALSE,
                      community.data = NULL)
    
    # store results for each species
    mod.data[mod.data$Species == x, paste0('nna_', a)] <- FPD$rare.nnd
    mod.data[mod.data$Species == x, paste0('mpa_', a)] <- FPD$rare.mpd
    
  }
}

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
             show.legend = FALSE,
             size = 3) +
  geom_line(aes(y = NND, color = 'NND'), 
            alpha = 0.4,
            show.legend = FALSE,
            size = 1.25) + 
  geom_point(aes(y = MPD, color = 'MPD'),
             alpha = 0.4,
             show.legend = FALSE,
             size = 3) +
  geom_line(aes(y = MPD, color = 'MPD'), 
            alpha = 0.4,
            show.legend = FALSE,
            size = 1.25) +
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
  annotate('text', x = 0.205, y = 0.09,
           label = 'Functional\n information \nonly',
           size = 5) +
  annotate('text', x = 0.93, y = 0.09,
           label = 'Phylogenetic \ninformation \nonly',
           size = 5) +
  annotate('text', x = 0.56, y = .11,
           label = 'Phylogenetic Scaling Parameter',
           size = 5) + 
  annotate('text', x = 0.69, y = .1075,
           label = '~italic(a)',
           parse = TRUE,
           size = 6) +
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

ggsave(filename = 'R2_A_all_OTB_For_Manuscript.png',
       path = '../Eco_Letters_Manuscript/Figures',
       height = 8,
       width = 12.5,
       units = 'in',
       dpi = 600)

ggsave(filename = "Figure_2_all_OTB.pdf",
       path = '../Eco_Letters_Manuscript/Figures',
       height = 8.5,
       width = 12.5,
       units = 'in',
       dpi = 600)
