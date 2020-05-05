# Script to produce all of the analyses and figures in the main text of the paper.

library(FunPhylo)
library(pez)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(tidyr)
library(stringr)
library(grid)
library(gridExtra)
library(png)
library(viridis)

data(tyson)

# do a little data munging and variable conversion

spp.list <- tyson$spp.list
for(i in 2:4) {
  spp.list[ ,i] <- as.numeric(spp.list[ ,i])
}

phylo <- tyson$phylo
communities <- tyson$communities
demo.data <- tyson$demo.data

# Get habitat info for every species in the community data set

temp <- demo.data %>%
  select(Species, Habitat) %>%
  rbind(
    data.frame(
      Species = c('Desmodium_perplexum',
                  'Geum_vernum',
                  'Teucrium_canadense',
                  'Symphoricarpos_orbiculatus'),
      Habitat = c('Grass', 'Forest', 'Forest', 'Forest')
    )
  )

com_hab <- vector('character',
                  length(communities$exotic_species))

for(i in seq_len(dim(communities)[1])) {
  
  foc_spp    <- communities$exotic_species[i]
  hab        <- temp$Habitat[temp$Species == foc_spp]
  com_hab[i] <- hab
  
}

communities$Habitat <- com_hab

# regional scale analysis-----------
# create Tyson scale phylo distance matrix
regionalTysonDists <- cophenetic(phylo) %>% sqrt()
regionalTysonDists <- data.frame(regionalTysonDists)


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

for(x in unique(exotics$Species)) {
  
  spp_hab <- exotics$Habitat[exotics$Species == x]
  
  spp_hab_regex <- gsub('; ', '|', spp_hab)
  
  spp_hab_ind   <- spp.list$Species[grepl(spp_hab_regex, spp.list$Habitat)]
  
  exotics[exotics$Species == x, 'MPD'] <- mean(regionalTysonDists[spp_hab_ind ,x],
                                               na.rm = TRUE)
  exotics[exotics$Species == x, 'NND'] <- min(regionalTysonDists[spp_hab_ind ,x],
                                              na.rm = TRUE)
  exotics[exotics$Species == x, 'Focal'] <- ifelse(x %in% demo.data$Species,
                                                   'Y', 
                                                   'N')
}


# local scale effect sizes, phylo only---------------
# phylogenetic only to start
# trait * phylogeny comes later

communities$MPD.Local <- NA_real_
communities$NND.Local <- NA_real_

All.Local.exotics <- subset(communities, exotic_species == '') %>%
  mutate(Focal = NA_real_)

demo.data$MEPPInv <- NA_real_
demo.data$MPD <- NA_real_
demo.data$NND <- NA_real_
demo.data$AWMPD <- NA_real_
demo.data$AWNND <- NA_real_
demo.data$Regional_MPD <- NA_real_
demo.data$Regional_NND <- NA_real_

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
# mpdLM <- lm(ESCR ~ MPD + CRBM, data = demo.data)
# nndLM <- lm(ESCR ~ NND + CRBM, data = demo.data)
# logAWmpdLM <- lm(ESCR ~ logAWMPD + CRBM, data = demo.data)
# logAWnndLM <- lm(ESCR ~ logAWNND + CRBM, data = demo.data)
# RegMPDLM <- lm(ESCR ~ Regional_MPD + CRBM, data = demo.data)
# RegNNDLM <- lm(ESCR ~ Regional_NND + CRBM, data = demo.data)

# Use raw effect sizes + variances as suggested by reviewer 1

demo.data$use_var <- demo.data$var_escr / sum(demo.data$var_escr)

mpdLM <- lm(ESCR ~ MPD + CRBM, 
            data = demo.data,
            weights = 1/use_var)

nndLM <- lm(ESCR ~ NND + CRBM, 
            data = demo.data,
            weights = 1 / use_var)

logAWmpdLM <- lm(ESCR ~ logAWMPD + CRBM, 
                 data = demo.data,
                 weights = 1 / use_var)

logAWnndLM <- lm(ESCR ~ logAWNND + CRBM, 
                 data = demo.data,
                 weights = 1 / use_var)

RegMPDLM <- lm(ESCR ~ Regional_MPD + CRBM, 
               data = demo.data,
               weights = 1 / use_var)

RegNNDLM <- lm(ESCR ~ Regional_NND + CRBM, 
               data = demo.data,
               weights = 1 / use_var)


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
write.csv(modelTable, '../Eco_Letters_Manuscript/Figures/Phylo_models_output.csv', 
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
                              y = ESCR)) +
  scale_x_continuous('MPD', 
                     breaks = seq(10,
                                  16,
                                  1.5),
                     limits = c(9.9, 16.1)) +
  scale_y_continuous('',
                     breaks = seq(0, 3.5, 1),
                     limits = c(-.5, 3.5)) +
  geom_point(alpha = .4, 
             aes(color = MEPPInv),
             size = 4) + 
  scale_color_manual(
    name = 'Status',
    labels = c("Non-invasive", "Invasive"),
    values = c('red','blue')
  ) +
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
                          aes(x = Magnitude, y = ESCR)) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE,
             size = 4) +
  scale_x_continuous('NND', 
                     breaks = seq(0,
                                  16,
                                  4),
                     limits = c(0, 16)) +
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
                             aes(x = Magnitude, y = ESCR)) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE,
             size = 4) +
  scale_x_continuous('log(abundance-weighted MPD)', 
                     breaks = seq(2.2,
                                  2.9,
                                  0.15),
                     limits = c(2.2, 2.9)) +
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
                             aes(x = Magnitude, y = ESCR)) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE,
             size = 4) +
  scale_x_continuous('log(abundance-weighted NND)', 
                     breaks = seq(-5,
                                  -0.5,
                                  0.75),
                     limits = c(-5, -0.5)) +
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
                          aes(x = Magnitude, y = ESCR)) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE,
             size = 4) +
  scale_x_continuous('MPD', 
                     breaks = seq(13.5,
                                  15.6,
                                  0.7),
                     limits = c(13.4, 15.7)) +
  scale_y_continuous('',
                     breaks = seq(0, 3.5, 1),
                     limits = c(-.5, 3.5)) +
  scale_color_manual(values = c('red','blue')) +
  geom_hline(yintercept = 0, linetype = 'dotted',
             alpha = 0.5,
             size = 2) + 
  plt.blank


reg.lrr.nnd.plt <- ggplot(data = filter(forPlot, Metric == 'Regional_NND'),
                          aes(x = Magnitude, y = ESCR)) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE,
             size = 4) +
  scale_x_continuous('NND', 
                     breaks = seq(1,
                                  11,
                                  2),
                     limits = c(1, 11)) +
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
  draw_plot(loc.lrr.mpd.plt + 
              theme(legend.position = 'top',
                    legend.direction = 'horizontal'),
            x = .1, y = .67,
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
  annotate('text', x = .52, y = .95,
           label = 'A', size = 5.5) +
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
           size = 5)

ggsave(filename = 'LRR_Regressions_Phylo_Only_For_Manuscript.png',
       path = '../Eco_Letters_Manuscript/Figures',
       height = 8.5,
       width = 12.5,
       units = 'in',
       dpi = 600)

ggsave(filename = "Figure_2.pdf",
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

demo.data <- arrange(demo.data, desc(ESCR))

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
  nnd.form <- as.formula(paste('ESCR ~ ', paste0('nna_', a), " + CRBM"))
  mpd.form <- as.formula(paste('ESCR ~ ', paste0('mpa_', a), " + CRBM"))
  
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
           x = .85, y = .85) +
  theme(
    panel.background = element_rect(fill = NA,
                                    color = 'black',
                                    size = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_text(size = 16,
                                    vjust = -0.1),
    axis.text.y      = element_text(size = 16,
                                    hjust = -0.1)
  )

ggdraw() +
  draw_plot(Fig,
            x = .1, y = .1,
            height = 0.88, width = 0.88) +
  annotate('text', x = 0.075, y = 0.55,
           label = 'Adjusted R-Squared',
           angle = 90, size = 5) +
  annotate('text', x = 0.19, y = 0.08,
           label = 'Functional\n information \nonly',
           size = 5) +
  annotate('text', x = 0.94, y = 0.08,
           label = 'Phylogenetic \ninformation \nonly',
           size = 5) +
  annotate('text', x = 0.56, y = .10,
           label = "paste('Phylogenetic Scaling Parameter ', italic(a), sep = '')",
           parse = TRUE,
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

ggsave(filename = 'R2_A_For_Manuscript.png',
       path = '../Eco_Letters_Manuscript/Figures',
       height = 8,
       width = 12.5,
       units = 'in',
       dpi = 600)

ggsave(filename = "Figure_3.pdf",
       path = '../Eco_Letters_Manuscript/Figures',
       height = 8.5,
       width = 12.5,
       units = 'in',
       dpi = 600)


## Conceptual figure - Figure 1!
detach('package:cowplot', unload = TRUE)

lrr_concept_data <- data.frame(
  Treatment   = c("Control", "Competitor Removal"),
  Significant = c(rep("paste('Within-species ', lambda, ' values\nnot significantly different')", 2),
                  rep("paste('Within-species ', lambda, ' values\nsignificantly different')",2)),
  lambda      = c(1.1, 1.2, 1.1, 1.6), 
  UpCI        = c(1.3, 1.5, 1.3, 1.9),
  LoCI        = c(0.92, 1.1, 0.92, 1.4)
)

lrr_panel <- ggplot(lrr_concept_data, 
                    aes(x = Treatment)) +
  geom_point(aes(y = lambda, 
                 color = Treatment),
             size = 4) + 
  geom_linerange(aes(ymax = UpCI,
                     ymin = LoCI,
                     color = Treatment),
                 size = 2) +
  geom_hline(yintercept = 1,
             alpha = 0.5,
             linetype = 'dotted',
             size = 2) +
  facet_wrap(~Significant) + 
  scale_y_continuous(expression(lambda)) +
  scale_color_manual(values = viridis::viridis(2, 
                                               begin = 0.3,
                                               end = 0.8)) + 
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.title.y = element_text(size = 35),
        axis.text.y  = element_text(size = 16),
        strip.background = element_blank(),
        strip.text.x     = element_blank(),
        panel.border     = element_rect(color = 'black',
                                        size = 2.5),
        legend.title     = element_text(size = 16),
        legend.text      = element_text(size = 14),
        plot.title       = element_text(size = 16)) + 
  ggtitle("B")


set.seed(1242152)
xx <- seq(1.1, 2.9, 0.3)
yy <- 1:2

sim_lrr_data <- data.frame(
  xx      = c(xx, 2.9, 2),
  sim_lrr = c(0.7 + -0.22 * xx + rnorm(7, 0, 0.1), 0, -0.1)
)


reg_concept_data <- data.frame(
  Novelty_Metric = yy,
  Text    = c("Competitive Interactions", 
              "Facilitative Interactions"),
  xmin_coord = c(1, 1),
  xmax_coord = c(3, 3),
  ymin_coord = c(0.03, -1),
  ymax_coord = c(1, -0.03),
  no_ef_line = 0
)


reg_panel <- ggplot(reg_concept_data,
                    aes(x = Novelty_Metric,
                        y = seq(-1, 1, length.out = 2)))  +
  geom_rect(aes(xmin = xmin_coord,
                xmax = xmax_coord,
                ymin = ymin_coord,
                ymax = ymax_coord,
                fill = Text),
            alpha = 0.6) + 
  geom_line(aes(y = no_ef_line,
                x = seq(1, 3, length.out = 2)),
            alpha = 0.9,
            size = 2,
            color = "red",
            linetype = 'dotted') +
  geom_point(data = sim_lrr_data,
             aes(x = xx,
                 y = sim_lrr),
             color = viridis::inferno(1, begin = 0.9, end = 1),
             size = 3) + 
  scale_y_continuous(
    expression(
      paste(
        frac(
          ln(lambda[CR] + 0.5), 
          ln(lambda[Control] + 0.5)
        ), 
        "    =   Effect Size of Competition"
      )
    )
  ) +
  scale_x_continuous("Distinctiveness Metric",
                     limits = c(1, 3),
                     breaks = c(1.2, 2.8),
                     labels = c("Less\nDistinct", "More\nDistinct")) + 
  scale_fill_manual("Interpretation of the direction of effect",
                    values = viridis::inferno(2,
                                              end = 0.3)) +
  theme(axis.title.y = element_text(size = 18,
                                    margin = margin(t = 0,
                                                    b = 0,
                                                    r = 10,
                                                    l = 5)),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16,
                                    margin = margin(t = 20,
                                                    b = 5,
                                                    l = 0,
                                                    r = 0)),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = NA,
                                        color = 'black',
                                        size = 2,
                                        linetype = 1),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        plot.title = element_text(hjust = 0,
                                  size = 16)
  ) + 
  ggtitle("C") + 
  guides(fill = guide_legend(title.position = 'top',
                              title.hjust = 0.5))


cr_concept_fig <- readPNG(
  '../Eco_Letters_Manuscript/Figures/Figure_1_Images/cr_concept_fig.png'
)

cr_use <- as.raster(cr_concept_fig)

cr_concept_grob <- rasterGrob(cr_use,
                              interpolate = TRUE)


cr_panel <- qplot(1:10, 1:10, geom = 'blank') + 
  annotation_custom(cr_concept_grob,
                    xmin=-Inf, 
                    xmax=Inf, 
                    ymin=-Inf, 
                    ymax=Inf) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.025,
                                  size = 16)) + 
  ggtitle("A")



pdf(file = "../Eco_Letters_Manuscript/Figures/Figure_1.pdf",
    width = 12,
    height = 9)
  grid.arrange(cr_panel,
               lrr_panel,
               reg_panel,
               layout_matrix = matrix(c(1, 3,
                                        2, 3), 
                                      byrow = TRUE,
                                      nrow = 2))
 
dev.off()
