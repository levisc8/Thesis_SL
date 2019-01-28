# Test for relationship between richness and ESCR
rm(list = ls())
library(FunPhylo)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(viridis)
library(ape)

data(tyson)

communities <- tyson$communities
demo.data <- tyson$demo.data
phylo <- tyson$phylo

Richness <- communities %>% group_by(exotic_species) %>%
            summarise(Richness = n()) %>% 
  filter(exotic_species %in% demo.data$Species) %>%
  cbind(demo.data, .)

Richness$Focal.Abundance <- NA
Richness$Other.Abundance <- NA
Richness$NN.Abundance <- NA

for(x in unique(Richness$Species)){
  focal.abund <- communities[communities$exotic_species == x &
                               communities$community == x, 'percentcover']
  
  other.abund <- filter(communities, exotic_species == x & 
                          community != x) %>%
    select(percentcover) %>% sum(na.rm = TRUE)

  com <- filter(communities, exotic_species == x & 
                 !is.na(percentcover)) %>% select(community) %>%
    unlist()
  
  phy <- drop.tip(phylo, setdiff(phylo$tip.label, com)) %>%
    cophenetic() %>%
    data.frame()
  
  diag(phy) <- NA
  nn <- which.min(phy[ ,x]) %>% 
    rownames(phy)[.]
  
  NN.abund <- filter(communities, exotic_species == x &
                       community == nn) %>% 
    select(percentcover) %>%
    unlist()
  
    
  Richness[Richness$Species == x, 'Focal.Abundance'] <- focal.abund
  Richness[Richness$Species == x, 'Other.Abundance'] <- other.abund
  Richness[Richness$Species == x, 'NN.Abundance'] <- NN.abund
  
}

# see whether log or normal richness is better for model
# par(mfrow = c(2,1))
# hist(Richness$Richness, breaks = 8)
# hist(log(Richness$Richness), breaks = 8)

# not much difference, going with normal richness
RichLM <- lm(ESCR2 ~ Richness, data = Richness)
summary(RichLM)

def_plot <- theme(plot.background = element_rect(fill = NA),
                  panel.background = element_rect(fill = NA),
                  axis.line = element_line(size = 1.25),
                  axis.text = element_text(size = 14),
                  axis.title = element_text(size = 16),
                  plot.margin = margin(20, 20, 20, 20),
                  legend.title = element_text(size = 16),
                  legend.text = element_text(size = 14),
                  legend.key = element_rect(fill = NA))

plt <- ggplot(Richness, aes(x = Richness, y = ESCR2)) + 
  def_plot +
  geom_point(aes(color = Habitat),
             size = 3.5,
             alpha = 0.7) + 
  ylab('Effect Size of Competition') +
  xlab('Site-level Species Richness')+
  scale_color_manual(values = c('black','orange','green'))
plt

ggsave('ESCR_by_Richness_Appendix.png',
       path = '../Eco_Letters_Manuscript/Figures',
       height = 7,
       width = 8,
       units = 'in',
       dpi = 600)

RichCRBM <- lm(CRBM ~ Richness, data = Richness)
summary(RichCRBM)

plt2 <- ggplot(Richness, aes(x = Richness, y = CRBM)) + 
  def_plot +  
  geom_point(aes(color = Habitat),
             size = 3.5,
             alpha = 0.7) + 
  ylab('Standardized Competitor Biomass Removed') +
  xlab('Site-level Species Richness') +
  scale_color_manual(values = c('black','orange','green')) +
  stat_smooth(method = 'lm', se = FALSE)
plt2

ggsave('CRBM_by_Richness_Appendix.png',
       path = '../Eco_Letters_Manuscript/Figures',
       height = 7,
       width = 8,
       units = 'in',
       dpi = 600)

RichAbund <- lm(Richness ~ Focal.Abundance, data = Richness)
summary(RichAbund)

plt3 <- ggplot(Richness, aes(x = Focal.Abundance, y = Richness)) +
  def_plot + 
  geom_point(aes(color = Habitat),
             size = 3.5,
             alpha = 0.7) + 
  ylab('Site-level Species Richness') +
  xlab('Abundance of Focal Species') +
  scale_color_manual(values = c('black','orange','green')) 

plt3

ggsave('Richness_Focal_Abundance_Appendix.png',
       path = '../Eco_Letters_Manuscript/Figures',
       height = 7,
       width = 8,
       units = 'in',
       dpi = 600)

AbundLM <- lm(Richness$Other.Abundance ~ log(Richness$Focal.Abundance))
summary(AbundLM)

# plt4 <- ggplot(Richness, aes(y = Other.Abundance,
  #                            x = log(Focal.Abundance))) +
  # theme_tufte() + 
  # geom_rangeframe() +
  # geom_point(aes(color = Habitat)) + 
  # xlab('Abundance of Focal Species') +
  # ylab('Abundance of Co-occurring Species') +
  # scale_color_manual(values = c('black','orange','green')) +
  # stat_smooth(method = 'lm', se = F)

# plt4

plt5 <- plt4 <- ggplot(Richness, aes(y = NN.Abundance,
                                     x = ESCR2)) +
  def_plot + 
  geom_point(aes(color = Habitat),
             size = 3,
             alpha = 0.7) + 
  xlab('Effect size of Competitor Removal') +
  ylab('Abundance of Phylogenetic Nearest Neighbor') +
  scale_color_manual(values = c('black','orange','green')) 

plt5

ggsave('NN_Abundance_ESCR_Appendix.png',
       path = '../Eco_Letters_Manuscript/Figures',
       height = 7,
       width = 8,
       units = 'in',
       dpi = 600)


outTable <- select(Richness, -c(exotic_species, Other.Abundance))

outTable$CRBM <- exp(outTable$CRBM)
outTable$Species <- gsub("_", " ", outTable$Species)

write.csv(outTable,
          '../Eco_Letters_Manuscript/Figures/Abundance_Richness_Table.csv')

# different tests of ESCR~Biomass ------------

load('C:/Users/sl13sise/Desktop/Tyson/Data and design/Biomass/R/CRBMMeans.rda')
Means[Means$Site2 == 'privet', 'Site2'] <- 'ligustrum'
Means[Means$Site2 == 'microthlaspi', 'Site2'] <- 'thlaspi'


for(x in unique(demo.data$Species)){
  y <- tolower(stringr::str_split(x, "_")[[1]][1])
  demo.data$RawCRBM[demo.data$Species == x] <- Means[Means$Site2 == y,
                                                     'FirstCut'] 
}

demo.data$RawCRBM <- unlist(demo.data$RawCRBM)
# get phylogeny and communities ready
phylo <- tyson$phylo
communities <- tyson$communities

#
# create Tyson scale phylo distance matrix
regionalTysonDists <- cophenetic(phylo) %>% data.frame()

diag(regionalTysonDists) <- NA

demo.data$MEPPInv <- NA
demo.data$MPD <- NA
demo.data$NND <- NA
demo.data$AWMPD <- NA
demo.data$AWNND <- NA
demo.data$Regional_MPD <- NA
demo.data$Regional_NND <- NA

for(x in unique(demo.data$Species)) {
  cat('Crunching data for species: ', x, '\n')
  # extract mpd and nnd for focal species with regional species pool
  demo.data[demo.data$Species == x, 'Regional_MPD'] <- mean(regionalTysonDists[ ,x],
                                                            na.rm = TRUE)
  demo.data[demo.data$Species == x, 'Regional_NND'] <- min(regionalTysonDists[ ,x],
                                                           na.rm = TRUE)
  
  # make local phylo and functional distance matrices. The functional
  # distance matrices are really just place holders, they won't be used at
  # all because a = 1. I have to make them though because the function I wrote 
  # to do this requires a functional distance matrix too because I'm an idiot.
  phyloMat <- make_local_phylo_dist(x, communities, phylo)
  funMat <- make_local_trait_dist(x, communities, 
                                  trait.data = tyson$traits,
                                  traits = names(tyson$traits)[-1],
                                  scale = 'scaledBYrange')
  
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

RawCrbmLM <- lm(ESCR2 ~ RawCRBM, data = demo.data)
summary(RawCrbmLM)

# unstandardized biomass
RAWmpdLM <- lm(ESCR2 ~ MPD + RawCRBM, data = demo.data)
RAWnndLM <- lm(ESCR2 ~ NND + RawCRBM, data = demo.data)
RAWlogAWmpdLM <- lm(ESCR2 ~ logAWMPD + RawCRBM, data = demo.data)
RAWlogAWnndLM <- lm(ESCR2 ~ logAWNND + RawCRBM, data = demo.data)
RAWRegMPDLM <- lm(ESCR2 ~ Regional_MPD + RawCRBM, data = demo.data)
RAWRegNNDLM <- lm(ESCR2 ~ Regional_NND + RawCRBM, data = demo.data)

message('unstandardized biomass')
summary(RAWmpdLM)
summary(RAWnndLM)
summary(RAWlogAWmpdLM)
summary(RAWlogAWnndLM)
summary(RAWRegMPDLM)
summary(RAWRegNNDLM)

# without Biomass

NBMmpdLM <- lm(ESCR2 ~ MPD, data = demo.data)
NBMnndLM <- lm(ESCR2 ~ NND, data = demo.data)
NBMlogAWmpdLM <- lm(ESCR2 ~ logAWMPD, data = demo.data)
NBMlogAWnndLM <- lm(ESCR2 ~ logAWNND, data = demo.data)
NBMRegMPDLM <- lm(ESCR2 ~ Regional_MPD, data = demo.data)
NBMRegNNDLM <- lm(ESCR2 ~ Regional_NND, data = demo.data)

message('without biomass')

summary(NBMmpdLM)
summary(NBMnndLM)
summary(NBMlogAWmpdLM)
summary(NBMlogAWnndLM)
summary(NBMRegMPDLM)
summary(NBMRegNNDLM)

# Now with residuals from asymmetric covariate
demo.data$AsymResids <- residuals(lm(ESCR2 ~ CRBM, data = demo.data))

RESmpdLM <- lm(AsymResids ~ MPD, data = demo.data)
RESnndLM <- lm(AsymResids ~ NND, data = demo.data)
RESlogAWmpdLM <- lm(AsymResids ~ logAWMPD, data = demo.data)
RESlogAWnndLM <- lm(AsymResids ~ logAWNND, data = demo.data)
RESRegMPDLM <- lm(AsymResids ~ Regional_MPD, data = demo.data)
RESRegNNDLM <- lm(AsymResids ~ Regional_NND, data = demo.data)

message('residuals from standard biomass')

summary(RESmpdLM)
summary(RESnndLM)
summary(RESlogAWmpdLM)
summary(RESlogAWnndLM)
summary(RESRegMPDLM)
summary(RESRegNNDLM)

demo.data$MEPPInv <- ifelse(demo.data$MEPPInv == 1,
                            "Invasive",
                            "Exotic")
# Plot ESCR2 ~ RawCRBM
Fig6 <- ggplot(demo.data, aes(x = RawCRBM, y = ESCR2)) +
  def_plot + 
  geom_point(aes(color = MEPPInv),
             size = 3,
             alpha = 0.4) +
  scale_color_manual('',
                     breaks = c('Invasive', 'Exotic'),
                     values = c('red','blue')) +
  scale_x_continuous('Raw competitor biomass') +
  scale_y_continuous('Effect size of competition') + 
  stat_smooth(method = 'lm', 
              formula = y ~ x,
              se = FALSE,
              color = 'grey',
              alpha = 0.4)

Fig6

# Make Figures S3.6 (resids(lm(ESCR2~CRBM)) ~ novelty) and
# S3.7 (ESCR2 ~ novelty, no biomass)
bad.mets <- c('AWMPD', 'AWNND')
library(tidyr)

forPlot <- gather(demo.data, Metric, Magnitude, MPD:logAWNND) %>% 
  filter(!Metric %in% bad.mets )
# forPlot$MEPPInv <- ifelse(forPlot$MEPPInv == 1, 'Y', 'N')


loc.lrr.mpd.plt <- ggplot(data = filter(forPlot, Metric == 'MPD'),
                          aes(x = Magnitude,
                              y = AsymResids)) +
  stat_smooth(formula = y ~ x,
              method = 'lm',
              se = FALSE,
              color = 'grey',
              alpha = .3) + 
  scale_x_continuous('MPD', 
                     breaks = seq(140,
                                  260,
                                  40),
                     limits = c(135, 260)) +
  scale_y_continuous('',
                     breaks = seq(-2, 2, 1),
                     limits = c(-2, 2)) +
  geom_point(alpha = .4, aes(color = MEPPInv),
             show.legend = FALSE) +
  scale_color_manual(values = c('red','blue'))

loc.lrr.nnd.plt <- ggplot(data = filter(forPlot, Metric == 'NND'),
                          aes(x = Magnitude, y = AsymResids)) +
  stat_smooth(formula = y ~ x,
              method = 'lm',
              se = FALSE,
              color = 'grey',
              alpha = .3) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE) +
  scale_x_continuous('NND', 
                     breaks = seq(0,
                                  240,
                                  60),
                     limits = c(0, 240)) +
  scale_y_continuous('',
                     breaks = seq(-2, 2, 1),
                     limits = c(-2, 2)) +
  scale_color_manual(values = c('red','blue')) 


loc.lrr.aw.mpd.plt <- ggplot(data = filter(forPlot, Metric == 'logAWMPD'),
                             aes(x = Magnitude, y = AsymResids)) +
  stat_smooth(formula = y ~ x,
              method = 'lm',
              se = FALSE,
              color = 'grey',
              alpha = .3) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE) +
  scale_x_continuous('log(abundance-weighted MPD)', 
                     breaks = seq(4.5,
                                  6,
                                  .3),
                     limits = c(4.5, 6)) +
  scale_y_continuous('',
                     breaks = seq(-2, 2, 1),
                     limits = c(-2, 2)) +
  scale_color_manual(values = c('red','blue')) 


loc.lrr.aw.nnd.plt <- ggplot(data = filter(forPlot, Metric == 'logAWNND'),
                             aes(x = Magnitude, y = AsymResids)) +
  stat_smooth(formula = y ~ x,
              method = 'lm',
              se = FALSE,
              color = 'grey',
              alpha = .3) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE) +
  scale_x_continuous('log(abundance-weighted NND)', 
                     breaks = seq(-3,
                                  2,
                                  .8),
                     limits = c(-3.1, 2)) +
  scale_y_continuous('',
                     breaks = seq(-2, 2, 1),
                     limits = c(-2, 2)) +
  scale_color_manual(values = c('red','blue')) 


reg.lrr.mpd.plt <- ggplot(data = filter(forPlot, Metric == 'Regional_MPD'),
                          aes(x = Magnitude, y = AsymResids)) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE) +
  scale_x_continuous('MPD', 
                     breaks = seq(200,
                                  250,
                                  10),
                     limits = c(200, 250)) +
  scale_y_continuous('',
                     breaks = seq(-2, 2, 1),
                     limits = c(-2, 2)) +
  scale_color_manual(values = c('red','blue')) 


reg.lrr.nnd.plt <- ggplot(data = filter(forPlot, Metric == 'Regional_NND'),
                          aes(x = Magnitude, y = AsymResids)) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE) +
  scale_x_continuous('NND', 
                     breaks = seq(0,
                                  100,
                                  20),
                     limits = c(0, 110)) +
  scale_y_continuous('',
                     breaks = seq(-2, 2, 1),
                     limits = c(-2, 2)) +
  scale_color_manual(values = c('red','blue')) 
# create blank canvas and get to drawing!

library(cowplot)
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
           label = 'Residuals from CR Biomass', size = 5,
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

ggsave('Residuals_by_Novelty_LRR_Appendix.png',
       path = '../Eco_Letters_Manuscript/Figures',
       height = 8.5,
       width = 12.5,
       units = 'in',
       dpi = 600)

# Test traits for ESCR ~ Trait relationship
trait.demo.data <- tyson$traits %>%
  setNames(c('Species', names(tyson$traits)[-1])) %>%
  filter(Species %in% demo.data$Species) %>%
  left_join(., demo.data, by = 'Species')

SLA.LM <- lm(ESCR2 ~ SLA*Habitat + CRBM, data = trait.demo.data)
summary(SLA.LM)

HT.LM <- lm(ESCR2 ~ Height*Habitat + CRBM, data = trait.demo.data)
summary(HT.LM)

Tough.LM <- lm(ESCR2 ~ Tough + CRBM, data = trait.demo.data)
summary(Tough.LM)
# Omit columns for which we only have a couple or have no observations
binary.columns <- c('Woody', 'Clonal', 'N_Fixer',
                    'Stemmed_Herb', 'Tree', 'Rosette',
                    'Shrub', 'Unassisted', 'Ballistic',
                    'EndoZoochory', 'Water')

for(i in unique(binary.columns)){
  message(paste0('\nAnova Results for trait: ', i))
  Anova <- aov(trait.demo.data$ESCR2 ~ unlist(trait.demo.data[ ,i]))
  print(summary(Anova))
}
