# Test for relationship between richness and ESCR
rm(list = ls())
library(FunPhylo)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(viridis)
library(ape)
library(glue)
library(purrr)

data(tyson)

communities <- tyson$communities
demo.data <- tyson$demo.data
phylo <- tyson$phylo
spp.list <- tyson$spp.list

Richness <- communities %>% group_by(exotic_species) %>%
  summarise(Richness = n()) %>% 
  filter(exotic_species %in% demo.data$Species) %>%
  cbind(demo.data, .)

Richness$Focal.Abundance <- NA
Richness$Other.Abundance <- NA
Richness$NN.Abundance <- NA
Richness$NND <- NA
Richness$MPD <- NA

for(x in unique(Richness$Species)){
  focal.abund <- communities[communities$exotic_species == x &
                               communities$community == x, 'percentcover']
  
  other.abund <- filter(communities, exotic_species == x & 
                          community != x) %>%
    select(percentcover) %>% 
    sum(na.rm = TRUE)
  
  com <- filter(communities, exotic_species == x & 
                  !is.na(percentcover)) %>% select(community) %>%
    unlist()
  
  phy <- drop.tip(phylo, setdiff(phylo$tip.label, com)) %>%
    cophenetic() %>%
    data.frame()
  
  diag(phy) <- NA
  nnd <- min(phy[ ,x], na.rm = TRUE)
  mpd <- mean(phy[ ,x], na.rm = TRUE)
  
  nn <- which.min(phy[ ,x]) %>% 
    rownames(phy)[.]
  
  NN.abund <- filter(communities, exotic_species == x &
                       community == nn) %>% 
    select(percentcover) %>%
    unlist()
  
  
  Richness[Richness$Species == x, 'Focal.Abundance'] <- focal.abund
  Richness[Richness$Species == x, 'Other.Abundance'] <- other.abund
  Richness[Richness$Species == x, 'NN.Abundance'] <- NN.abund
  Richness[Richness$Species == x, 'NND'] <- nnd
  Richness[Richness$Species == x, 'MPD'] <- mpd
}

# see whether log or normal richness is better for model
# par(mfrow = c(2,1))
# hist(Richness$Richness, breaks = 8)
# hist(log(Richness$Richness), breaks = 8)

# not much difference, going with normal richness
RichLM <- lm(ESCR ~ Richness, data = Richness)
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

plt <- ggplot(Richness, aes(x = Richness, y = ESCR)) + 
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

ggsave('Fig_S3_1.pdf',
       path = '../Eco_Letters_Manuscript/SI/Figures',
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

ggsave('Fig_S3_2.pdf',
       path = '../Eco_Letters_Manuscript/SI/Figures',
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

ggsave('Fig_S3_3.pdf',
       path = '../Eco_Letters_Manuscript/SI/Figures',
       height = 7,
       width = 8,
       units = 'in',
       dpi = 600)


AbundLM <- lm(Richness$Other.Abundance ~ Richness$Focal.Abundance)
summary(AbundLM)

plt4 <- ggplot(Richness, aes(y = Other.Abundance,
                             x = Focal.Abundance)) +
  def_plot + 
  geom_point(aes(color = Habitat),
             size = 3,
             alpha = 0.7) +
  xlab('Abundance of Focal Species') +
  ylab('Abundance of Co-occurring Species') +
  scale_color_manual(values = c('black','orange','green')) +
  stat_smooth(method = 'lm', se = FALSE, size = 1)

plt4

ggsave('All_Abundance_Focal_Abundance_Appendix.png',
       path = '../Eco_Letters_Manuscript/Figures',
       height = 7,
       width = 8,
       units = 'in',
       dpi = 600)

ggsave('Fig_S3_4.pdf',
       path = '../Eco_Letters_Manuscript/SI/Figures',
       height = 7,
       width = 8,
       units = 'in',
       dpi = 600)


summary(lm(NN.Abundance ~ ESCR, data = Richness))


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
regionalTysonDists <- cophenetic(phylo) %>% sqrt()
regionalTysonDists <- data.frame(regionalTysonDists)

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
  
  spp_hab <- spp.list$Habitat[spp.list$Species == x]
  
  spp_hab_regex <- gsub('; ', '|', spp_hab)
  
  spp_hab_ind   <- spp.list$Species[grepl(spp_hab_regex, spp.list$Habitat)]
  
  demo.data[demo.data$Species == x, 'Regional_MPD'] <- mean(regionalTysonDists[spp_hab_ind , x],
                                                            na.rm = TRUE)
  demo.data[demo.data$Species == x, 'Regional_NND'] <- min(regionalTysonDists[spp_hab_ind , x],
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
                      p = 2, abundance.weighted = TRUE,
                      community.data = communities)
  # Unweighted data
  UWDat <- rarefy_FPD(x, phylo.mat = phyloMat,
                      fun.mat = funMat,
                      n.rare = 11, a = 1,
                      p = 2, abundance.weighted = FALSE,
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

RawCrbmLM <- lm(ESCR ~ RawCRBM, data = demo.data)
summary(RawCrbmLM)

# unstandardized biomass
RAWmpdLM <- lm(ESCR ~ MPD + RawCRBM, data = demo.data)
RAWnndLM <- lm(ESCR ~ NND + RawCRBM, data = demo.data)
RAWlogAWmpdLM <- lm(ESCR ~ logAWMPD + RawCRBM, data = demo.data)
RAWlogAWnndLM <- lm(ESCR ~ logAWNND + RawCRBM, data = demo.data)
RAWRegMPDLM <- lm(ESCR ~ Regional_MPD + RawCRBM, data = demo.data)
RAWRegNNDLM <- lm(ESCR ~ Regional_NND + RawCRBM, data = demo.data)

message('unstandardized biomass')
summary(RAWmpdLM)
summary(RAWnndLM)
summary(RAWlogAWmpdLM)
summary(RAWlogAWnndLM)
summary(RAWRegMPDLM)
summary(RAWRegNNDLM)

# without Biomass

NBMmpdLM <- lm(ESCR ~ MPD, data = demo.data)
NBMnndLM <- lm(ESCR ~ NND, data = demo.data)
NBMlogAWmpdLM <- lm(ESCR ~ logAWMPD, data = demo.data)
NBMlogAWnndLM <- lm(ESCR ~ logAWNND, data = demo.data)
NBMRegMPDLM <- lm(ESCR ~ Regional_MPD, data = demo.data)
NBMRegNNDLM <- lm(ESCR ~ Regional_NND, data = demo.data)

message('without biomass')

summary(NBMmpdLM)
summary(NBMnndLM)
summary(NBMlogAWmpdLM)
summary(NBMlogAWnndLM)
summary(NBMRegMPDLM)
summary(NBMRegNNDLM)

# Now with residuals from asymmetric covariate
demo.data$AsymResids <- residuals(lm(ESCR ~ CRBM, data = demo.data))

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
                            "Naturalized")
# Plot ESCR ~ RawCRBM
Fig6 <- ggplot(demo.data, aes(x = RawCRBM, y = ESCR)) +
  def_plot + 
  geom_point(aes(color = MEPPInv),
             size = 3,
             alpha = 0.4) +
  scale_color_manual('',
                     breaks = c('Invasive', 'Naturalized'),
                     values = c('red','blue')) +
  scale_x_continuous('Raw competitor biomass') +
  scale_y_continuous('Effect size of competition') + 
  stat_smooth(method = 'lm', 
              formula = y ~ x,
              se = FALSE,
              color = 'grey',
              alpha = 0.4,
              size = 1.5)

Fig6

ggsave('ESCR_RawCRBM_Appendix.png',
       path = '../Eco_Letters_Manuscript/Figures',
       height = 7,
       width = 8,
       units = 'in',
       dpi = 600)

ggsave('Fig_S3_5.pdf',
       path = '../Eco_Letters_Manuscript/SI/Figures',
       height = 7,
       width = 8,
       units = 'in',
       dpi = 600)


# Make Figures S3.6 (resids(lm(ESCR~CRBM)) ~ novelty) and
# S3.7 (ESCR ~ novelty, no biomass)
bad.mets <- c('AWMPD', 'AWNND')
library(tidyr)

forPlot <- gather(demo.data, Metric, Magnitude, MPD:logAWNND) %>% 
  filter(!Metric %in% bad.mets )
# forPlot$MEPPInv <- ifelse(forPlot$MEPPInv == 1, 'Y', 'N')


loc.lrr.mpd.plt <- ggplot(data = filter(forPlot, Metric == 'MPD'),
                          aes(x = Magnitude,
                              y = AsymResids)) + 
  geom_hline(yintercept = 0, linetype = 'dotted',
             alpha = 0.5,
             size = 2)  + 
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             # show.legend = FALSE,
             size = 4) +
  stat_smooth(formula = y ~ x,
              method = 'lm',
              se = FALSE,
              color = 'grey',
              alpha = .3) + 
  scale_x_continuous('MPD', 
                     breaks = seq(10,
                                  16,
                                  1.5),
                     limits = c(9.9, 16.1)) +
  scale_y_continuous('',
                     breaks = seq(-1, 2, 1),
                     limits = c(-1.1, 2)) +
  scale_color_manual(
    name = 'Status',
    labels = c("Naturalized", "Invasive"),
    values = c('red','blue')
  ) +
  geom_hline(yintercept = 0, linetype = 'dotted',
             alpha = 0.5,
             size = 2) 


loc.lrr.nnd.plt <- ggplot(data = filter(forPlot, Metric == 'NND'),
                          aes(x = Magnitude, y = AsymResids)) + 
  geom_hline(yintercept = 0, linetype = 'dotted',
             alpha = 0.5,
             size = 2) +
  stat_smooth(formula = y ~ x,
              method = 'lm',
              se = FALSE,
              color = 'grey',
              alpha = .3) +
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
                     breaks = seq(-1, 2, 1),
                     limits = c(-1.1, 2)) +
  scale_color_manual(values = c('red','blue'))


loc.lrr.aw.mpd.plt <- ggplot(data = filter(forPlot, Metric == 'logAWMPD'),
                             aes(x = Magnitude, y = AsymResids)) + 
  geom_hline(yintercept = 0, linetype = 'dotted',
             alpha = 0.5,
             size = 2) +
  stat_smooth(formula = y ~ x,
              method = 'lm',
              se = FALSE,
              color = 'grey',
              alpha = .3) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE, 
             size = 4) +
  scale_x_continuous('log(abundance-weighted MPD)', 
                     breaks = seq(2.2,
                                  2.9,
                                  .15),
                     limits = c(2.2, 2.9)) +
  scale_y_continuous('',
                     breaks = seq(-1, 2, 1),
                     limits = c(-1.1, 2.3)) +
  scale_color_manual(values = c('red','blue'))

loc.lrr.aw.nnd.plt <- ggplot(data = filter(forPlot, Metric == 'logAWNND'),
                             aes(x = Magnitude, y = AsymResids)) + 
  geom_hline(yintercept = 0, linetype = 'dotted',
             alpha = 0.5,
             size = 2)  +
  stat_smooth(formula = y ~ x,
              method = 'lm',
              se = FALSE,
              color = 'grey',
              alpha = .3) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE,
             size = 4) +
  scale_x_continuous('log(abundance-weighted NND)', 
                     breaks = seq(-5,
                                  -0.5,
                                  .8),
                     limits = c(-5, -0.5)) +
  scale_y_continuous('',
                     breaks = seq(-1, 2, 1),
                     limits = c(-1.1, 2.3)) +
  scale_color_manual(values = c('red','blue'))



reg.lrr.mpd.plt <- ggplot(data = filter(forPlot, Metric == 'Regional_MPD'),
                          aes(x = Magnitude, y = AsymResids)) + 
  geom_hline(yintercept = 0, linetype = 'dotted',
             alpha = 0.5,
             size = 2)  +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE,
             size = 4) +
  scale_x_continuous('MPD', 
                     breaks = seq(13.5,
                                  15.6,
                                  0.7),
                     limits = c(13, 15.7)) +
  scale_y_continuous('',
                     breaks = seq(-1, 2, 1),
                     limits = c(-1.1, 2.3)) +
  scale_color_manual(values = c('red','blue'))


reg.lrr.nnd.plt <- ggplot(data = filter(forPlot, Metric == 'Regional_NND'),
                          aes(x = Magnitude, y = AsymResids)) + 
  geom_hline(yintercept = 0, linetype = 'dotted',
             alpha = 0.5,
             size = 2) +
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
                     breaks = seq(-1, 2, 1),
                     limits = c(-1.1, 2.3)) +
  scale_color_manual(values = c('red','blue')) 

# create blank canvas and get to drawing!

library(cowplot)
ggdraw() +
  draw_plot(loc.lrr.mpd.plt + 
              theme(legend.position = 'top',
                    legend.direction = 'horizontal'),
            x = .1, y = .67,
            width = .45, height = .33) + 
  draw_plot(loc.lrr.nnd.plt, 
            x = .55, y = .67,
            width = .45, height = .33) +
  draw_plot(loc.lrr.aw.mpd.plt, 
            x = .1, y = .34,
            width = .45, height = .33) +
  draw_plot(loc.lrr.aw.nnd.plt, 
            x = .55, y = .34,
            width = .45, height = .33) +
  draw_plot(reg.lrr.mpd.plt, 
            x = .1, y = 0,
            width = .45, height = .33) +
  draw_plot(reg.lrr.nnd.plt, 
            x = .55, y = 0,
            width = .45, height = .33) +
  annotate('text',
           x = .03, y = .55,
           label = 'Effect size of competition (Regression using residuals)',
           size = 6, angle = 90) +
  annotate('text', x = .52, y = .95,
           label = 'A', size = 5.5) +
  annotate('text', x = .97, y = .97, label = 'B',
           size = 5) +
  annotate('text', x = .52, y = .65, label = 'C',
           size = 5) +
  annotate('text', x = .97, y = .65, label = 'D',
           size = 5) +
  annotate('text', x = .52, y = .3, label = 'E',
           size = 5) +
  annotate('text', x = .97, y = .3, label = 'F',
           size = 5) +
  annotate('text', x = .10, y = .785, label = 'Local',
           size = 5) +
  annotate('text', x = .10, y = .46, label = 'Local',
           size = 5)

ggplot2::ggsave('Residuals_by_Novelty_LRR_Appendix.png',
                path = '../Eco_Letters_Manuscript/Figures',
                height = 8.5,
                width = 12.5,
                units = 'in',
                dpi = 600)

ggsave('Fig_S3_7.pdf',
       path = '../Eco_Letters_Manuscript/SI/Figures',
       height = 7,
       width = 8,
       units = 'in',
       dpi = 600)


# No Biomass in regressions

loc.lrr.mpd.plt <- ggplot(data = filter(forPlot, Metric == 'MPD'),
                          aes(x = Magnitude,
                              y = ESCR)) + 
  geom_hline(yintercept = 0, linetype = 'dotted',
             alpha = 0.5,
             size = 2)  + 
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             # show.legend = FALSE,
             size = 4) +
  stat_smooth(formula = y ~ x,
              method = 'lm',
              se = FALSE,
              color = 'grey',
              alpha = .3) + 
  scale_x_continuous('MPD', 
                     breaks = seq(10,
                                  16,
                                  1.5),
                     limits = c(9.9, 16.1)) +
  scale_y_continuous('',
                     breaks = seq(0, 3.5, 1),
                     limits = c(-0.5, 3.5)) +
  scale_color_manual(
    name = 'Status',
    labels = c("Naturalized", "Invasive"),
    values = c('red','blue')
  ) +
  geom_hline(yintercept = 0, linetype = 'dotted',
             alpha = 0.5,
             size = 2) 


loc.lrr.nnd.plt <- ggplot(data = filter(forPlot, Metric == 'NND'),
                          aes(x = Magnitude, y = ESCR)) + 
  geom_hline(yintercept = 0, linetype = 'dotted',
             alpha = 0.5,
             size = 2) +
  stat_smooth(formula = y ~ x,
              method = 'lm',
              se = FALSE,
              color = 'grey',
              alpha = .3) +
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
  scale_color_manual(values = c('red','blue'))


loc.lrr.aw.mpd.plt <- ggplot(data = filter(forPlot, Metric == 'logAWMPD'),
                             aes(x = Magnitude, y = ESCR)) + 
  geom_hline(yintercept = 0, linetype = 'dotted',
             alpha = 0.5,
             size = 2) +
  stat_smooth(formula = y ~ x,
              method = 'lm',
              se = FALSE,
              color = 'grey',
              alpha = .3) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE, 
             size = 4) +
  scale_x_continuous('log(abundance-weighted MPD)', 
                     breaks = seq(2.2,
                                  2.9,
                                  .15),
                     limits = c(2.2, 2.9)) +
  scale_y_continuous('',
                     breaks = seq(0, 3.5, 1),
                     limits = c(-0.5, 3.5)) +
  scale_color_manual(values = c('red','blue'))

loc.lrr.aw.nnd.plt <- ggplot(data = filter(forPlot, Metric == 'logAWNND'),
                             aes(x = Magnitude, y = ESCR)) + 
  geom_hline(yintercept = 0, linetype = 'dotted',
             alpha = 0.5,
             size = 2)  +
  stat_smooth(formula = y ~ x,
              method = 'lm',
              se = FALSE,
              color = 'grey',
              alpha = .3) +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE,
             size = 4) +
  scale_x_continuous('log(abundance-weighted NND)', 
                     breaks = seq(-5,
                                  -0.5,
                                  .8),
                     limits = c(-5, -0.5)) +
  scale_y_continuous('',
                     breaks = seq(0, 3.5, 1),
                     limits = c(-0.5, 3.5)) +
  scale_color_manual(values = c('red','blue'))



reg.lrr.mpd.plt <- ggplot(data = filter(forPlot, Metric == 'Regional_MPD'),
                          aes(x = Magnitude, y = ESCR)) + 
  geom_hline(yintercept = 0, linetype = 'dotted',
             alpha = 0.5,
             size = 2)  +
  geom_point(aes(color = MEPPInv),
             alpha = .4,
             show.legend = FALSE,
             size = 4) +
  scale_x_continuous('MPD', 
                     breaks = seq(13.5,
                                  15.6,
                                  0.7),
                     limits = c(13, 15.7)) +
  scale_y_continuous('',
                     breaks = seq(0, 3.5, 1),
                     limits = c(-0.5, 3.5)) +
  scale_color_manual(values = c('red','blue'))


reg.lrr.nnd.plt <- ggplot(data = filter(forPlot, Metric == 'Regional_NND'),
                          aes(x = Magnitude, y = ESCR)) + 
  geom_hline(yintercept = 0, linetype = 'dotted',
             alpha = 0.5,
             size = 2) +
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
                     limits = c(-0.5, 3.5)) +
  scale_color_manual(values = c('red','blue')) 

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
  annotate('text',
           x = .03, y = .55,
           label = 'Effect size of competition (LRR without biomass)',
           size = 6, angle = 90) +
  annotate('text', x = .52, y = .95, 
           label = 'A',
           size = 5.5) +
  annotate('text', x = .97, y = .97, label = 'B',
           size = 5) +
  annotate('text', x = .52, y = .65, label = 'C',
           size = 5) +
  annotate('text', x = .97, y = .65, label = 'D',
           size = 5) +
  annotate('text', x = .52, y = .3, label = 'E',
           size = 5) +
  annotate('text', x = .97, y = .3, label = 'F',
           size = 5) +
  annotate('text', x = .10, y = .785, label = 'Local',
           size = 5) +
  annotate('text', x = .10, y = .46, label = 'Local',
           size = 5) + 
  annotate('text', x = .09, y = .115, label = 'Regional',
           size = 5) 

ggplot2::ggsave('LRR_Novelty_NoBiomass_Appendix.png',
                path = '../Eco_Letters_Manuscript/Figures',
                height = 8.5,
                width = 12.5,
                units = 'in',
                dpi = 600)

ggsave('Fig_S3_6.pdf',
       path = '../Eco_Letters_Manuscript/SI/Figures',
       height = 7,
       width = 8,
       units = 'in',
       dpi = 600)



# Test traits for ESCR ~ Trait relationship. Analyses retained in the script, but
# no figures for the appendix. Our replication is just too low to have any power.
# But now you know we tried! ;)

trait.demo.data <- tyson$traits %>%
  setNames(c('Species', names(tyson$traits)[-1])) %>%
  filter(Species %in% demo.data$Species) %>%
  left_join(., demo.data, by = 'Species')

SLA.LM <- lm(ESCR ~ SLA*Habitat + CRBM, data = trait.demo.data)
summary(SLA.LM)

HT.LM <- lm(ESCR ~ Height*Habitat + CRBM, data = trait.demo.data)
summary(HT.LM)

Tough.LM <- lm(ESCR ~ Tough + CRBM, data = trait.demo.data)
summary(Tough.LM)



# Create figure s1.2 (r2~a for conserved traits)
traits <- c('Height', 'WoodDens', "Stemmed_Herb", 
            "Tree", "Rosette", "Vine", 
            "SubShrub", "Shrub", "Elongated_Leafy_Rhizomatous", "N_Fixer")
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

ggsave(filename = 'R2_A_Conserved_Traits_Appendix.png',
       path = '../Eco_Letters_Manuscript/Figures',
       height = 8,
       width = 12.5,
       units = 'in',
       dpi = 600)

ggsave(filename = 'Fig_S1_1.pdf',
       path = '../Eco_Letters_Manuscript/SI/Figures',
       height = 8,
       width = 12.5,
       units = 'in',
       dpi = 600)



# Trait Coverage ~ Habitat

library(rlang)
library(tidyr)

# Remove native species from community
natives <- c("Desmodium_perplexum", 
             'Geum_vernum', 
             "Symphoricarpos_orbiculatus",
             "Teucrium_canadense")

communities <- tyson$communities %>% filter(!exotic_species %in% natives)
traits <- tyson$traits
habs <- tyson$demo.data %>%
  select(Species, Habitat)

trait_com <- left_join(communities, traits, by = c("community" = "Species.Name"))

# Create catch all for the dummy variables we created
disp_traits <- c("Ballistic", "EndoZoochory", 
                 'ExoZoochory', 'Hoarding', 
                 "Myrmecochory", "Subterranean",
                 "Unassisted", "Wind",
                 "Water")
grow_forms <- c("Elongated_Leafy_Rhizomatous",
                "Rosette", "Shrub", "SubShrub",
                "Stemmed_Herb", "Tree",
                "Vine", "N_Fixer")

# Join together everything and generate long form data for ggplot
trait_hab_com <- left_join(trait_com, habs, by = c("exotic_species" = "Species")) %>%
  gather(key = "Trait", value = 'value', -c(exotic_species:Invasive, Habitat)) %>%
  mutate(trait_type = case_when(
    Trait %in% disp_traits ~ "Dispersal",
    Trait %in% grow_forms ~ "Growth_Form",
    TRUE ~ Trait
  ))

trait_hab_com$trait_type <- gsub('\\.', '_', trait_hab_com$trait_type)
trait_hab_com$trait_type <- gsub('WoodDens', 'Wood_Density', trait_hab_com$trait_type)

# summary table for ggplot!
sum_tib <- trait_hab_com %>%
  group_by(trait_type, Habitat) %>%
  summarise(
    N = n(),
    spp_cover = sum(!is.na(value))/N,
    tot_cover = sum(percentcover, na.rm = TRUE),
    cover = sum(percentcover[!is.na(value)], na.rm = TRUE)/tot_cover
  )

trait_tib <- sum_tib %>% 
  group_by(trait_type) %>%
  summarise(Average_Coverage = mean(cover))

write.csv(trait_tib, 
          file = '../Eco_Letters_Manuscript/SI/Figures/total_trait_coverage.csv',
          row.names = FALSE)

ggplot(sum_tib, aes(x = trait_type, y = cover)) + 
  geom_jitter(aes(color = Habitat),
              size = 3,
              height = 0,
              width = 0.1) +
  theme_minimal() +
  ylab("Proportion of community sampled") + 
  xlab("Trait") + 
  theme(
    axis.text = element_text(size = 15),
    axis.title.x = element_text(size = 19,
                                margin = margin(t = 20,
                                                r = 0,
                                                l = 0,
                                                b = 0)),
    axis.title.y = element_text(size = 18,
                                margin = margin(t = 0,
                                                r = 20,
                                                l = 5,
                                                b = 0))
  ) +
  scale_y_continuous(breaks = seq(-0.25, 1.25, 0.25),
                     limits = c(-0.02, 1.02))

ggsave(filename = 'Trait_Habitat_coverage.png',
       path = '../Eco_Letters_Manuscript/SI/Figures',
       height = 9,
       width = 16,
       units = 'in',
       dpi = 600)


# Various combinations of traits only that produce reasonably
# high R^2s for the appendix. ----------------------

trait.data <- tyson$traits

traits_of_interest <- c('SLA', 'Height', 
                        'Growth_Form', 
                        'Flower_Period')

picks <- seq_along(traits_of_interest)[-1]

n_comb <- function(n, r) {
  factorial(n) / (factorial(r) * factorial(n - r))
}

n_combos <- 0

for(i in picks) {
  
  n_combos <- n_combos + n_comb(length(traits_of_interest), i)
  
}

all_combos <- character(n_combos)

it <- 1

for(i in seq_len(length(traits_of_interest))[-1]) {
  
  temp <- combn(traits_of_interest, i, simplify = FALSE) %>% 
    map_chr( ~paste(.x, collapse = ', '))
  
  ll <- length(temp)
  
  all_combos[it:(it + ll - 1)] <- temp
  
  it <- it + ll
}

# This particular combination will not work as there
# are too many missing values in our data set.

rm_ind <- which(all_combos == 'Height, Flower_Period')
all_combos <- all_combos[-rm_ind]
n_combos <- n_combos - 1

# Armed with our combos, now to get extracting novelty values...

demo.data$MEPPInv <- NA
demo.data$MPD <- NA
demo.data$NND <- NA
demo.data$AWMPD <- NA
demo.data$AWNND <- NA
demo.data$Regional_MPD <- NA
demo.data$Regional_NND <- NA


out <- expand.grid(demo.data$Species,
                   all_combos,
                   stringsAsFactors = FALSE)

names(out) <- c('species', "trait_combo")

out <- mutate(out,
              ESCR    = rep(demo.data$ESCR, n_combos),
              CRBM    = rep(demo.data$CRBM, n_combos),
              MPD     = NA_real_,
              NND     = NA_real_,
              AWMPD   = NA_real_,
              AWNND   = NA_real_,
              reg_mpd = NA_real_,
              reg_nnd = NA_real_,
              MEPPInv = NA_integer_)

for(i in seq_along(all_combos)) {
  
  # Prep vector of traits for usage. retain the collapsed form to put on the
  # plot title though.
  
  traits <- all_combos[i]
  
  use_trait_names <- strsplit(traits, split = ', ') %>%
    unlist()
  
  cat('Crunching data for trait combo: ', traits, '\n')
  
  
  if("Growth_Form" %in% use_trait_names) {
    
    use_trait_names <- c(use_trait_names, grow_forms)
    use_trait_names <- use_trait_names[use_trait_names != 'Growth_form']
    
  }
  
  if("Flower_Period" %in% use_trait_names) {
    use_trait_names <- gsub("Flower_Period", "Flower.Period", use_trait_names)
  }
  
  reg_traits <- make_regional_trait_dist(trait.data, use_trait_names) %>%
    as.matrix()
  
  diag(reg_traits) <- NA_real_
  
  
  for(x in unique(demo.data$Species)) {
    
    # Need to use lines 84-88 like code (master script) to segregate by habitat
    # for regional novelty here
    
    spp_hab <- spp.list$Habitat[spp.list$Species == x]
    
    spp_hab_regex <- gsub('; ', '|', spp_hab)
    
    spp_hab_ind   <- spp.list$Species[grepl(spp_hab_regex, spp.list$Habitat)] %>%
      .[. %in% dimnames(reg_traits)[[1]]]
    
    
    # extract mpd and nnd for focal species with regional species pool
    out[out$species == x, 'reg_mpd'] <- mean(reg_traits[spp_hab_ind ,x],
                                             na.rm = TRUE)
    out[out$species == x, 'reg_nnd'] <- min(reg_traits[spp_hab_ind ,x],
                                            na.rm = TRUE)
    
    use_com <- filter(communities, exotic_species == x &
                        community %in% phylo$tip.label)
    
    # make local phylo and functional distance matrices. The functional
    # distance matrices are really just place holders, they won't be used at
    # all because a = 1. I have to make them though because the function I wrote 
    # to do this requires a functional distance matrix too because I'm an idiot.
    phyloMat <- make_local_phylo_dist(x, use_com, phylo)
    funMat <- make_local_trait_dist(x, use_com, 
                                    trait.data = trait.data,
                                    traits = use_trait_names,
                                    scale = 'scaledBYrange')
    
    AWDat <- rarefy_FPD(x, phylo.mat = phyloMat,
                        fun.mat = funMat,
                        n.rare = 11, a = 0,
                        p = 2, abundance.weighted = TRUE,
                        community.data = communities)
    # Unweighted data
    UWDat <- rarefy_FPD(x, phylo.mat = phyloMat,
                        fun.mat = funMat,
                        n.rare = 11, a = 0,
                        p = 2, abundance.weighted = FALSE,
                        community.data = NULL)  
    
    # extract all data, as well as MEPP invasive rating for focal species
    out[out$species == x & out$trait_combo == traits, 'MPD'] <- UWDat$rare.mpd
    out[out$species == x & out$trait_combo == traits, 'NND'] <- UWDat$rare.nnd
    out[out$species == x & out$trait_combo == traits, 'AWMPD'] <- log(AWDat$rare.mpd)
    out[out$species == x & out$trait_combo == traits, 'AWNND'] <- log(AWDat$rare.nnd)
    out[out$species == x & out$trait_combo == traits, 
        'MEPPInv'] <- communities[communities$exotic_species == x &
                                    communities$community == x,
                                  'Invasive']
  }
  
  temp_dat <- filter(out, trait_combo == traits) %>%
    mutate(Invasive = ifelse(MEPPInv == 1, "Yes", "No"))
  models <- list(
    temp_nnd     = lm(ESCR ~ NND + CRBM, data = temp_dat),
    temp_mpd     = lm(ESCR ~ MPD + CRBM, data = temp_dat),
    temp_awnnd   = lm(ESCR ~ AWNND + CRBM, data = temp_dat),
    temp_awmpd   = lm(ESCR ~ AWMPD + CRBM, data = temp_dat),
    temp_reg_nnd = lm(ESCR ~ reg_nnd + CRBM, data = temp_dat),
    temp_reg_mpd = lm(ESCR ~ reg_mpd + CRBM, data = temp_dat)
  )
  temp_dat$hline <- 0
  
  x2 <- seq(min(temp_dat$CRBM), 
            max(temp_dat$CRBM), 
            length.out = 14)
  x1 <- seq(min(temp_dat$NND),
            max(temp_dat$NND),
            length.out = 14)
  
  nnd_preds <- predict(models$temp_nnd, 
                       data.frame(NND = x1,
                                  CRBM = x2),
                       type = 'response',
                       interval = 'confidence',
                       se.fit = TRUE)$fit %>%
    data.frame %>% 
    cbind(., x1)
  
  x1 <- seq(min(temp_dat$MPD),
            max(temp_dat$MPD),
            length.out = 14)
  
  mpd_preds <- predict(models$temp_mpd, 
                       data.frame(MPD = x1,
                                  CRBM = x2),
                       type = 'response',
                       interval = 'confidence',
                       se.fit = TRUE)$fit %>%
    data.frame %>% 
    cbind(., x1)
  
  x1 <- seq(min(temp_dat$AWNND),
            max(temp_dat$AWNND),
            length.out = 14)
  
  awnnd_preds <- predict(models$temp_awnnd, 
                         data.frame(AWNND = x1,
                                    CRBM = x2),
                         type = 'response',
                         interval = 'confidence',
                         se.fit = TRUE)$fit %>%
    data.frame %>% 
    cbind(., x1)
  
  x1 <- seq(min(temp_dat$AWMPD),
            max(temp_dat$AWMPD),
            length.out = 14)
  
  awmpd_preds <- predict(models$temp_awmpd, 
                         data.frame(AWMPD = x1,
                                    CRBM = x2),
                         type = 'response',
                         interval = 'confidence',
                         se.fit = TRUE)$fit %>%
    data.frame %>% 
    cbind(., x1)
  
  x1 <- seq(min(temp_dat$reg_mpd),
            max(temp_dat$reg_mpd),
            length.out = 14)
  
  reg_mpd_preds <- predict(models$temp_reg_mpd, 
                           data.frame(reg_mpd = x1,
                                      CRBM = x2),
                           type = 'response',
                           interval = 'confidence',
                           se.fit = TRUE)$fit %>%
    data.frame %>% 
    cbind(., x1)
  
  x1 <- seq(min(temp_dat$reg_nnd),
            max(temp_dat$reg_nnd),
            length.out = 14)
  
  reg_nnd_preds <- predict(models$temp_reg_nnd, 
                           data.frame(reg_nnd = x1,
                                      CRBM = x2),
                           type = 'response',
                           interval = 'confidence',
                           se.fit = TRUE)$fit %>%
    data.frame %>% 
    cbind(., x1)
  
  # Alphas for the trend lines on the plots. These will determine whether or not
  # the predicted values are displayed without having to insert some voodoo into
  # the calls to ggplot()
  
  alphas <- lapply(models, function(x) {
    f <- summary(x)$fstatistic
    p <- pf(f[1], f[2], f[3], lower.tail = FALSE)
    attributes(p) <- NULL
    if(p < 0.05) {
      return(0.8)
    } else {
      return(0.01)
    }
  })
  
  plt.blank <- theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_rect(fill = NA,
                                                     size = 1.25,
                                                     color = 'black'),
                     axis.title = element_text(size = 18),
                     axis.line = element_line(size = 1.3),
                     axis.text = element_text(size = 16))
  
  loc.lrr.mpd.plt <- ggplot(data = temp_dat,
                            aes(x = MPD,
                                y = ESCR)) +
    scale_x_continuous('MPD', 
                       breaks = round(
                         seq(min(temp_dat$MPD) - 0.05,
                             max(temp_dat$MPD) + 0.05,
                             length.out = 5),
                         3
                       ),
                       limits = c(min(temp_dat$MPD), max(temp_dat$MPD))) +
    scale_y_continuous('',
                       breaks = seq(0, 3.5, 1),
                       limits = c(-.5, 3.5)) +
    geom_point(alpha = .4, 
               aes(color = Invasive),
               size = 4) + 
    scale_color_manual(
      name = 'Status',
      labels = c("Naturalized", "Invasive"),
      values = c('red','blue')
    ) +
    geom_hline(yintercept = 0, linetype = 'dotted',
               alpha = 0.5,
               size = 2) +
    geom_line(data = mpd_preds, 
              aes(x = x1, 
                  y = fit),
              color = 'black',
              alpha = alphas$temp_mpd,
              size = .8)  + 
    plt.blank 
  
  
  loc.lrr.nnd.plt <- ggplot(data = temp_dat,
                            aes(x = NND,
                                y = ESCR)) +
    scale_x_continuous('NND', 
                       breaks = round(
                         seq(min(temp_dat$NND) - 0.05,
                             max(temp_dat$NND) + 0.05,
                             length.out = 5),
                         3
                       ),
                       limits = c(min(temp_dat$NND), max(temp_dat$NND))) +
    scale_y_continuous('',
                       breaks = seq(0, 3.5, 1),
                       limits = c(-.5, 3.5)) +
    geom_point(alpha = .4, 
               aes(color = Invasive),
               size = 4,
               show.legend = FALSE) + 
    scale_color_manual(
      name = 'Status',
      labels = c("Naturalized", "Invasive"),
      values = c('red','blue')
    ) +
    geom_hline(yintercept = 0, linetype = 'dotted',
               alpha = 0.5,
               size = 2) +
    geom_line(data = nnd_preds, 
              aes(x = x1, 
                  y = fit),
              color = 'black',
              alpha = alphas$temp_nnd,
              size = .8)  + 
    plt.blank +
    ggtitle(traits)
  
  
  
  
  loc.lrr.aw.mpd.plt <- ggplot(data = temp_dat,
                               aes(x = AWMPD,
                                   y = ESCR)) +
    scale_x_continuous('Abundance Weighted MPD', 
                       breaks = round(
                         seq(min(temp_dat$AWMPD) - 0.05,
                             max(temp_dat$AWMPD) + 0.05,
                             length.out = 5),
                         3
                       ),
                       limits = c(min(temp_dat$AWMPD), max(temp_dat$AWMPD))) +
    scale_y_continuous('',
                       breaks = seq(0, 3.5, 1),
                       limits = c(-.5, 3.5)) +
    geom_point(alpha = .4, 
               aes(color = Invasive),
               size = 4,
               show.legend = FALSE) + 
    scale_color_manual(
      name = 'Status',
      labels = c("Naturalized", "Invasive"),
      values = c('red','blue')
    ) +
    geom_hline(yintercept = 0, linetype = 'dotted',
               alpha = 0.5,
               size = 2) +
    geom_line(data = awmpd_preds, 
              aes(x = x1, 
                  y = fit),
              color = 'black',
              alpha = alphas$temp_awmpd,
              size = .8)  + 
    plt.blank 
  
  loc.lrr.aw.nnd.plt <- ggplot(data = temp_dat,
                               aes(x = AWNND,
                                   y = ESCR)) +
    scale_x_continuous('Abundance Weighted NND', 
                       breaks = round(
                         seq(min(temp_dat$AWNND) - 0.05,
                             max(temp_dat$AWNND) + 0.05,
                             length.out = 5), 
                         3
                       ),
                       limits = c(min(temp_dat$AWNND), max(temp_dat$AWNND))) +
    scale_y_continuous('',
                       breaks = seq(0, 3.5, 1),
                       limits = c(-.5, 3.5)) +
    geom_point(alpha = .4, 
               aes(color = Invasive),
               size = 4,
               show.legend = FALSE) + 
    scale_color_manual(
      name = 'Status',
      labels = c("Naturalized", "Invasive"),
      values = c('red','blue')
    ) +
    geom_hline(yintercept = 0, linetype = 'dotted',
               alpha = 0.5,
               size = 2) +
    geom_line(data = awnnd_preds, 
              aes(x = x1, 
                  y = fit),
              color = 'black',
              alpha = alphas$temp_awnnd,
              size = .8)  + 
    plt.blank 
  
  
  reg.lrr.mpd.plt <- ggplot(data = temp_dat,
                            aes(x = reg_mpd,
                                y = ESCR)) +
    scale_x_continuous('Regional Scale MPD', 
                       breaks = round(
                         seq(min(temp_dat$reg_mpd) - 0.05,
                             max(temp_dat$reg_mpd) + 0.05,
                             length.out = 5),
                         3
                       ),
                       limits = c(min(temp_dat$reg_mpd), 
                                  max(temp_dat$reg_mpd))) +
    scale_y_continuous('',
                       breaks = seq(0, 3.5, 1),
                       limits = c(-.5, 3.5)) +
    geom_point(alpha = .4, 
               aes(color = Invasive),
               size = 4,
               show.legend = FALSE) + 
    scale_color_manual(
      name = 'Status',
      labels = c("Naturalized", "Invasive"),
      values = c('red','blue')
    ) +
    geom_hline(yintercept = 0, linetype = 'dotted',
               alpha = 0.5,
               size = 2) +
    geom_line(data = reg_mpd_preds, 
              aes(x = x1, 
                  y = fit),
              color = 'black',
              alpha = alphas$temp_reg_mpd,
              size = .8)  + 
    plt.blank 
  
  
  reg.lrr.nnd.plt <- ggplot(data = temp_dat,
                            aes(x = reg_nnd,
                                y = ESCR)) +
    scale_x_continuous('Regional Scale NND', 
                       breaks = round(
                         seq(min(temp_dat$reg_nnd) - 0.005,
                             max(temp_dat$reg_nnd) + 0.005,
                             length.out = 5),
                         3
                       ),
                       limits = c(min(temp_dat$reg_nnd), max(temp_dat$reg_nnd))) +
    scale_y_continuous('',
                       breaks = seq(0, 3.5, 1),
                       limits = c(-.5, 3.5)) +
    geom_point(alpha = .4, 
               aes(color = Invasive),
               size = 4,
               show.legend = FALSE) + 
    scale_color_manual(
      name = 'Status',
      labels = c("Naturalized", "Invasive"),
      values = c('red','blue')
    ) +
    geom_hline(yintercept = 0, linetype = 'dotted',
               alpha = 0.5,
               size = 2) +
    geom_line(data = reg_nnd_preds, 
              aes(x = x1, 
                  y = fit),
              color = 'black',
              alpha = alphas$temp_reg_nnd,
              size = .8)  + 
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
  
  
  
  fp_name <- gsub(', ', '_', traits)
  
  ggsave(filename = glue('figure_1_{fp_name}.png'),
         path = '../Eco_Letters_Manuscript/Figures/trait_appendix_plots',
         height = 8.5,
         width = 12.5,
         units = 'in',
         dpi = 600)
}

