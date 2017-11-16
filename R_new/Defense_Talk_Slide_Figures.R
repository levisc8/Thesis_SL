# Figures and misc stuff for defense slides

rm(list = ls(all = T))
library(FunPhylo)
library(pez)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(cowplot)

data(tyson)

# figure to emphasize Kummerowia and Ligustrum differences
forPlot <- tyson$demo.data %>% filter(Species == 'Ligustrum_obtusifolium' |
                                        Species == 'Kummerowia_striata')

communities <- tyson$communities
phylo <- tyson$phylo


for(x in unique(forPlot$Species)) {
  phylomat <- make_local_phylo_dist(x, communities, phylo)
  diag(phylomat) <- NA

  forPlot[forPlot$Species == x, 'MPD'] <- mean(phylomat[ ,x],
                                               na.rm = TRUE)
  
  forPlot[forPlot$Species == x, 'NND'] <- min(phylomat[ ,x],
                                              na.rm = TRUE)
  
  
}

# make ggplot Figures
KstraPlot <- ggplot(filter(forPlot, Species == 'Kummerowia_striata'),
                    aes(x = MPD, y = ESCR2)) +
  geom_point(color = 'blue', alpha = 0.4, 
             show.legend = FALSE, size = 4) +
  scale_x_continuous('MPD', 
                     breaks = seq(140, 280,40),
                     limits = c(139, 281)) + 
  scale_y_continuous('Effect size of competition',
                     breaks = seq(0,3,1),
                     limits = c(-0.5, 3.1)) + 
  annotate('text', x = 152.5, y = 2.9,
           label = 'Kummerowia striata',
           size = 5)
  
KstraPlot

LigPlot <- KstraPlot + 
  geom_point(data = filter(forPlot, Species == 'Ligustrum_obtusifolium'),
             aes(x = MPD, y = ESCR2),
             size = 4,
             color =  'blue',
             alpha = 0.4) + 
  annotate('text', x = 256.5, y = 0.2, label = 'Ligustrum obtusifolium',
           size = 5) + 
  geom_hline(yintercept = 0, linetype = 'dotted')
LigPlot
