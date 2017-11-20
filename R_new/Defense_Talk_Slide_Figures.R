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
           size = 5) +
  geom_hline(yintercept = 0, linetype = 'dotted')

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

# create color coded Tyson scale phylogeny figure
library(ggtree)
library(viridis)
spp.list <- tyson$spp.list
spp.list$Species <- gsub('\\.', '-', spp.list$Species)
forPhylo <- filter(spp.list, Monocot == 0) %>% as_tibble()
demo.data <- tyson$demo.data

# add group labels. Will be important for color coding with ggtree
forPhylo$PhyloGroup <- NA
forPhylo$PhyloGroup <- ifelse(forPhylo$Exotic == 0, 'Native',
                              ifelse(forPhylo$Invasive == 0, 'Exotic', 'Invasive'))

forPhylo$PhyloGroup[forPhylo$Species %in% unique(demo.data$Species)] <- 'Focal_Species'
forPhylo <- forPhylo[!duplicated(forPhylo$Species) & 
                       forPhylo$Species %in% phylo$tip.label, ] # remove duplicated species 

# create grouping list for plotting
groups <- list(Natives = forPhylo$Species[forPhylo$PhyloGroup == 'Native'],
               Exotics = forPhylo$Species[forPhylo$PhyloGroup == 'Exotic'],
               Invasives = forPhylo$Species[forPhylo$PhyloGroup == 'Invasive'],
               Focal.Species = forPhylo$Species[forPhylo$PhyloGroup == 'Focal_Species'])

# add group variables to the phylogeny
forPlotting <- groupOTU(phylo, groups)

# select colors from plotting. the hexcode values come from a call to
# viridis::viridis(4, alpha = 0.9, direction = -1) and removing the 3rd
# (the contrast wasn't as good as I'd hoped).
colors <- c("#FDE725E6", "#35B779E6", "red", "#440154E6")

# Create the ggplot object and then annotate it with cowplot below
Tree <- ggtree(forPlotting, layout = 'circular', 
               aes(color = group)) + 
  theme(legend.position = 'right') + 
  scale_color_manual('',
                     breaks = c('Natives',
                                'Focal.Species',
                                'Invasives',
                                'Exotics'),
                     labels = c('Native',
                                'Focal Species',
                                'Invasive',
                                'Exotic/Non-Invasive'),
                     values = colors)
# Using cowplot so I can offset the title with it's really simple interface 
ggdraw() + 
  draw_plot(Tree, x = 0, y = 0, height = 1, width = 1) + 
  annotate('text', x = .825, y = .65,  
           label = 'Phylogeny of the \nTRC Species Pool',
           size = 6)

# Save the output and have a look!
# ggsave('../Figures/TRC_Phylo.png', height = 8, width = 8, units = 'in')
