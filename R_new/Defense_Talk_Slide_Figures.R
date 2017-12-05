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


# Figures to demonstrate hypotheses------------
# novelty and invasive class
x <- runif(80, 0, 200) %>% sort() # MPD
y <- c(rbinom(40, 1, 0.2), rbinom(40, 1, 0.8)) # Invasiveness

positiveGLM <- glm(y ~ x, family = binomial())
plot(y ~ x, type = 'n',
     xlab = 'Novelty',
     ylab = 'Invasive Status',
     axes = FALSE,
     cex.lab = 1.5)
box(lwd = 1)
lines(x, predict(positiveGLM, type = 'response'))

y <- c(rbinom(40, 1, 0.8), rbinom(40, 1, 0.2)) # Invasiveness

negativeGLM <- glm(y ~ x, family = binomial())
plot(y ~ x, type = 'n',
     xlab = 'Novelty',
     ylab = 'Invasive Status',
     axes = FALSE,
     cex.lab = 1.5)
box(lwd = 1)
lines(x, predict(negativeGLM, type = 'response'))

# Linear LRR Simulation
x <- runif(14, 0, 240)
y1 <- 5 + -.025 * x + rnorm(14, 0, 1.5) # phylogeny only. adding more error for demonstration purposes
PhyloLM <- lm(y1 ~ x)
PhyloPreds <- predict(PhyloLM, se.fit = TRUE,
                      interval = 'confidence',
                      level = 0.95)

forPlot <- data.frame(Novelty = rep(x, 3),
                      ESCR = rep(y1, 3),
                      Type = c(rep('Phylo', 14),
                               rep('Poor-Phylo',14),
                               rep('Good-Phylo',14)),
                      UpCI = c(PhyloPreds$fit[ ,1] + PhyloPreds$se.fit,
                               PhyloPreds$fit[ ,1] + 0.3 * PhyloPreds$se.fit, # less error for poorly conserved traits
                               PhyloPreds$fit[ ,1] + 0.75 * PhyloPreds$se.fit), # slightly less error for well conserved traits
                      LoCI = c(PhyloPreds$fit[ ,1] - PhyloPreds$se.fit,
                               PhyloPreds$fit[ ,1] - 0.3 * PhyloPreds$se.fit,
                               PhyloPreds$fit[ ,1] - 0.75 * PhyloPreds$se.fit),
                      stringsAsFactors = FALSE)

plt.blank <-  theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    panel.border = element_rect(size = 2,
                                                fill = NA),
                    axis.line = element_line(colour = "black"),
                    axis.title = element_text(size = 20))

Plt1 <- ggplot(data = filter(forPlot, Type == 'Phylo'),
               aes(x = Novelty, y = ESCR)) + 
  plt.blank +
  geom_ribbon(aes(ymin = LoCI, 
                  ymax = UpCI),
              fill = 'red',
              alpha = 0.6) +
  geom_abline(aes(slope = coef(PhyloLM)[2],
                  intercept = coef(PhyloLM)[1])) +
  scale_y_continuous('Effect Size of Competition',
                     breaks = NULL,
                     labels = NULL) +
  scale_x_continuous('Novelty',
                     breaks = NULL,
                     labels = NULL)

Plt1

Plt2 <- ggplot(data = filter(forPlot, Type == 'Good-Phylo'),
               aes(x = Novelty, y = ESCR)) + 
  plt.blank +
  geom_ribbon(aes(ymin = LoCI, 
                  ymax = UpCI),
              fill = 'grey40',
              alpha = 0.6) +
  geom_ribbon(data = filter(forPlot, Type == 'Phylo'),
              aes(ymin = LoCI,
                  ymax = UpCI),
              fill = 'red',
              alpha = 0.2) +
  geom_abline(aes(slope = coef(PhyloLM)[2],
                  intercept = coef(PhyloLM)[1])) +
  scale_y_continuous('Effect Size of Competition',
                     breaks = NULL,
                     labels = NULL) +
  scale_x_continuous('Novelty',
                     breaks = NULL,
                     labels = NULL)
Plt2

Plt3 <- ggplot(data = filter(forPlot, Type == 'Poor-Phylo'),
               aes(x = Novelty, y = ESCR)) + 
  plt.blank +
  geom_ribbon(aes(ymin = LoCI,
                  ymax = UpCI),
              fill = 'black',
              alpha = 0.8) + 
  geom_ribbon(data = filter(forPlot, Type == 'Good-Phylo'),
              aes(ymin = LoCI, 
                  ymax = UpCI),
              fill = 'grey40',
              alpha = 0.5) +
  geom_ribbon(data = filter(forPlot, Type == 'Phylo'),
              aes(ymin = LoCI,
                  ymax = UpCI),
              fill = 'red',
              alpha = 0.25) +
  geom_abline(aes(slope = coef(PhyloLM)[2],
                  intercept = coef(PhyloLM)[1])) +
  scale_y_continuous('Effect Size of Competition',
                     breaks = NULL,
                     labels = NULL) +
  scale_x_continuous('Novelty',
                     breaks = NULL,
                     labels = NULL)
Plt3

phylo <- tyson$phylo
kstra <- filter(communities, exotic_species == 'Kummerowia_striata')

kstraPhylo <- drop.tip(phylo, setdiff(phylo$tip.label, kstra$community))

plot(kstraPhylo, show.tip.label = FALSE)
