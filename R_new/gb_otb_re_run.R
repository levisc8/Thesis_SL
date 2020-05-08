# Re-run analyses with Smith Trees to see if any difference is made

# Script to produce all of the analyses and figures in the main text of the paper.

library(FunPhylo)
library(pez)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(grid)
library(gridExtra)
library(tidyr)
library(stringr)

data(tyson)

# do a little data munging and variable conversion

phylo <- tyson$phylo_gb

spp.list <- tyson$spp.list %>%
  filter(Species %in% phylo$tip.label)

for(i in 2:4){
  spp.list[ ,i] <- as.numeric(spp.list[ ,i])
}


communities <- tyson$communities %>%
  filter(community %in% phylo$tip.label)
demo.data <- tyson$demo.data

# Make sure all species in our trait data can be found in the phylogney
traits <- tyson$traits %>%
  filter(Species.Name %in% phylo$tip.label)

# Get habitat info for every species in the community data set

temp <- demo.data %>%
  dplyr::select(Species, Habitat) %>%
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
  
  
  exotics[exotics$Species == x, 'MPD'] <- mean(regionalTysonDists[, x],
                                               na.rm = TRUE)
  exotics[exotics$Species == x, 'NND'] <- min(regionalTysonDists[ , x],
                                              na.rm = TRUE)
  exotics[exotics$Species == x, 'Focal'] <- ifelse(x %in% demo.data$Species,
                                                   'Y', 
                                                   'N')
}


# local scale effect sizes, phylo only---------------
# phylogenetic only to start
# trait * phylogeny comes later

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
 
  phyloMat <- make_local_phylo_dist(x, communities, phylo)
  funMat <- make_local_trait_dist(x, communities, 
                                  trait.data = traits,
                                  traits = names(traits)[-1],
                                  scale = 'scaledBYrange')
  
  Loc.Mat <- phyloMat
  diag(Loc.Mat) <- NA  
  
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

genera <- vapply(demo.data$Species,
                 function(x) strsplit(x, "_")[[1]][1],
                 character(1L))

paths <- paste('../Data/bootstrap_lambdas/', genera,'_lambdas.rds', sep = "")

all_lams <- lapply(paths, function(x) readRDS(x))
names(all_lams) <- genera

# First regression, then initialize output data

escr_reg <- function(demo_data, novelty_data, iteration) {
  
  escr <- lapply(demo_data, 
                 function(x, it) {
                   log(x$lambda_cr[it] + 0.5) - log(x$lambda_c[it] + 0.5)
                 },
                 it = iteration) %>%
    unlist()
  
  reg_data <- cbind(escr, novelty_data)
  
  model_mpd    <- lm(escr ~ MPD + CRBM, data = reg_data)
  model_nnd    <- lm(escr ~ NND + CRBM, data = reg_data)
  model_awmpd  <- lm(escr ~ logAWMPD + CRBM, data = reg_data)
  model_awnnd  <- lm(escr ~ logAWNND + CRBM, data = reg_data)
  model_re_mpd <- lm(escr ~ Regional_MPD + CRBM, data = reg_data)
  model_re_nnd <- lm(escr ~ Regional_NND + CRBM, data = reg_data)
  
  betas_mpd    <- c(coef(model_mpd), summary(model_mpd)$adj.r.squared)
  betas_nnd    <- c(coef(model_nnd), summary(model_nnd)$adj.r.squared)
  betas_awmpd  <- c(coef(model_awmpd), summary(model_awmpd)$adj.r.squared)
  betas_awnnd  <- c(coef(model_awnnd), summary(model_awnnd)$adj.r.squared)
  betas_re_mpd <- c(coef(model_re_mpd), summary(model_re_mpd)$adj.r.squared)
  betas_re_nnd <- c(coef(model_re_nnd), summary(model_re_nnd)$adj.r.squared)
  
  out <- list(mpd = betas_mpd,
              nnd = betas_nnd,
              awmpd = betas_awmpd,
              awnnd = betas_awnnd,
              reg_mpd = betas_re_mpd,
              reg_nnd = betas_re_nnd)
  
  return(out)
  
}

nov_data <- demo.data %>%
  dplyr::select(Species, CRBM, MPD:logAWNND)

output <- list(
  mpd = data.frame(int = numeric(1001L),
                   nov = numeric(1001L),
                   bio = numeric(1001L),
                   r2  = numeric(1001L),
                   met = "Small Grain MPD",
                   stringsAsFactors = FALSE),
  nnd = data.frame(int = numeric(1001L),
                   nov = numeric(1001L),
                   bio = numeric(1001L),
                   r2  = numeric(1001L),
                   met = "Small Grain NND",
                   stringsAsFactors = FALSE),
  awmpd = data.frame(int = numeric(1001L),
                     nov = numeric(1001L),
                     bio = numeric(1001L),
                     r2  = numeric(1001L),
                     met = "Small Grain AW-MPD",
                     stringsAsFactors = FALSE),
  awnnd = data.frame(int = numeric(1001L),
                     nov = numeric(1001L),
                     bio = numeric(1001L),
                     r2  = numeric(1001L),
                     met = "Small Grain AW-NND",
                     stringsAsFactors = FALSE),
  reg_mpd = data.frame(int = numeric(1001L),
                       nov = numeric(1001L),
                       bio = numeric(1001L),
                       r2  = numeric(1001L),
                       met = "Large Grain MPD",
                       stringsAsFactors = FALSE),
  reg_nnd = data.frame(int = numeric(1001L),
                       nov = numeric(1001L),
                       bio = numeric(1001L),
                       r2  = numeric(1001L),
                       met = "Large Grain NND",
                       stringsAsFactors = FALSE)
)

for(i in 1:1001) {
  
  temp_betas <- escr_reg(all_lams, nov_data, i)
  
  output$mpd[i, 1:4] <- temp_betas$mpd
  output$nnd[i, 1:4] <- temp_betas$nnd
  output$awmpd[i, 1:4] <- temp_betas$awmpd
  output$awnnd[i, 1:4] <- temp_betas$awnnd
  output$reg_mpd[i, 1:4] <- temp_betas$reg_mpd
  output$reg_nnd[i, 1:4] <- temp_betas$reg_nnd
  
}

for_plot <- do.call("rbind", output)

reg_fig_data <- lapply(output,
                       function(x) { 
                         
                         obs <- t(x[1, 1:4]) %>%
                           as.data.frame()
                         sorting <- x[2:1001, 1:3]
                         to_bind <- sorting[order(sorting$nov), ]
                         to_bind <- t(to_bind[c(25, 975), 1:3])
                         
                         sorting_r2 <- x[2:1001, 4]
                         r2_bind <- sorting_r2[order(sorting_r2)]
                         r2_bind <- t(r2_bind[c(25, 975)])
                         
                         to_bind <- rbind(to_bind, r2_bind)
                         
                         out <- data.frame(to_bind) %>%
                           cbind(obs, .) %>%
                           setNames(
                             c("Value", "Lo_CI", "Up_CI")
                           ) %>%
                           mutate(Coefficient = factor(c("Intercept",
                                                         "Distinctiveness",
                                                         "Biomass",
                                                         "paste(R[adj]^2)"),
                                                       levels = c("Intercept",
                                                                  "Distinctiveness",
                                                                  "Biomass",
                                                                  "paste(R[adj]^2)"),
                                                       ordered = TRUE))
                         
                         return(out)
                         
                       }) %>%
  do.call("rbind", args = .) %>%
  mutate(Metric = c(
    rep("Small Grain MPD", 4),
    rep("Small Grain NND", 4),
    rep("Small Grain AW-MPD", 4),
    rep("Small Grain AW-NND", 4),
    rep("Large Grain MPD", 4),
    rep("Large Grain NND", 4)
  )) %>%
  filter(! Coefficient %in% c("Intercept", "Biomass"))

for_plot$Metric <- factor(for_plot$met,
                          levels = c(
                            "Small Grain NND",
                            "Small Grain MPD",
                            "Small Grain AW-NND",
                            "Small Grain AW-MPD",
                            "Large Grain NND",
                            "Large Grain MPD"
                          ),
                          ordered = TRUE)

r2s   <- filter(reg_fig_data, Coefficient != "Distinctiveness")

beta_hists <- ggplot(for_plot, aes(x = nov)) +
  geom_histogram(bins = 75) +
  facet_wrap(~ Metric,
             scales = "free_x",
             nrow = 3,
             ncol = 2) + 
  theme_bw() +
  geom_vline(xintercept = 0, color = 'red', size = 1.5, linetype = 'dashed')  + 
  labs(x = "Coefficient Value",
       y = "Count") + 
  theme(strip.text = element_text(size = 14),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13))

r2_plot <- ggplot(r2s, aes(x = Metric, y = Value)) +
  geom_point(size = 3) +
  geom_linerange(aes(ymin = Lo_CI, ymax = Up_CI), size = 1.25) +
  theme_bw() +
  scale_y_continuous(limits = c(-0.05, 1.05)) +
  coord_flip() +
  labs(x = "",
       y = expression(paste(R[adj]^2))) +
  theme(axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 13))

png(filename = "../Eco_Letters_Manuscript/Figures/Figure_2_gb_OTB.png",
    height = 8.5,
    width = 12.5,
    units = 'in',
    res = 600)

grid.arrange(beta_hists, r2_plot, nrow = 1, ncol = 2)

dev.off()

pdf(file = "../Eco_Letters_Manuscript/Figures/Figure_2_gb_OTB.pdf",
    height = 8.5,
    width = 13.5)

grid.arrange(beta_hists, r2_plot, nrow = 1, ncol = 2)

dev.off()

#  write table to csv so we can include in paper

write.csv(reg_fig_data, 
          file = '../Eco_Letters_Manuscript/Figures/Phylo_models_output_gb_OTB.csv', 
          na = "",
          row.names = FALSE)


# functional phylogenetic regressions-------------
# next, determine best models using latest version of invasives FPD
# https://sam-levin.shinyapps.io/Invasives_FPD/

# ----- leaving RStudio Beep Boop, report back w/ best models-----

# ----- Beep Boop and we're back. Best model was Phylogeny + SLA, Ht, Flower period,
# and leaf toughness with an a-value of 0.3-0.4 (best a's vary due to rarefying, 
# but the max R^2 is always ~0.8 and it's always in this range). 

demo.data <- arrange(demo.data, desc(ESCR))

traits <- c('Height', 'SLA', 'Tough', 'Flower.Period')
trait.data <- tyson$traits %>%
  filter(Species.Name %in% phylo$tip.label)


# Remove monocots for this part because we have no functional
# trait information for them

monocots <- tyson$spp.list %>%
  filter(Monocot == 1)

communities <- filter(communities, ! community %in% monocots$Species)

a_seq <- seq(0, 1, .025)
R2dat <- data.frame(A = a_seq,
                    NND = rep(NA_real_, length(a_seq)),
                    MPD = rep(NA_real_, length(a_seq)))
mod.data <- demo.data

mod.data[ , paste0('nna_', a_seq)] <- NA_real_
mod.data[ , paste0('mpa_', a_seq)] <- NA_real_

for(x in unique(demo.data$Species)){
  cat('Calculating FPD for species: ', x, '\n')
  
  
  use_com <- filter(communities, exotic_species == x &
                      community %in% phylo$tip.label)
  
  # phylo and functional distance matrices
  phylo.mat <- make_local_phylo_dist(x, 
                                     use_com, 
                                     phylo)
  fun.mat <- make_local_trait_dist(x, 
                                   use_com, 
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

ggsave(filename = 'R2_A_gb_OTB_For_Manuscript.png',
       path = '../Eco_Letters_Manuscript/Figures',
       height = 8,
       width = 12.5,
       units = 'in',
       dpi = 600)

ggsave(filename = "Figure_2_gb_OTB.pdf",
       path = '../Eco_Letters_Manuscript/Figures',
       height = 8.5,
       width = 12.5,
       units = 'in',
       dpi = 600)
