# Script to bootstrap regressions to handle uncertainty in ESCR.

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
  
  # If we don't have habitat info for some reason, just use all of the species
  # This occurs for the Iris genus and couple other non-focal species, so will
  # not affect broader results anyway. 
  
  if(all(is.na(spp_hab_ind))) {
    spp_hab_ind <- seq(1, dim(exotics)[1], by = 1)
  }
  
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
  
  betas_mpd    <- coef(lm(escr ~ MPD + CRBM, data = reg_data))
  betas_nnd    <- coef(lm(escr ~ NND + CRBM, data = reg_data))
  betas_awmpd  <- coef(lm(escr ~ logAWMPD + CRBM, data = reg_data))
  betas_awnnd  <- coef(lm(escr ~ logAWNND + CRBM, data = reg_data))
  betas_re_mpd <- coef(lm(escr ~ Regional_MPD + CRBM, data = reg_data))
  betas_re_nnd <- coef(lm(escr ~ Regional_NND + CRBM, data = reg_data))
  
  out <- list(mpd = betas_mpd,
              nnd = betas_nnd,
              awmpd = betas_awmpd,
              awnnd = betas_awnnd,
              reg_mpd = betas_re_mpd,
              reg_nnd = betas_re_nnd)
  
  return(out)
  
}

nov_data <- demo.data %>%
  select(Species, CRBM, MPD:logAWNND)

output <- list(
  mpd = data.frame(int = numeric(1001L),
                   nov = numeric(1001L),
                   bio = numeric(1001L),
                   met = "MPD",
                   stringsAsFactors = FALSE),
  nnd = data.frame(int = numeric(1001L),
                   nov = numeric(1001L),
                   bio = numeric(1001L),
                   met = "NND",
                   stringsAsFactors = FALSE),
  awmpd = data.frame(int = numeric(1001L),
                     nov = numeric(1001L),
                     bio = numeric(1001L),
                     met = "AW-MPD",
                     stringsAsFactors = FALSE),
  awnnd = data.frame(int = numeric(1001L),
                     nov = numeric(1001L),
                     bio = numeric(1001L),
                     met = "AW-NND",
                     stringsAsFactors = FALSE),
  reg_mpd = data.frame(int = numeric(1001L),
                       nov = numeric(1001L),
                       bio = numeric(1001L),
                       met = "Regional MPD",
                       stringsAsFactors = FALSE),
  reg_nnd = data.frame(int = numeric(1001L),
                       nov = numeric(1001L),
                       bio = numeric(1001L),
                       met = "Regional NND",
                       stringsAsFactors = FALSE)
)

for(i in 1:1001) {
  
  temp_betas <- escr_reg(all_lams, nov_data, i)
  
  output$mpd[i, 1:3] <- temp_betas$mpd
  output$nnd[i, 1:3] <- temp_betas$nnd
  output$awmpd[i, 1:3] <- temp_betas$awmpd
  output$awnnd[i, 1:3] <- temp_betas$awnnd
  output$reg_mpd[i, 1:3] <- temp_betas$reg_mpd
  output$reg_nnd[i, 1:3] <- temp_betas$reg_nnd
  
}

for_plot <- do.call("rbind", output)

ggplot(for_plot, aes(x = nov)) + 
  geom_histogram() +
  facet_wrap(~met, scales = "free_x") + 
  geom_vline(xintercept = 0,
             color = 'red') +
  theme_bw()

ggsave(filename = "../Eco_Letters_Manuscript/Figures/Figure_2_hist.png")

reg_fig_data <- lapply(output,
                       function(x) { 
                         obs <- t(x[1, 1:3]) %>%
                           as.data.frame()
                         sorting <- x[2:1001, 1:3]
                         to_bind <- sorting[order(sorting$nov), ]
                         to_bind <- t(to_bind[c(25, 975), 1:3])
                         
                         out <- data.frame(to_bind) %>%
                         cbind(obs, .) %>%
                           setNames(
                             c("Value", "Lo_CI", "Up_CI")
                           ) %>%
                           mutate(Coefficient = factor(c("Intercept",
                                                         "Distinctiveness",
                                                         "Biomass"),
                                                       levels = c("Intercept",
                                                                  "Distinctiveness",
                                                                  "Biomass"),
                                                       ordered = TRUE))
                         
                         return(out)
                         
                       }) %>%
  do.call("rbind", args = .) %>%
  mutate(Metric = c(
    rep("MPD", 3),
    rep("NND", 3),
    rep("AW-MPD", 3),
    rep("AW-NND", 3),
    rep("Regional MPD", 3),
    rep("Regional NND", 3)
  ))




ggplot(reg_fig_data, aes(x = Coefficient, y = Value)) +
  facet_wrap(~Metric,
             scales = 'free_y') + 
  geom_point(size = 2) +
  geom_linerange(aes(ymin = Lo_CI,
                     ymax = Up_CI),
                 size = 1.25) + 
  theme_bw() +
  geom_hline(yintercept = 0, color = 'red') 

ggsave(filename = "../Eco_Letters_Manuscript/Figures/Figure_2_coef_barplot.png")


