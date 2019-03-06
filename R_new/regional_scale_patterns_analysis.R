

# mixed effects model with novelty metric and scale
# as interaction terms. We'll then perform a likelihood ratio test with
# successively more complicated models and see which is preferred. The key point 
# here is whether the model with scale * novelty is better, and whether that interaction
# is significant. If not, then we can conclude that novelty isn't necessarily important
# at either scale

exotics$Scale <- 'Regional'
All.Local.exotics$Scale <- 'Local'

# prepare regional data
exoticsForJoin <- exotics %>%
  select(Species,Exotic, Invasive,
         MPD, NND, Focal, Scale)

# prepare local data and join with regional
allExotics <- All.Local.exotics %>%
  select(community, alien, Invasive, MPD.Local, 
         NND.Local, Focal, Scale) %>%
  setNames(c('Species', 'Exotic', 'Invasive',
             'MPD', 'NND', 'Focal', 'Scale')) %>%
  rbind(exoticsForJoin)
allExotics$Scale <- as.factor(allExotics$Scale)

# Fit the models!
FullScaleRegMPD <- glmer(Invasive ~ MPD * Scale + (1|Species),
                         data = allExotics, 
                         family = binomial)
MidScaleRegMPD <- glmer(Invasive ~ MPD + Scale + (1|Species),
                        data = allExotics, 
                        family = binomial)
NullScaleRegMPD <- glmer(Invasive ~ MPD + (1|Species),
                         data = allExotics, 
                         family = binomial)

FullScaleRegNND <- glmer(Invasive ~ NND * Scale + (1|Species),
                         data = allExotics, 
                         family = binomial)
MidScaleRegNND <- glmer(Invasive ~ NND + Scale + (1|Species),
                        data = allExotics, 
                        family = binomial)
NullScaleRegNND <- glmer(Invasive ~ NND + (1|Species),
                         data = allExotics, 
                         family = binomial)

# I get warnings about identifiability and convergence with 4/6 models.
# Re-fit all of the models as they will not be comparable with variables on
# different scales.

allExotics <- mutate(allExotics, 
                     MPDScaled = scale(MPD, center = TRUE, scale = TRUE),
                     NNDScaled = scale(NND, center = TRUE, scale = TRUE))


FullScaleRegMPDRescale <- glmer(Invasive ~ MPDScaled * Scale + (1|Species),
                                data = allExotics, 
                                family = binomial)
MidScaleRegMPDRescale <- glmer(Invasive ~ MPDScaled + Scale + (1|Species),
                               data = allExotics, 
                               family = binomial)
NullScaleRegMPDRescale <- glmer(Invasive ~ MPDScaled + (1|Species),
                                data = allExotics, 
                                family = binomial)

FullScaleRegNNDRescale <- glmer(Invasive ~ NNDScaled * Scale + (1|Species),
                                data = allExotics, 
                                family = binomial)
MidScaleRegNNDRescale <- glmer(Invasive ~ NNDScaled + Scale + (1|Species),
                               data = allExotics, 
                               family = binomial)
NullScaleRegNNDRescale <- glmer(Invasive ~ NNDScaled + (1|Species),
                                data = allExotics, 
                                family = binomial)

# Well, now only 1 is having issues. I'll try running it with more iterations
# next
MidScaleRegNNDRescale <- glmer(Invasive ~ NNDScaled + Scale + (1|Species),
                               data = allExotics, 
                               family = binomial,
                               control = glmerControl(optCtrl = list(maxfun = 1e6)))

# Nothing changes. going to try restarting from the end points of the last one.
NewStart <- getME(MidScaleRegNNDRescale, c('theta', 'fixef'))
MidScaleRegNNDRescale_2 <- update(MidScaleRegNNDRescale,
                                  start = NewStart,
                                  control = glmerControl(optCtrl = list(maxfun = 1e6)))


# Onward to the summaries
message('Mixed effects models for MPD patterns\n')
summary(FullScaleRegMPDRescale)
summary(MidScaleRegMPDRescale)
summary(NullScaleRegMPDRescale)

# Based on the summaries, it seems pretty clear that scale makes no difference
# for MPD

message('Mixed effects models for NND patterns\n')

summary(FullScaleRegNNDRescale)
summary(MidScaleRegNNDRescale_2)
summary(NullScaleRegNNDRescale)

# In this case, it is not so clear from the summaries. On to the likelihood
# ratio tests

NndLrt <- anova(NullScaleRegNNDRescale, 
                MidScaleRegNNDRescale_2,
                FullScaleRegNNDRescale)

MpdLrt <- anova(NullScaleRegMPDRescale,
                MidScaleRegMPDRescale,
                FullScaleRegMPDRescale)

NndLrt
MpdLrt

# Regional scale logistic regressions - for the manuscript text
logisticModelTable <- data.frame(ModelType = c(rep('Nearest Neighbor Distance', 5),
                                               rep('Mean Pairwise Distance', 3)),
                                 rbind(tidy(FullScaleRegNNDRescale),
                                       tidy(NullScaleRegMPDRescale))) %>%
  map_if(is.numeric, ~round(.x, 4)) %>%
  data.frame() %>%
  setNames(c('Model Type', 'Parameter', 'Estimate', 'Std. Error', 'Statistic',
             'p-value', 'Group')) 


logisticModelTable$`p-value` <- lapply(logisticModelTable$`p-value`,
                                       FUN = function (x) add_stars(x)) %>%
  unlist()


write.csv(logisticModelTable,
          file = '../Eco_Letters_Manuscript/Figures/Logistic_regression_outputs.csv',
          na = "", row.names = FALSE)

# Regional scale logistic regressions - for the appendix

AppendixModelTable <- data.frame(ModelType = c(rep('Nearest Neighbor Distance', 7),
                                               rep('Mean Pairwise Distance', 9)),
                                 rbind(tidy(NullScaleRegNNDRescale),
                                       tidy(MidScaleRegNNDRescale_2),
                                       tidy(MidScaleRegMPDRescale),
                                       tidy(FullScaleRegMPDRescale)))%>%
  map_if(is.numeric, ~round(.x, 4)) %>%
  data.frame() %>%
  setNames(c('Model Type', 'Parameter', 'Estimate', 'Std. Error', 'Statistic',
             'p-value', 'Group')) 


AppendixModelTable$`p-value` <- lapply(AppendixModelTable$`p-value`,
                                       FUN = function (x) add_stars(x)) %>%
  unlist()

write.csv(AppendixModelTable,
          file = '../Eco_Letters_Manuscript/Figures/All_logistic_regression_outputs_for_Appendix.csv',
          na = "", row.names = FALSE)
