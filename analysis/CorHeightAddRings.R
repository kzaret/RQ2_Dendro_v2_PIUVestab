### Coring height corrections:  
# regressions of ring counts on height along stem for harvested cross section.  
# To be used to estimate rings not captured due to coring height of cores.

library(mlmRev) # for mlm tutorial (https://mc-stan.org/users/documentation/case-studies/tutorial_rstanarm.html)
library(lme4)
library(rstanarm)
library(dplyr)
library(bayesplot)
library(here)
bayesplot_theme_set(bayesplot::theme_default())
options(mc.cores = parallel::detectCores() - 1)


#############################################################################

## Poisson GLMM

# Data:  
# harvested PIUV from all patches
# pot_outlier = 1 if had been removed when fitting linear regressions and/or 
# given abnormal ring counts in consecutive rounds.

# Re Bog forest outliers:  
# For linear regressions, had removed BF0802 entirely [9 cross sections] due to 
# nonlinear growth pattern that differs from the rest of the harvested trees.  
# It’s also one of the older BF trees. 
# Removed cross-sections 06, 05 and 04 from BF0902 
# (bc ring counts of 92, 91, 90 and 90 in a row and notes re 04 being rotten 
# and may want to leave out). 
# Removed ff bc high on stem and no cross sections in between:  
# BF0303_110to112, BF0305_100to102. Removed BF0202 [5 cross sections]. 
# Removed BF0301_110to112.

harv <- read.csv(here("data","RC_Height_cross_sections.csv"), header = TRUE)

colnames(harv)
head(harv)

# The poisson family function defaults to using the log link
help(priors, package = 'rstanarm')  
vignette("priors", package = 'rstanarm')

# Model includes intercept for fixed effects, 
# but no intercept for random effect of individual harvested tree.  
# Note:  could add dummy data such that each tree has a basal CS (height = 0) 
# and add. rings = 0 for each basal CS.
harv_glmer <- stan_glmer(AddRings ~ Height_RC + Height_RC:Patch + (-1 + Height_RC | Sapling), 
                         data = harv, family = poisson(link = "identity"), 
                         iter = 2000) 

# errors @ inter = default of 2000: 
# Warning messages: 
#1: Bulk Effective Samples Size (ESS) is too low, 
# indicating posterior means and medians may be unreliable.
# Running the chains for more iterations may help. See http://mc-stan.org/misc/warnings.html#bulk-ess 
# 2: Tail Effective Samples Size (ESS) is too low, indicating posterior variances 
# and tail quantiles may be unreliable. Running the chains for more iterations may help. 
# See http://mc-stan.org/misc/warnings.html#tail-ess

# Time for iter = 3000; 18:08 - 18:49; only got the 1st of the above errors; 
# summary is the same; removing outliers didn't change the results.

prior_summary(harv_glmer)
print(harv_glmer, digits=3)

yrep_harv <- posterior_predict(harv_glmer)
n_sims <- nrow(yrep_harv)
subset <-sample(n_sims, 100)
ppc_dens_overlay(harv$AddRings, yrep_harv)

?posterior_predict
?ppc_dens_overlay


# Predict outcome (additional ring counts) for new values of predictors 
# (height, patch, sapling)

cores <- read.csv(here("data","PIUV_CoredProcessed.csv"), header = TRUE)

#select & mutate predictor columns to match those used in the model
cores.nd <- cores %>%
  filter(Patch2 != "Cushion") %>% #remove cores from Cushion patch
  select(Patch2, Individual, Cor_Height_cm) %>%
  mutate(Patch = Patch2, Sapling = Individual, Height_RC = Cor_Height_cm, .keep = "unused")

colnames(cores.nd)
dim(cores.nd)

ynew <- posterior_predict (harv_glmer, newdata=cores.nd)



###############################################################################

###plots

plot(harv$Height_RC, harv$AddRings)

harv %>% ggplot(aes(Height_RC, AddRings)) +
  geom_point() + 
  facet_grid(Patch~.) + theme_bw() +
  xlab("Height on Stem (cm)") + ylab("No. Additional Rings")


#Scatterplots AddRings ~ Height on Stem, patch, symbolized by sapling
harv %>% filter(Patch=="Forest") %>%
ggplot(aes(Height_RC, AddRings)) + geom_point(aes(color=Sapling)) +
  xlab("Height on Stem (cm)") + ylab("No. Add. Rings to Base")

#Scatterplots AddRings ~ Height on Stem, plot, symbolized by sapling
harv %>% filter(Plot=="BF03") %>%
  ggplot(aes(Height_RC, AddRings)) + geom_point(aes(color=Sapling)) +
  xlab("Height on Stem (cm)") + ylab("No. Add. Rings to Base")

#log(AddRings)
harv_log <- harv %>% filter(AddRings>0) %>%
  mutate(AddRings_log = log(AddRings))

head(harv_log)

harv_log %>% filter(Plot=="F03") %>%
  ggplot(aes(Height_RC, AddRings_log)) + geom_point(aes(color=Sapling)) +
  xlab("Height on Stem (cm)") + ylab("No. Add. Rings to Base")

################################################################################

###Poisson GLMs

#link = log
poi1 <- stan_glm(AddRings ~ Height_RC, family = poisson(link="log"), data = harv)
print(poi1, digit=3)

plot(harv$Height_RC, harv$AddRings)
curve(exp(coef(poi1)[1] + coef(poi1)[2]*x), add=TRUE)

yrep_1 <- posterior_predict(poi1)
n_sims <- nrow(yrep_1)
subset <-sample(n_sims, 100)
ppc_dens_overlay(harv$AddRings, yrep_1)

#link = identity
poi2 <- stan_glm(AddRings ~ Height_RC, family = poisson(link="identity"), data = harv)
print(poi2, digit=3)

plot(harv$Height_RC, harv$AddRings)
curve(coef(poi2)[1] + coef(poi2)[2]*x, add=TRUE)

yrep_2 <- posterior_predict(poi2)
n_sims <- nrow(yrep_2)
subset <-sample(n_sims, 100)
ppc_dens_overlay(harv$AddRings, yrep_2)


#link = identity, patch as random effect
poi3 <- stan_glm(AddRings ~ Height_RC + Height_RC:Patch, family = poisson(link="identity"), data = harv)
print(poi3, digit=3)

plot(harv$Height_RC, harv$AddRings)
curve(coef(poi2)[1] + coef(poi2)[2]*x, add=TRUE)

yrep_2 <- posterior_predict(poi2)
n_sims <- nrow(yrep_2)
subset <-sample(n_sims, 100)
ppc_dens_overlay(harv$AddRings, yrep_2)

