### Coring height corrections:  
# regressions of ring counts on height along stem for harvested cross section.  
# To be used to estimate rings not captured due to coring height of cores.

#---------------------------------------------------------------------------
# SETUP
#---------------------------------------------------------------------------

library(mlmRev) # for mlm tutorial (https://mc-stan.org/users/documentation/case-studies/tutorial_rstanarm.html)
library(lme4)
library(rstanarm)
library(dplyr)
library(bayesplot)
library(shinystan)
library(ggplot2)
library(matrixStats)
library(here)
bayesplot_theme_set(bayesplot::theme_default())
options(mc.cores = parallel::detectCores(logical = FALSE) - 1)
if(.Platform$OS.type == "windows") options(device = windows)


#---------------------------------------------------------------------------
# POISSON GLMM
#---------------------------------------------------------------------------

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
                         chains = getOption("mc.cores"), iter = 3000, warmup = 1000) 

prior_summary(harv_glmer)
print(harv_glmer, digits=3)
summary(harv_glmer)
cbind(rstan::get_elapsed_time(harv_glmer$stanfit), 
      total = rowSums(rstan::get_elapsed_time(harv_glmer$stanfit)))

# Marginal posterior predictive density
yrep_harv <- posterior_predict(harv_glmer)
indx <- sample(nrow(yrep_harv), 100)
ppc_dens_overlay(harv$AddRings, yrep_harv[indx,])

# Marginal posterior predictive density grouped by patch
ppc_dens_overlay_grouped(harv$AddRings, yrep_harv[indx,], group = harv$Patch)

# Normal QQ plot of tree-level random slope point estimates, grouped by patch
grp_slopes <- as.matrix(harv_glmer, regex_pars = "b")

colMedians(grp_slopes) %>% data.frame() %>% setNames("slope") %>% 
  mutate(Patch = factor(tapply(harv$Patch, harv$Sapling, unique))) %>% 
  ggplot(aes(sample = slope)) + stat_qq(size = 2) + geom_qq_line() +
  theme_bw() + facet_wrap(vars(Patch), ncol = 2)

# Linear predictor (expectation) and posterior predictive density 
# as a function of height and patch
# Use "new" level of grouping factor so PPD marginalizes over tree-level variance
# (predictions are for a "random tree")
fitdata <- expand.grid(AddRings = 0, Height_RC = 1:round(max(harv$Height_RC),-1), 
                       Patch = unique(harv$Patch), Sapling = "0")  
fit_linpred <- posterior_linpred(harv_glmer, newdata = fitdata, re.form = NA)
fit_linpred_stats <- colQuantiles(fit_linpred, probs = c(0.025, 0.5, 0.975)) %>% 
  as.data.frame() %>% rename(c(lo = `2.5%`,med = `50%`,up = `97.5%`))
  
fit_ppd <- posterior_predict(harv_glmer, newdata = fitdata, re.form = NA)

#---------------------------------------------------------------------------
# Predict outcome (additional ring counts) for cored trees 
#---------------------------------------------------------------------------

cores <- read.csv(here("data","PIUV_CoredProcessed.csv"), header = TRUE)

# select & mutate predictor columns to match those used in the model
cores_nd <- cores %>%
  filter(Patch2 != "Cushion") %>% # remove cores from Cushion patch
  select(Patch2, Individual, Cor_Height_cm) %>%
  mutate(Patch = Patch2, Sapling = Individual, Height_RC = Cor_Height_cm, .keep = "unused")

colnames(cores_nd)
dim(cores_nd)

ppd_cores <- posterior_predict(harv_glmer, newdata = cores_nd)


#---------------------------------------------------------------------------
# FIGURES
#---------------------------------------------------------------------------

# Additional rings vs. height, all data
plot(harv$Height_RC, harv$AddRings)

# Additional rings vs. height, grouped by patch
# Overlay posterior distribution (median and 95% credible interval) of linear predictor
save_plot <- TRUE
if(save_plot) {
  png(filename=here("analysis", "results", "harv_GLMM_fits.png"),
                  width=7, height=7, units="in", res=300, type="cairo-png")
} else { dev.new(width=13,height=8.5) }

harv %>% ggplot(aes(Height_RC, AddRings)) +
  geom_ribbon(aes(ymin = lo, ymax = up), data = cbind(fitdata, fit_linpred_stats),
              fill = "lightgray", alpha = 0.7) +
  geom_line(aes(Height_RC, med), data = cbind(fitdata, fit_linpred_stats),
            color = "darkgray", lwd = 1) +
  geom_point(shape = 16, alpha = 0.6, size = 2) + 
  xlab("Height on stem (cm)") + ylab("Additional rings") + 
  theme_bw() + facet_wrap(vars(Patch), ncol = 2)

if(save_plot) dev.off()

# Scatterplots AddRings ~ Height on Stem, grouped by patch, symbolized by sapling
harv %>% ggplot(aes(Height_RC, AddRings)) + geom_point(aes(color=Sapling)) +
  xlab("Height on Stem (cm)") + ylab("No. Add. Rings to Base") + 
  facet_wrap(vars(Patch), ncol = 1) + theme_bw()

# Scatterplots AddRings ~ Height on Stem, grouped by plot, symbolized by sapling
harv %>% filter(Plot=="BF03") %>%
  ggplot(aes(Height_RC, AddRings)) + geom_point(aes(color=Sapling)) +
  xlab("Height on Stem (cm)") + ylab("No. Add. Rings to Base")

# log(AddRings)
harv_log <- harv %>% filter(AddRings>0) %>%
  mutate(AddRings_log = log(AddRings))

head(harv_log)

harv_log %>% filter(Plot=="F03") %>%
  ggplot(aes(Height_RC, AddRings_log)) + geom_point(aes(color=Sapling)) +
  xlab("Height on Stem (cm)") + ylab("No. Add. Rings to Base")

################################################################################

### Poisson GLMs

# link = log
poi1 <- stan_glm(AddRings ~ Height_RC, family = poisson(link="log"), data = harv)
print(poi1, digit=3)

plot(harv$Height_RC, harv$AddRings)
curve(exp(coef(poi1)[1] + coef(poi1)[2]*x), add=TRUE)

yrep_1 <- posterior_predict(poi1)
n_sims <- nrow(yrep_1)
subset <-sample(n_sims, 100)
ppc_dens_overlay(harv$AddRings, yrep_1)

# link = identity
poi2 <- stan_glm(AddRings ~ Height_RC, family = poisson(link="identity"), data = harv)
print(poi2, digit=3)

plot(harv$Height_RC, harv$AddRings)
curve(coef(poi2)[1] + coef(poi2)[2]*x, add=TRUE)

yrep_2 <- posterior_predict(poi2)
n_sims <- nrow(yrep_2)
subset <-sample(n_sims, 100)
ppc_dens_overlay(harv$AddRings, yrep_2)


# link = identity, patch as fixed effect
poi3 <- stan_glm(AddRings ~ Height_RC + Height_RC:Patch, family = poisson(link="identity"), data = harv)
print(poi3, digit=3)

plot(harv$Height_RC, harv$AddRings)
curve(coef(poi2)[1] + coef(poi2)[2]*x, add=TRUE)

yrep_2 <- posterior_predict(poi2)
n_sims <- nrow(yrep_2)
subset <-sample(n_sims, 100)
ppc_dens_overlay(harv$AddRings, yrep_2)

