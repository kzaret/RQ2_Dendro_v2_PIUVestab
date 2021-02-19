### Coring height corrections:  
# regressions of ring counts on height along stem for harvested cross section.  
# To be used to estimate rings not captured due to coring height of cores.

#===========================================================================
# SETUP
#===========================================================================

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

#===========================================================================
# DATA
#===========================================================================

## Harvested PIUV from all patches
harv <- read.csv(here("data","RC_Height_cross_sections.csv"), header = TRUE)

# Notes:  
# harv includes dummy data such that each tree has a basal CS (height = 0) 
# and add. rings = 0 for each basal CS.
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

# Change nonzero ring counts at Height_RC==0 to zeros (affects one row)
# Add first-differenced section height and ring counts (for V.2 of GLMM)
harv <- harv %>% group_by(Sapling) %>% 
  mutate(AddRings = replace(AddRings, Height_RC==0 & AddRings > 0, 0),
         DiffHeight = Height_RC - lag(Height_RC),
         DiffRings = pmax(lag(RC) - RC, 0)) %>%  # change 3 small negative values to 0
  ungroup() %>% as.data.frame()

## Cored PIUV from all patches
cores <- read.csv(here("data","PIUV_CoredProcessed.csv"), header = TRUE)

# select & mutate predictor columns to match those used in the model
cores_nd <- cores %>%
  filter(Patch2 != "Cushion") %>% # remove cores from Cushion patch
  select(Patch2, Individual, Cor_Height_cm) %>%
  rename(Patch = Patch2, Sapling = Individual, Height_RC = Cor_Height_cm)

#===========================================================================
# POISSON GLMM
#===========================================================================

#---------------------------------------------------------------------------
# VERSION 1 
# Model additional rings relative to root-shoot boundary
# as a linear (proportional) function of section height
#---------------------------------------------------------------------------

# The poisson family function defaults to using the log link but we need the
# identity link to preserve the linear relationship b/w height and ring count
help(priors, package = 'rstanarm')  
vignette("priors", package = 'rstanarm')

# Model includes intercept for fixed effects, 
# but no intercept for random effect of individual harvested tree.  
harv_glmer1 <- stan_glmer(AddRings ~ Height_RC + Height_RC:Patch + (-1 + Height_RC | Sapling), 
                          data = harv, family = poisson(link = "identity"), 
                          chains = getOption("mc.cores"), iter = 3000, warmup = 1000) 

prior_summary(harv_glmer1)
print(harv_glmer1, digits=3)
summary(harv_glmer1)
cbind(rstan::get_elapsed_time(harv_glmer1$stanfit), 
      total = rowSums(rstan::get_elapsed_time(harv_glmer1$stanfit)))

# Marginal posterior predictive density
yrep <- posterior_predict(harv_glmer1)
indx <- sample(nrow(yrep), 100)
ppc_dens_overlay(harv$AddRings, yrep[indx,])

# Marginal posterior predictive density grouped by patch
ppc_dens_overlay_grouped(harv$AddRings, yrep[indx,], group = harv$Patch)

# Fitted vs. observed, grouped by patch
ppc_scatter_avg_grouped(harv$AddRings, yrep[indx,], group = harv$Patch) +
  geom_abline(intercept = 0, slope = 1)

# Normal QQ plot of tree-level random slope point estimates, grouped by patch
grp_slope <- as.matrix(harv_glmer1, regex_pars = "b")

colMedians(grp_slope) %>% data.frame() %>% setNames("slope") %>% 
  mutate(Patch = factor(tapply(harv$Patch, harv$Sapling, unique))) %>% 
  ggplot(aes(sample = slope)) + stat_qq(size = 2) + geom_qq_line() +
  theme_bw() + facet_wrap(vars(Patch), ncol = 2)

# Linear predictor (expectation) and posterior predictive density as fn of height and patch
# Use "new" level of grouping factor so PPD marginalizes over tree-level variance
# (predictions are for a "random tree")
# posterior_linpred() and posterior_predict() throw an lme4-derived error if predicted
# values go negative, so we have to add the effects of the tree-level random slopes by hand
fitdata <- expand.grid(AddRings = 0, Height_RC = 1:round(max(harv$Height_RC),-1), 
                       Patch = unique(harv$Patch), Sapling = "0") 

# linear predictor ignoring tree-level slopes (hyper-means only)
fit_linpred0 <- posterior_linpred(harv_glmer1, newdata = fitdata, re.form = NA)
# posterior draws of hyper-SD
sigma_slope <- as.matrix(harv_glmer1, regex_pars = "Sigma")
# posterior draws of a random slope from the hyperdistribution
bnew <- rnorm(nrow(sigma_slope), 0, sigma_slope)
# add tree-level effects to hyper-mean linear predictor
fit_linpred <- fit_linpred0 + tcrossprod(bnew, fitdata$Height_RC)
# truncate negative predictions (approx 0.01% of total draws)
fit_linpred <- pmax(fit_linpred, 0)
# posterior median and credible interval of linear predictor
fit_linpred_stats <- colQuantiles(fit_linpred, probs = c(0.025, 0.5, 0.975)) %>% 
  as.data.frame() %>% rename(c(lo = `2.5%`, med = `50%`, up = `97.5%`))

# posterior predictive distribution
fit_ppd <- matrix(rpois(length(fit_linpred), fit_linpred), nrow = nrow(fit_linpred))
# posterior median and credible interval of posterior predictive distribution
fit_ppd_stats <- colQuantiles(fit_ppd, probs = c(0.025, 0.5, 0.975)) %>% 
  as.data.frame() %>% rename(c(lo = `2.5%`, med = `50%`, up = `97.5%`))

## FIGURES

# Additional rings vs. height, grouped by patch
# Overlay posterior distribution (median and 95% credible interval) of linear predictor
save_plot <- TRUE
if(save_plot) {
  png(filename=here("analysis", "results", "harv_GLMM_fits.png"),
      width=7, height=7, units="in", res=300, type="cairo-png")
} else { dev.new() }

harv %>% ggplot(aes(Height_RC, AddRings, group = Sapling)) +
  geom_ribbon(aes(ymin = lo, ymax = up), data = cbind(fitdata, fit_linpred_stats),
              fill = "gray", alpha = 0.9) +
  geom_ribbon(aes(ymin = lo, ymax = up), data = cbind(fitdata, fit_ppd_stats),
              fill = "gray", alpha = 0.5) +
  geom_line(aes(Height_RC, med), data = cbind(fitdata, fit_linpred_stats),
            color = "darkgray", lwd = 1) +
  geom_line(alpha = 0.4) + geom_point(shape = 1, alpha = 0.5, size = 2) + 
  xlab("Height on stem (cm)") + ylab("Additional rings") + 
  theme_bw() + facet_wrap(vars(Patch), ncol = 2)

if(save_plot) dev.off()


#---------------------------------------------------------------------------
# VERSION 2 
# Model the difference in ring count from one section to the next
# as a mean (intercept) only, scaling for the difference in height
# as a multiplicative offset (additive log-offset)
# Implies the same proportional relationship b/w section height
# and total additional rings as V.1: the mean DiffRing / cm is equivalent
# to the slope w.r.t. Height_RC
# But incremental ring counts are a priori independent; 
# also we can use the log link (and lognormal random effects)
#---------------------------------------------------------------------------

# Intercept grouped by sapling
harv_glmer2 <- stan_glmer(DiffRings ~ Patch + (1 | Sapling), offset = log(DiffHeight), 
                          family = poisson(link = "log"), data = harv, na.action = na.omit,  
                          chains = getOption("mc.cores"), iter = 2000, warmup = 1000) 

prior_summary(harv_glmer2)
print(harv_glmer2, digits=2)
summary(harv_glmer2)
cbind(rstan::get_elapsed_time(harv_glmer2$stanfit), 
      total = rowSums(rstan::get_elapsed_time(harv_glmer2$stanfit)))

# Marginal posterior predictive density
yrep <- posterior_predict(harv_glmer2)
indx <- sample(nrow(yrep), 100)
ppc_dens_overlay(as.vector(na.omit(harv$DiffRings)), yrep[indx,])

# Rootogram of marginal posterior predictive density
ppc_rootogram(as.vector(na.omit(harv$DiffRings)), yrep[indx,])

# Marginal posterior predictive density grouped by patch
ppc_dens_overlay_grouped(as.vector(na.omit(harv$DiffRings)), yrep[indx,], 
                         group = harv$Patch[!is.na(harv$DiffRings)])

# Normal QQ plot of tree-level random intercept point estimates, grouped by patch
grp_intercept <- as.matrix(harv_glmer2, regex_pars = "b")

colMedians(grp_intercept) %>% data.frame() %>% setNames("intercept") %>% 
  mutate(Patch = factor(tapply(harv$Patch, harv$Sapling, unique))) %>% 
  ggplot(aes(sample = intercept)) + stat_qq(size = 2) + geom_qq_line() +
  theme_bw() + facet_wrap(vars(Patch), ncol = 2)



#---------------------------------------------------------------------------
# Predict outcome (additional ring counts) for cored trees 
#---------------------------------------------------------------------------

ppd_cores <- posterior_predict(harv_glmer, newdata = cores_nd)


#===========================================================================
# ADDITIONAL FIGURES
#===========================================================================

# Additional rings vs. height, all data
plot(harv$Height_RC, harv$AddRings)

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



#===========================================================================
# ADDITIONAL CODE
#===========================================================================

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

