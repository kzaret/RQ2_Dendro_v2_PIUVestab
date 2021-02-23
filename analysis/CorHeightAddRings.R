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
# COUNT GLMMs
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

# Additional rings vs. height, grouped by patch
# Overlay posterior distribution (median and 95% credible interval) of linear predictor
save_plot <- TRUE
if(save_plot) {
  png(filename=here("analysis", "results", "harv_GLMM_fits.png"),
      width=7, height=7, units="in", res=300, type="cairo-png")
} else dev.new()

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

# Poisson
# Intercept grouped by sapling and plot
harv_glmer2_pois <- stan_glmer(DiffRings ~ Patch + (1 | Sapling) + (1 | Plot), 
                               offset = log(DiffHeight), family = poisson(link = "log"),
                               data = harv, na.action = na.omit,  
                               chains = getOption("mc.cores"), iter = 3000, warmup = 1000) 

prior_summary(harv_glmer2_pois)
print(harv_glmer2_pois, digits=2)
summary(harv_glmer2_pois)
cbind(rstan::get_elapsed_time(harv_glmer2_pois$stanfit), 
      total = rowSums(rstan::get_elapsed_time(harv_glmer2_pois$stanfit)))

# Negative binomial
# Intercept grouped by sapling and plot
harv_glmer2_nb <- stan_glmer(DiffRings ~ Patch + (1 | Sapling) + (1 | Plot), 
                             offset = log(DiffHeight), family = neg_binomial_2(link = "log"),
                             prior_aux = exponential(0.1),
                             data = harv, na.action = na.omit,  
                             chains = getOption("mc.cores"), iter = 3000, warmup = 1000) 

prior_summary(harv_glmer2_nb)
print(harv_glmer2_nb, digits=2)
summary(harv_glmer2_pois)
cbind(rstan::get_elapsed_time(harv_glmer2_nb$stanfit), 
      total = rowSums(rstan::get_elapsed_time(harv_glmer2_nb$stanfit)))

# Compare Poisson and negative binomial models using LOO 
# to estimate out-of-sample predictive performance
loo_pois <- loo(harv_glmer2_pois)
loo_pois

loo_nb <- loo(harv_glmer2_nb)
loo_nb

loo_compare(loo_pois, loo_nb)

## FIGURES

# Choose your fighter
mod <- harv_glmer2_nb

# Histograms of data and draws from marginal PPD
# Note that the offset argument is apparently needed, 
# contrary to help(posterior_predict)
yrep <- posterior_predict(mod, offset = mod$offset)
indx <- sample(nrow(yrep), 100)
ppc_hist(as.vector(na.omit(harv$DiffRings)), yrep[indx[1:3],]) + ggtitle(mod$family$family)

# Rootogram of marginal posterior predictive density
ppc_rootogram(as.vector(na.omit(harv$DiffRings)), yrep[indx,]) + ggtitle(mod$family$family)

# Posterior predictive check: mean
ppc_stat_grouped(as.vector(na.omit(harv$DiffRings)), yrep[indx,], 
                 group = harv$Patch[!is.na(harv$DiffRings)], stat = mean) + 
  ggtitle(mod$family$family)

# Posterior predictive check: SD
ppc_stat_grouped(as.vector(na.omit(harv$DiffRings)), yrep[indx,], 
                 group = harv$Patch[!is.na(harv$DiffRings)], stat = sd) + 
  ggtitle(mod$family$family)

# Posterior predictive check: dispersion
ppc_stat_grouped(as.vector(na.omit(harv$DiffRings)), yrep[indx,], 
                 group = harv$Patch[!is.na(harv$DiffRings)], 
                 stat = function(x) var(x)/mean(x)) + 
  guides(fill = guide_legend(title = "V(Y)/E(Y)", title.vjust = 5)) + 
  ggtitle(mod$family$family)

# Posterior predictive check: proportion of zeros
ppc_stat_grouped(as.vector(na.omit(harv$DiffRings)), yrep[indx,], 
                 group = harv$Patch[!is.na(harv$DiffRings)],
                 stat = function(x) mean(x == 0)) +
  guides(fill = guide_legend(title = "Proportion zeros", title.vjust = 5)) + 
  ggtitle(mod$family$family)

# Normal QQ plot of tree-level random intercept point estimates, grouped by patch
grp_intercept <- as.data.frame(mod, regex_pars = "Sapling:") %>% 
  select(-contains("Sigma")) %>% as.matrix()

colMedians(grp_intercept) %>% data.frame() %>% setNames("intercept") %>% 
  mutate(Patch = factor(tapply(harv$Patch, harv$Sapling, unique))) %>% 
  ggplot(aes(sample = intercept)) + stat_qq(size = 2) + geom_qq_line() +
  theme_bw() + facet_wrap(vars(Patch), ncol = 2) + ggtitle(mod$family$family)

# Normal QQ plot of plot-level random intercept point estimates, grouped by patch
grp_intercept <- as.data.frame(mod, regex_pars = "Plot:") %>% 
  select(-contains("Sigma")) %>% as.matrix()

colMedians(grp_intercept) %>% data.frame() %>% setNames("intercept") %>% 
  mutate(Patch = factor(tapply(harv$Patch, harv$Plot, unique))) %>% 
  ggplot(aes(sample = intercept)) + stat_qq(size = 2) + geom_qq_line() +
  theme_bw() + facet_wrap(vars(Patch), ncol = 2) + ggtitle(mod$family$family)

# Dotplots of rings/cm by sapling
# Overlay violin plots of PPD
save_plot <- TRUE
if(save_plot) {
  png(filename=here("analysis", "results", paste0("harv_GLMMv2_", mod$family$family, "_fits.png")),
      width=7, height=7, units="in", res=300, type="cairo-png")
} else dev.new()

dat <- data.frame(iter = rep(indx, ncol(yrep)),
                  Patch = rep(mod$glmod$fr$Patch, each = length(indx)),
                  Sapling = rep(mod$glmod$fr$Sapling, each = length(indx)),
                  DiffHeight = rep(exp(mod$offset), each = length(indx)),
                  DiffRings = as.vector(yrep[indx,]))

harv %>% group_by(Sapling) %>% 
  ungroup() %>% as.data.frame() %>%
  ggplot(aes(x = Sapling, y = DiffRings / DiffHeight)) +
  geom_violin(aes(x = Sapling, y = DiffRings / DiffHeight, group = Sapling), 
              data = dat, color = "darkgray", fill = "darkgray", alpha = 0.8) +
  geom_jitter(shape = 16, alpha = 0.5, size = 1, width = 0.1, height = 0) +
  xlab("Sapling") + ylab("Additional rings / cm") +
  theme_bw() + theme(axis.text.x = element_blank()) + 
  facet_wrap(vars(Patch), scales = "free_x") + 
  ggtitle(mod$family$family)

if(save_plot) dev.off()

# Scatterplot of rings/cm vs section height
# Any residual relationship?
harv %>% arrange(Patch, Plot, Sapling, Height_RC) %>% group_by(Sapling) %>% 
  ungroup() %>% as.data.frame() %>%
  ggplot(aes(x = Height_RC, y = DiffRings / DiffHeight, group = Sapling)) +
  geom_point(shape = 1, alpha = 0.5) + geom_line(alpha = 0.4) +
  xlab("Section height") + ylab("Additional rings / cm") +
  theme_bw() + facet_wrap(vars(Patch), scales = "free") + 
  ggtitle(mod$family$family)


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

