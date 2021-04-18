### Coring height corrections:  
# regressions of ring counts on height along stem for harvested cross section.  
# To be used to estimate rings not captured due to coring height of cores.

#===========================================================================
# SETUP
#===========================================================================

library(mlmRev) # for mlm tutorial (https://mc-stan.org/users/documentation/case-studies/tutorial_rstanarm.html)
library(matrixStats)
library(dplyr)
library(tidyr)
library(forcats)
library(rstanarm)
library(bayesplot)
library(shinystan)
library(truncnorm)
library(ggplot2)
library(ggridges)
library(here)
options(mc.cores = parallel::detectCores(logical = FALSE) - 1)
if(.Platform$OS.type == "windows") options(device = windows)

#===========================================================================
# DATA
#===========================================================================

## Harvested PIUV from all patches
# includes dummy data such that each tree has a basal CS (height = 0) 
# and add. rings = 0 for each basal CS.
harv_raw <- read.csv(here("data","RC_Height_cross_sections.csv"), header = TRUE)

# Change nonzero ring counts at height==0 to zeros (affects one row)
# Add first-differenced section height and ring counts (for V.2 of GLMM)
harv <- harv_raw %>% group_by(Sapling) %>% 
  rename(harvest_year = Harv_Yr, site = Site, order = Order, patch = Patch, plot = Plot,
         tree = Sapling, round = Round, height = Height_RC, rings = RC, establishment_year = Estab_Yr,
         priority = Priority, add_rings = AddRings, pot_outlier = Pot_outlier) %>% 
  mutate(add_rings = replace(add_rings, height==0 & add_rings > 0, 0),
         diff_height = height - lag(height),
         diff_rings = pmax(lag(rings) - rings, 0)) %>%  # change 3 small negative values to 0
  ungroup() %>% as.data.frame()

## Cored PIUV from all patches
cores_raw <- read.csv(here("data","PIUV_CoredProcessed.csv"), header = TRUE)

# select & mutate predictor columns to match those used in the model
# compute height measurement error SD
# remove Cushion patch, R0X and NA plots, and cores w/ remove == 1 (nonrandomly sampled or damaged) 
cores <- cores_raw %>%
  select(Patch2, Plot, Year_coll, Individual, Cor_Height_cm, Status, Outer_rings, 
         ORW, Remove, min_age, TR_Count, height_error_u95_cm) %>%
  rename(patch = Patch2, plot = Plot, year = Year_coll, tree = Individual, 
         height = Cor_Height_cm, status = Status, outer_rings = Outer_rings, orw = ORW,
         remove = Remove, ring_count = TR_Count, height_error_u95 = height_error_u95_cm) %>% 
  mutate(height_SD = height_error_u95/qnorm(0.975)) %>% 
  filter(patch != "Cushion" & plot != "R0X" & !is.na(plot) & remove == 0) 

#===========================================================================
# COUNT GLMMs
#
# Model change in ring count as linear (proportional) function of height 
# along stem
# On log link scale, this is equivalent to an intercept plus log(height)
# as an offset
# Models include fixed effect of patch plus tree- and plot-level 
# random effects on intercept
#===========================================================================

# Info about rstanarm priors
help(priors, package = 'rstanarm')  
vignette("priors", package = 'rstanarm')

#---------------------------------------------------------------------------
# VERSION 1 
# Model total additional rings relative to root-shoot boundary
#
# NOTE: STATISTICALLY INVALID because cumulative ring counts within a tree
# are mathematically nonindependent. Use incremental counts (V2) instead.
#
# Exclude "artificial zeros" at zero height b/c this is not a valid
# log-offset, and zero origin is already implied by proportionality
#---------------------------------------------------------------------------

# Poisson likelihood (log link is default)
harv_glmer1 <- stan_glmer(add_rings ~ patch + (1 | tree) + (1 | plot), 
                          offset = log(height), family = poisson, 
                          data = harv, subset = height > 0,
                          chains = getOption("mc.cores"), iter = 5000, warmup = 1000) 

prior_summary(harv_glmer1)
print(harv_glmer1, digits=2)
summary(harv_glmer1)
cbind(rstan::get_elapsed_time(harv_glmer1$stanfit), 
      total = rowSums(rstan::get_elapsed_time(harv_glmer1$stanfit)))

## FIGURES ##

# Simulate draws from posterior predictive distribution
# Note that the offset argument is apparently needed, 
# contrary to help(posterior_predict)
yrep <- posterior_predict(harv_glmer1, offset = harv_glmer1$offset)
indx <- sample(nrow(yrep), 1000)

# Histograms of data and draws from marginal PPD
ppc_hist(harv_glmer1$y, yrep[indx[1:3],]) + ggtitle(harv_glmer1$family$family)

# Rootogram of marginal posterior predictive density
ppc_rootogram(harv_glmer1$y, yrep[indx,]) + ggtitle(harv_glmer1$family$family)

# Posterior predictive check: mean
ppc_stat_grouped(harv_glmer1$y, yrep[indx,], group = harv_glmer1$glmod$fr$patch, stat = mean) + 
  ggtitle(harv_glmer1$family$family)

# Posterior predictive check: SD
ppc_stat_grouped(harv_glmer1$y, yrep[indx,], group = harv_glmer1$glmod$fr$patch, stat = sd) + 
  ggtitle(harv_glmer1$family$family)

# Posterior predictive check: dispersion
ppc_stat_grouped(harv_glmer1$y, yrep[indx,], group = harv_glmer1$glmod$fr$patch, 
                  stat = function(x) var(x)/mean(x)) + 
  guides(fill = guide_legend(title = "V(Y)/E(Y)", title.vjust = 5)) + 
  ggtitle(harv_glmer1$family$family)

# Posterior predictive check: proportion of zeros
ppc_stat_grouped(harv_glmer1$y, yrep[indx,], group = harv_glmer1$glmod$fr$patch,
                 stat = function(x) mean(x == 0)) +
  guides(fill = guide_legend(title = "Proportion zeros", title.vjust = 5)) + 
  ggtitle(harv_glmer1$family$family)

# Normal QQ plot of tree-level random intercept point estimates, grouped by patch
ranef(harv_glmer1)$tree %>% rename(intercept = `(Intercept)`) %>% 
  mutate(patch = factor(tapply(harv$patch, harv$tree, unique))) %>% 
  ggplot(aes(sample = intercept)) + stat_qq(size = 2) + geom_qq_line() +
  theme_bw() + theme(panel.grid = element_blank()) +
  facet_wrap(vars(patch), ncol = 2) + ggtitle(harv_glmer1$family$family)

# Normal QQ plot of plot-level random intercept point estimates, grouped by patch
ranef(harv_glmer1)$plot %>% rename(intercept = `(Intercept)`) %>% 
  mutate(patch = factor(tapply(harv$patch, harv$plot, unique))) %>% 
  ggplot(aes(sample = intercept)) + stat_qq(size = 2) + geom_qq_line() +
  theme_bw() +  theme(panel.grid = element_blank()) +
  facet_wrap(vars(patch), ncol = 2) + ggtitle(harv_glmer1$family$family)


#---------------------------------------------------------------------------
# VERSION 2 
# Model the difference in ring count from one section to the next
# Unlike total additional rings, incremental ring counts are 
# a priori independent
#---------------------------------------------------------------------------

# Poisson likelihood (log link is default)
harv_glmer2_pois <- stan_glmer(diff_rings ~ patch + (1 | tree) + (1 | plot), 
                               offset = log(diff_height), family = poisson,
                               data = harv, na.action = na.omit,  
                               chains = getOption("mc.cores"), iter = 5000, warmup = 1000) 

prior_summary(harv_glmer2_pois)
print(harv_glmer2_pois, digits=2)
summary(harv_glmer2_pois)
cbind(rstan::get_elapsed_time(harv_glmer2_pois$stanfit), 
      total = rowSums(rstan::get_elapsed_time(harv_glmer2_pois$stanfit)))

# Negative binomial likelihood (log link is default)
harv_glmer2_nb <- stan_glmer(diff_rings ~ patch + (1 | tree) + (1 | plot), 
                             offset = log(diff_height), family = neg_binomial_2,
                             prior_aux = exponential(0.1),
                             data = harv, na.action = na.omit,  
                             chains = getOption("mc.cores"), iter = 5000, warmup = 1000) 

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

## FIGURES ## 

# Choose your fighter
mod <- harv_glmer2_pois

# Simulate draws from posterior predictive distribution
# Note that the offset argument is apparently needed, 
# contrary to help(posterior_predict)
yrep <- posterior_predict(mod, offset = mod$offset)
indx <- sample(nrow(yrep), 1000)

# Histograms of data and draws from marginal PPD
ppc_hist(as.vector(na.omit(harv$diff_rings)), yrep[indx[1:3],]) + ggtitle(mod$family$family)

# Rootogram of marginal posterior predictive density
ppc_rootogram(as.vector(na.omit(harv$diff_rings)), yrep[indx,]) + ggtitle(mod$family$family)

# Posterior predictive check: mean
ppc_stat_grouped(mod$y, yrep[indx,], group = mod$glmod$fr$patch, stat = mean) + 
  ggtitle(mod$family$family)

# Posterior predictive check: SD
ppc_stat_grouped(mod$y, yrep[indx,], group = mod$glmod$fr$patch, stat = sd) + 
  ggtitle(mod$family$family)

# Posterior predictive check: dispersion
ppc_stat_grouped(mod$y, yrep[indx,], group = mod$glmod$fr$patch, stat = function(x) var(x)/mean(x)) + 
  guides(fill = guide_legend(title = "V(Y)/E(Y)", title.vjust = 5)) + ggtitle(mod$family$family)

# Posterior predictive check: proportion of zeros
ppc_stat_grouped(mod$y, yrep[indx,], group = mod$glmod$fr$patch, stat = function(x) mean(x == 0)) +
  guides(fill = guide_legend(title = "Proportion zeros", title.vjust = 5)) + ggtitle(mod$family$family)

# Normal QQ plot of tree-level random intercept point estimates, grouped by patch
ranef(mod)$tree %>% rename(intercept = `(Intercept)`) %>% 
  mutate(patch = factor(tapply(harv$patch, harv$tree, unique))) %>% 
  ggplot(aes(sample = intercept)) + stat_qq(size = 2) + geom_qq_line() +
  theme_bw() +  theme(panel.grid = element_blank()) +
  facet_wrap(vars(patch), ncol = 2) + ggtitle(mod$family$family)

# Normal QQ plot of plot-level random intercept point estimates, grouped by patch
ranef(mod)$plot %>% rename(intercept = `(Intercept)`) %>% 
  mutate(patch = factor(tapply(harv$patch, harv$plot, unique))) %>% 
  ggplot(aes(sample = intercept)) + stat_qq(size = 2) + geom_qq_line() +
  theme_bw() +  theme(panel.grid = element_blank()) +
  facet_wrap(vars(patch), ncol = 2) + ggtitle(mod$family$family)

# Dotplots of rings/cm by sapling
# Overlay violin plots of PPD
save_plot <- TRUE
dev.new()

dat <- data.frame(iter = rep(indx, ncol(yrep)),
                  patch = rep(mod$glmod$fr$patch, each = length(indx)),
                  tree = rep(mod$glmod$fr$tree, each = length(indx)),
                  diff_height = rep(exp(mod$offset), each = length(indx)),
                  diff_rings = as.vector(yrep[indx,])) %>% 
  mutate(tree = fct_reorder(.f = tree, .x = diff_rings / diff_height, .fun = mean))

harv %>% ggplot(aes(x = tree, y = diff_rings / diff_height)) +
  geom_violin(aes(x = tree, y = diff_rings / diff_height), 
              data = dat, color = "darkgray", fill = "darkgray", alpha = 0.8) +
  geom_jitter(shape = 16, alpha = 0.5, size = 1, width = 0.1, height = 0) +
  xlab("Sapling") + ylab("Additional rings / cm") +
  theme_bw() + theme(axis.text.x = element_blank()) + 
  facet_wrap(vars(patch), scales = "free_x") + 
  ggtitle(mod$family$family)

if(save_plot) 
  ggsave(filename=here("analysis", "results", paste0("harv_GLMMv2_", mod$family$family, "_fits.png")),
         width=7, height=7, units="in", dpi=300, type="cairo-png")

# Extrapolate model fitted to incremental ring counts to predict total additional rings
# from root-shoot boundary

# Inverse-link transformed linear predictor (expectation) and PPD as fn of height and patch
# Use "new" level of grouping factors to marginalize over plot- and tree-level variance
# (predictions are for a "random tree in a random plot")
predata <- expand.grid(diff_rings = 0, diff_height = 1:round(max(harv$height),-1),
                       patch = unique(harv$patch), plot = "0", tree = "0")
# transformed linear predictor
pre_epred <- posterior_epred(mod, newdata = predata, offset = log(predata$diff_height))
# posterior median and credible interval of transformed linear predictor
pre_epred_stats <- colQuantiles(pre_epred, probs = c(0.025, 0.5, 0.975)) %>%
  as.data.frame() %>% rename(c(lo = `2.5%`, med = `50%`, up = `97.5%`))
# posterior predictive distribution
pre_ppd <- posterior_predict(mod, newdata = predata, offset = log(predata$diff_height))
# posterior median and credible interval of posterior predictive distribution
pre_ppd_stats <- colQuantiles(pre_ppd, probs = c(0.025, 0.5, 0.975)) %>%
  as.data.frame() %>% rename(c(lo = `2.5%`, med = `50%`, up = `97.5%`))

# Cumulative additional rings vs. height, grouped by patch
# Overlay posterior distribution (median and 95% credible interval) of expectation and PPD
save_plot <- TRUE
dev.new()

harv %>% ggplot(aes(height, add_rings, group = tree)) +
  geom_ribbon(aes(x = diff_height, ymin = lo, ymax = up), inherit.aes = FALSE,
              data = cbind(predata, pre_epred_stats), fill = "gray", alpha = 0.9) +
  geom_ribbon(aes(x = diff_height, ymin = lo, ymax = up), inherit.aes = FALSE,
              data = cbind(predata, pre_ppd_stats), fill = "gray", alpha = 0.5) +
  geom_line(aes(diff_height, med), data = cbind(predata, pre_epred_stats),
            color = "darkgray", lwd = 1) +
  geom_line(alpha = 0.4) + geom_point(shape = 1, alpha = 0.5, size = 2) + 
  xlab("Height on stem (cm)") + ylab("Additional rings") + 
  theme_bw() + facet_wrap(vars(patch), ncol = 2) + 
  ggtitle(paste("GLMMv2", mod$family$family))

if(save_plot) 
  ggsave(filename=here("analysis", "results", paste0("harv_GLMMv2_", mod$family$family, "_predict_AddRings.png")),
         width=7, height=7, units="in", dpi=300, type="cairo-png")


#===========================================================================
# CORING-HEIGHT CORRECTIONS INCLUDING HEIGHT MEASUREMENT ERROR
#
# Generate coring height from zero-truncated normal measurement error prior
# => Unconditional posterior of expected additional ring counts
# (from root-shoot boundary) 
# Plots without harvested saplings have plot-level coefficients drawn
# from the hyperdistribution
#===========================================================================

# Choose your fighter
mod <- harv_glmer2_nb

# Generate prior draws of height with measurement error
# and posterior draws of expected additional rings / cm
# Adjust using measurement error prior as height offset
height_prior <- t(sapply(1:nsamples(mod), rtruncnorm, n = 1, mean = cores$height, 
                         sd = cores$height_SD, a = 0))
add_rings_cm <- posterior_epred(mod, newdata = cores, offset = 0)
add_rings <- add_rings_cm * height_prior

# Attach posterior draws to data
add_rings_height <- cores %>% cbind(as.data.frame(t(add_rings))) %>% 
  pivot_longer(cols = starts_with("V"), values_to = "add_rings_height", 
               names_to = "iter", names_prefix = "V", 
               names_transform = list(iter = as.numeric)) %>% as.data.frame()

# Save posterior draws and stanfit objects
save(list = c("harv_glmer2_pois", "harv_glmer2_nb", "add_rings_height"), 
     file = here("analysis","results","coring_height_correction.RData"))


#---------------------------------------------------------------------------
# FIGURES 
#---------------------------------------------------------------------------

# Credible intervals of posterior distribution of additional rings
save_plot <- TRUE
dev.new()

mcmc_intervals_data(add_rings, prob = 0.8, prob_outer = 0.95) %>% 
  mutate(patch = cores$patch, tree = fct_reorder(cores$tree, .x = m, .fun = identity)) %>% 
  ggplot(aes(x = tree, y = m)) +
  geom_linerange(aes(ymin = l, ymax = h), size = 1.5, color = "darkgray") +
  geom_linerange(aes(ymin = ll, ymax = hh), size = 0.5, color = "darkgray") +
  geom_point(pch = 16, size = 1) +
  xlab("Tree") + ylab("Additional rings at root-shoot boundary") + 
  theme_bw(base_size = 16) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        panel.grid = element_blank()) + 
  facet_wrap(vars(patch), scales = "free")

if(save_plot) 
  ggsave(filename=here("analysis", "results", "coring-height-correction_intervals.png"),
         width=7, height=7, units="in", dpi=300, type="cairo-png")

# Joyplots of posterior distribution of additional rings
##   No words could explain, no actions determine
##   Just watching the trees and the leaves as they fall
save_plot <- TRUE
dev.new()

add_rings_height %>%
  mutate(tree = fct_reorder(.f = tree, .x = add_rings_height, .fun = median)) %>% 
  ggplot(aes(x = add_rings_height, y = tree, height = stat(density))) + 
  geom_density_ridges(color = "white", fill = "black") +
  scale_x_log10(limits = c(1,NA), expand = expansion(0,0)) + 
  scale_y_discrete(expand = expansion(mult = c(0.02,0.08))) +
  xlab("Additional rings at root-shoot boundary") + ylab("Tree") + theme_bw(base_size = 16) + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(), 
        panel.background = element_rect(fill = "black", color = "black"), panel.border = element_blank(),
        strip.background = element_rect(color = "black", fill = "black"),
        strip.text = element_text(color = "white")) +
  facet_wrap(vars(patch), scales = "free")

if(save_plot) 
  ggsave(filename=here("analysis", "results", "coring-height-correction_joyplots.png"),
         width=7, height=7, units="in", dpi=300, type="cairo-png")





#===========================================================================
# ADDITIONAL FIGURES
#===========================================================================

# Additional rings vs. height, all data
plot(harv$height, harv$add_rings)

# Scatterplots add_rings ~ Height on Stem, grouped by patch, symbolized by sapling
harv %>% ggplot(aes(height, add_rings)) + geom_point(aes(color=tree)) +
  xlab("Height on Stem (cm)") + ylab("No. Add. Rings to Base") + 
  facet_wrap(vars(patch), ncol = 1) + theme_bw()

# Scatterplots add_rings ~ Height on Stem, grouped by plot, symbolized by sapling
harv %>% filter(plot=="BF03") %>%
  ggplot(aes(height, add_rings)) + geom_point(aes(color=tree)) +
  xlab("Height on Stem (cm)") + ylab("No. Add. Rings to Base")

# log(add_rings)
harv_log <- harv %>% filter(add_rings>0) %>%
  mutate(AddRings_log = log(add_rings))

head(harv_log)

harv_log %>% filter(plot=="F03") %>%
  ggplot(aes(height, AddRings_log)) + geom_point(aes(color=tree)) +
  xlab("Height on Stem (cm)") + ylab("No. Add. Rings to Base")



#===========================================================================
# ADDITIONAL CODE
#===========================================================================

### Poisson GLMs

# link = log
poi1 <- stan_glm(add_rings ~ height, family = poisson(link="log"), data = harv)
print(poi1, digit=3)

plot(harv$height, harv$add_rings)
curve(exp(coef(poi1)[1] + coef(poi1)[2]*x), add=TRUE)

yrep_1 <- posterior_predict(poi1)
n_sims <- nrow(yrep_1)
subset <-sample(n_sims, 100)
ppc_dens_overlay(harv$add_rings, yrep_1)

# link = identity
poi2 <- stan_glm(add_rings ~ height, family = poisson(link="identity"), data = harv)
print(poi2, digit=3)

plot(harv$height, harv$add_rings)
curve(coef(poi2)[1] + coef(poi2)[2]*x, add=TRUE)

yrep_2 <- posterior_predict(poi2)
n_sims <- nrow(yrep_2)
subset <-sample(n_sims, 100)
ppc_dens_overlay(harv$add_rings, yrep_2)


# link = identity, patch as fixed effect
poi3 <- stan_glm(add_rings ~ height + height:patch, family = poisson(link="identity"), data = harv)
print(poi3, digit=3)

plot(harv$height, harv$add_rings)
curve(coef(poi2)[1] + coef(poi2)[2]*x, add=TRUE)

yrep_2 <- posterior_predict(poi2)
n_sims <- nrow(yrep_2)
subset <-sample(n_sims, 100)
ppc_dens_overlay(harv$add_rings, yrep_2)

