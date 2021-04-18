## UNCERTAINTY ANALYSIS OF DUNCAN RINGS-TO-PITH CORRECTIONS
##
## Estimate uncertainty in Duncan (1989) adjustments for cores that missed the pith  
## Posterior of rings-to-pith, given sample of 3 innermost ring widths, will be used
## as a component of informative prior on tree age


#===========================================================================
# SETUP
#===========================================================================

## @knitr setup
library(matrixStats)
library(dplyr)
library(forcats)
library(rstanarm)
library(bayesplot)
library(shinystan)
library(ggplot2)
library(ggridges)
library(here)
options(mc.cores = parallel::detectCores(logical = FALSE) - 1)
if(.Platform$OS.type == "windows") options(device = windows)
## @knitr

#===========================================================================
# DATA
#===========================================================================

## @knitr data
# 3 innermost ring widths for pithless-cored trees in all plots
inner_rings_raw <- read.csv(here("data","PithlessTCs_InnerRingWidths.csv"), header = TRUE)
inner_rings <- inner_rings_raw %>% 
  rename(patch = Patch, tree = Core_ID, ring = Ring_ID, width = Ring.Width) %>% 
  mutate(ring = substring(ring, nchar(ring)))

# arc length (L) and height (h) of first full ring
# remove trees not present in inner_rings, those with unknown plot ("RLARGE01", "RLARGE02"),
# and those in plot "R0X" ("R0XT01a")
duncan_raw <- read.csv(here("data","Duncan_estimates.csv"), header = TRUE, na.strings = c("NA","#DIV/0!"))
duncan <- duncan_raw %>% rename(tree = `Core..`, h = H, width_1st_full_ring = X1st.full.ring,
                                r = Missing.Distance.To.Pith..DTP.,
                                width_inner3_rings = Width.innermost.3.rings..WI3R.,
                                mean_width_inner3_rings = Mean.WI3R,
                                rings_to_pith = Estimated...rings.to.pith..RTP.,
                                subtract_if_not_arching = If.last.vis.ring.not.arching..extra.rings.to.subtract,
                                rings_to_pith_adj = Adj.rings.to.pith,
                                duplicate = Adj.rings.to.pith.values,
                                downward_curvature = downward.curvature) %>% 
  mutate(patch = inner_rings$patch[match(tree, inner_rings$tree)], .before = tree) %>% 
  select(-duplicate) %>% filter(!is.na(patch) & !(tree %in% c("R0XT01a", "RLARGE01", "RLARGE02"))) 
## @knitr


#===========================================================================
# HIERARCHICAL LINEAR MODELS FOR RING WIDTH
#===========================================================================

# Global intercept + fixed effect of patch
# Tree-level random intercept
# Lognormal likelihood (log-transform data)
## @knitr fit_lmer
duncan_lmer <- stan_lmer(log(width) ~ patch + (1 | tree), data = inner_rings, 
                         chains = getOption("mc.cores"), iter = 5000, warmup = 1000)

## @knitr print_fit_lmer
print(duncan_lmer,2)
## @knitr
summary(duncan_lmer)

## DIAGNOSTIC PLOTS ##

# Simulate draws from posterior predictive distribution
## @knitr sim_ppd
yrep <- posterior_predict(duncan_lmer)
indx <- sample(nrow(yrep), 2000)
## @knitr

# Posterior predictive check: density overlay (note: takes a while)
## @knitr ppd_dens_overlay_grouped
ppc_dens_overlay_grouped(duncan_lmer$y, yrep[indx[1:500],], group = duncan_lmer$glmod$fr$patch)
## @knitr

# Posterior predictive check: histograms
## @knitr ppc_hist
ppc_hist(duncan_lmer$y, yrep[indx[1:3],])
## @knitr

# Posterior predictive check: mean
## @knitr ppc_mean_grouped
ppc_stat_grouped(duncan_lmer$y, yrep[indx,], group = duncan_lmer$glmod$fr$patch, stat = mean)
## @knitr

# Posterior predictive check: SD
## @knitr ppc_sd_grouped
ppc_stat_grouped(duncan_lmer$y, yrep[indx,], group = duncan_lmer$glmod$fr$patch, stat = sd)
## @knitr

# Normal QQ plot of tree-level random intercept point estimates, grouped by patch
## @knitr qqnorm_intercept_patch
ranef(duncan_lmer)$tree %>% rename(intercept = `(Intercept)`) %>% 
  mutate(patch = factor(tapply(duncan$patch, duncan$tree, unique))) %>% 
  ggplot(aes(sample = intercept)) + stat_qq(size = 2) + geom_qq_line() +
  theme_bw() + facet_wrap(vars(patch), ncol = 2)
## @knitr

# Normal QQ plot of observation-level residuals point estimates, grouped by patch
## @knitr qqnorm_resid_patch
data.frame(duncan_lmer$glmod$fr, resid = resid(duncan_lmer)) %>% 
  ggplot(aes(sample = resid)) + stat_qq(size = 2) + geom_qq_line() +
  theme_bw() + facet_wrap(vars(patch), ncol = 2)
## @knitr


#===========================================================================
# DUNCAN RINGS-TO-PITH CORRECTIONS
#
# Compute Duncan rings-to-pith correction for each tree as r/d_hat,
# where the missing radius (r) is already calculated as r = L^2/(8*h) + h/2
# and d_hat is the median inner ring width (more robust than mean for lognormally
# distributed data) estimated by the hierarchical model
# Repeat for each draw from the posterior
#===========================================================================

# median of lognormal is exp(mu)
# if wanted mean instead, would need exp(mu + 0.5*sigma^2)
## @knitr calc_rtp
d_hat <- exp(posterior_linpred(duncan_lmer, newdata = duncan))
rtp <- sweep(1/d_hat, 2, duncan$r, "*")  # keep continuous version for plotting

# Attach posterior draws to data
rings_to_pith <- data.frame(duncan[rep(1:nrow(duncan), each = nrow(rtp)), c("patch","tree")],
                            iter = rep(1:nrow(rtp), nrow(duncan)), 
                            rings_to_pith = as.vector(rtp))
rownames(rings_to_pith) <- 1:nrow(rings_to_pith)
## @knitr

# Save posterior draws and stanfit objects
save(list = c("duncan_lmer", "rings_to_pith"), 
     file = here("analysis","results","duncan_rings-to-pith.RData"))


#===========================================================================
# FIGURES
#===========================================================================

# Dotplots of ring width by core
# Overlay violin plots of PPD
save_plot <- TRUE
dev.new()

## @knitr duncan_lmer_ppd_violins
dat <- data.frame(iter = rep(indx, ncol(yrep)),
                  patch = rep(duncan_lmer$glmod$fr$patch, each = length(indx)),
                  tree = rep(duncan_lmer$glmod$fr$tree, each = length(indx)),
                  width = as.vector(exp(yrep[indx,]))) %>% 
  mutate(tree = fct_reorder(.f = tree, .x = width, .fun = mean))
  
inner_rings %>% ggplot(aes(x = tree, y = width)) +
  geom_violin(aes(x = tree, y = width),
              data = dat, color = "darkgray", fill = "darkgray", alpha = 0.8) +
  geom_jitter(shape = 16, alpha = 0.5, size = 1, width = 0.1, height = 0) +
  scale_y_log10() + xlab("Tree") + ylab("Inner ring width (cm)") + theme_bw(base_size = 16) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
  facet_wrap(vars(patch), scales = "free_x")
## @knitr

if(save_plot) 
  ggsave(filename=here("analysis", "results", "duncan_lmer_ppd_violins.png"),
         width=7, height=7, units="in", dpi=300, type="cairo-png")

# Credible intervals of posterior distribution of rings-to-pith estimates
save_plot <- TRUE
dev.new()

## @knitr duncan_rings-to-pith_intervals
mcmc_intervals_data(rtp, prob = 0.8, prob_outer = 0.95) %>% 
  mutate(patch = duncan$patch, tree = fct_reorder(duncan$tree, .x = m, .fun = identity)) %>% 
  ggplot(aes(x = tree, y = m)) +
  geom_linerange(aes(ymin = l, ymax = h), size = 1.5, color = "darkgray") +
  geom_linerange(aes(ymin = ll, ymax = hh), size = 0.5, color = "darkgray") +
  geom_point(pch = 16, size = 1) +
  xlab("Tree") + ylab("Estimated rings to pith") + 
  theme_bw(base_size = 16) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        panel.grid = element_blank()) + 
  facet_wrap(vars(patch), scales = "free")
## @knitr

if(save_plot) 
  ggsave(filename=here("analysis", "results", "duncan_rings-to-pith_intervals.png"),
         width=7, height=7, units="in", dpi=300, type="cairo-png")

# Joyplots of posterior distribution of rings-to-pith estimates
##   No words could explain, no actions determine
##   Just watching the trees and the leaves as they fall
save_plot <- TRUE
dev.new()

## @knitr duncan_rings-to-pith_joyplots
data.frame(iter = rep(nrow(rtp), ncol(rtp)),
           patch = rep(duncan$patch, each = nrow(rtp)),
           tree = rep(duncan$tree, each = nrow(rtp)),
           rtp = as.vector(rtp)) %>%
  mutate(tree = fct_reorder(.f = tree, .x = rtp, .fun = median)) %>% 
  ggplot(aes(x = rtp, y = tree, height = stat(density))) + 
  geom_density_ridges(color = "white", fill = "black") +
  scale_x_log10(expand = expansion(0,0)) + scale_y_discrete(expand = expansion(mult = c(0.02,0.08))) +
  xlab("Estimated rings to pith") + ylab("Tree") + theme_bw(base_size = 16) + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(), 
        panel.background = element_rect(fill = "black", color = "black"), panel.border = element_blank(),
        strip.background = element_rect(color = "black", fill = "black"),
        strip.text = element_text(color = "white")) +
  facet_wrap(vars(patch), scales = "free")
## @knitr

if(save_plot) 
  ggsave(filename=here("analysis", "results", "duncan_rings-to-pith_joyplots.png"),
         width=7, height=7, units="in", dpi=300, type="cairo-png")


         