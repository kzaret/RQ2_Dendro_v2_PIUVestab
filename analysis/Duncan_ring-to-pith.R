## UNCERTAINTY ANALYSIS OF DUNCAN RING-TO-PITH CORRECTIONS
##
## Estimate uncertainty in Duncan (1989) adjustments for cores that missed the pith  
## Posterior of rings-to-pith, given sample of 3 innermost ring widths, will be used
## as a component of informative prior on tree age


#===========================================================================
# SETUP
#===========================================================================

library(dplyr)
library(forcats)
library(rstanarm)
library(bayesplot)
library(shinystan)
library(ggplot2)
library(ggridges)
library(matrixStats)
library(here)
bayesplot_theme_set(bayesplot::theme_default())
options(mc.cores = parallel::detectCores(logical = FALSE) - 1)
if(.Platform$OS.type == "windows") options(device = windows)

#===========================================================================
# DATA
#===========================================================================

# 3 innermost ring widths for pithless-cored trees in all plots
inner_rings <- read.csv(here("data","PithlessTCs_InnerRingWidths.csv"), header = TRUE)
inner_rings <- inner_rings %>% rename(patch = Patch, tree = Core_ID, ring = Ring_ID, width = Ring.Width) %>% 
  mutate(ring = substring(ring, nchar(ring)))

# arc length (L) and height (h) of first full ring
duncan <- read.csv(here("data","Duncan_estimates.csv"), header = TRUE, na.strings = c("NA","#DIV/0!"))
duncan <- duncan %>% rename(tree = `Core..`, h = H, width_1st_full_ring = X1st.full.ring,
                            r = Missing.Distance.To.Pith..DTP.,
                            width_inner3_rings = Width.innermost.3.rings..WI3R.,
                            mean_width_inner3_rings = Mean.WI3R,
                            rings_to_pith = Estimated...rings.to.pith..RTP.,
                            subtract_if_not_arching = If.last.vis.ring.not.arching..extra.rings.to.subtract,
                            rings_to_pith_adj = Adj.rings.to.pith,
                            duplicate = Adj.rings.to.pith.values,
                            downward_curvature = downward.curvature) %>% select(-duplicate) %>% 
  mutate(patch = inner_rings$patch[match(tree, inner_rings$tree)], .before = tree) %>% 
  filter(!is.na(patch)) # temporarily (?) throw out trees not present in inner_rings


#===========================================================================
# HIERARCHICAL LINEAR MODELS FOR RING WIDTH
#===========================================================================

# Global intercept + fixed effect of patch
# Tree-level random intercept
# Lognormal likelihood (log-transform data)
duncan_lmer <- stan_lmer(log(width) ~ patch + (1 | tree), data = inner_rings, 
                         chains = getOption("mc.cores"), iter = 5000, warmup = 1000)

print(duncan_lmer,2)
summary(duncan_lmer)

## Diagnostic plots ##

# Simulate draws from posterior predictive distribution
yrep <- posterior_predict(duncan_lmer)
indx <- sample(nrow(yrep), 2000)

# Posterior predictive check: density overlay
ppc_dens_overlay_grouped(duncan_lmer$y, yrep[indx[1:500],], group = duncan_lmer$glmod$fr$patch)

# Posterior predictive check: histograms
ppc_hist(duncan_lmer$y, yrep[indx[1:3],])

# Posterior predictive check: mean
ppc_stat_grouped(duncan_lmer$y, yrep[indx,], group = duncan_lmer$glmod$fr$patch, stat = mean)

# Posterior predictive check: SD
ppc_stat_grouped(duncan_lmer$y, yrep[indx,], group = duncan_lmer$glmod$fr$patch, stat = sd)

# Normal QQ plot of tree-level random intercept point estimates, grouped by patch
ranef(duncan_lmer)$tree %>% rename(intercept = `(Intercept)`) %>% 
  mutate(patch = factor(tapply(duncan$patch, duncan$tree, unique))) %>% 
  ggplot(aes(sample = intercept)) + stat_qq(size = 2) + geom_qq_line() +
  theme_bw() + facet_wrap(vars(patch), ncol = 2)

# Normal QQ plot of observation-level residuals point estimates, grouped by patch
data.frame(duncan_lmer$glmod$fr, resid = resid(duncan_lmer)) %>% 
  ggplot(aes(sample = resid)) + stat_qq(size = 2) + geom_qq_line() +
  theme_bw() + facet_wrap(vars(patch), ncol = 2)


#===========================================================================
# DUNCAN RING-TO-PITH CORRECTIONS
#===========================================================================

# Compute Duncan ring-to-pith correction for each tree as r/d_hat,
# where the missing radius (r) is already calculated as r = L^2/(8*h) + h/2
# and d_hat is the median inner ring width (more robust than mean for lognormally
# distributed data) estimated by the hierarchical model
# Repeat for each draw from the posterior

# median of lognormal is exp(mu)
# if wanted mean instead, would need exp(mu + 0.5*sigma^2)
d_hat <- exp(posterior_linpred(duncan_lmer, newdata = duncan))
rtp <- sweep(1/d_hat, 2, duncan$r, "*")  # keep continuous version for plotting
rings_to_pith <- round(rtp)


#===========================================================================
# FIGURES
#===========================================================================

# Dotplots of ring width by core
# Overlay violin plots of PPD
save_plot <- TRUE
if(save_plot) {
  png(filename=here("analysis", "results", "duncan_lmer_fits.png"),
      width=7, height=7, units="in", res=300, type="cairo-png")
} else dev.new()

dat <- data.frame(iter = rep(indx, ncol(yrep)),
                  patch = rep(duncan_lmer$glmod$fr$patch, each = length(indx)),
                  tree = rep(duncan_lmer$glmod$fr$tree, each = length(indx)),
                  width = as.vector(exp(yrep[indx,]))) %>% 
  mutate(tree = fct_reorder(.f = tree, .x = width, .fun = mean))
  
inner_rings %>% ggplot(aes(x = tree, y = width)) +
  geom_violin(aes(x = tree, y = width, group = tree),
              data = dat, color = "darkgray", fill = "darkgray", alpha = 0.8) +
  geom_jitter(shape = 16, alpha = 0.5, size = 1, width = 0.1, height = 0) +
  scale_y_log10() + xlab("Tree") + ylab("Inner ring width (cm)") +
  theme_bw(base_size = 16) + theme(axis.text.x = element_blank()) + 
  facet_wrap(vars(patch), scales = "free_x")

if(save_plot) dev.off()

# Violin plots of posterior distribution of ring-to-pith estimates
save_plot <- TRUE
if(save_plot) {
  png(filename=here("analysis", "results", "duncan_rings-to-pith_violins.png"),
      width=7, height=7, units="in", res=300, type="cairo-png")
} else dev.new()

data.frame(iter = rep(indx, ncol(rtp)),
           patch = rep(duncan$patch, each = length(indx)),
           tree = rep(duncan$tree, each = length(indx)),
           rtp = as.vector(rtp[indx,])) %>%
  mutate(tree = fct_reorder(.f = tree, .x = rtp, .fun = median)) %>% 
  ggplot(aes(x = tree, y = rtp, group = tree)) +
  geom_violin(color = "darkgray", fill = "darkgray", alpha = 0.8) +
  xlab("Tree") + ylab("Estimated rings to pith") +
  theme_bw(base_size = 16) + theme(axis.text.x = element_blank()) + 
  facet_wrap(vars(patch), scales = "free_x")

if(save_plot) dev.off()

# Joyplots of posterior distribution of ring-to-pith estimates
##   No words could explain, no actions determine
##   Just watching the trees and the leaves as they fall
save_plot <- TRUE
if(save_plot) {
  png(filename=here("analysis", "results", "duncan_rings-to-pith_joyplots.png"),
      width=7, height=7, units="in", res=300, type="cairo-png")
} else dev.new()

data.frame(iter = rep(indx, ncol(rtp)),
           patch = rep(duncan$patch, each = length(indx)),
           tree = rep(duncan$tree, each = length(indx)),
           rtp = as.vector(rtp[indx,])) %>%
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

if(save_plot) dev.off()


         