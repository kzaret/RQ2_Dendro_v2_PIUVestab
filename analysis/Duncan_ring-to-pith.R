## UNCERTAINTY ANALYSIS OF DUNCAN RING-TO-PITH CORRECTIONS
##
## Estimate uncertainty in Duncan (1989) adjustments for cores that missed the pith  
## Posterior of rings-to-pith, given sample of 3 innermost ring widths, will be used
## as a component of informative prior on tree age


#===========================================================================
# SETUP
#===========================================================================

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

# 3 innermost ring widths for pithless-cored trees in all plots
duncan <- read.csv(here("data","PithlessTCs_InnerRingWidths.csv"), header = TRUE)
duncan <- duncan %>% rename(patch = Patch, tree = Core_ID, ring = Ring_ID, width = Ring.Width) %>% 
  mutate(ring = substring(ring, nchar(ring)))


#===========================================================================
# HIERARCHICAL LINEAR MODELS FOR RING WIDTH
#===========================================================================

# Global intercept + fixed effect of patch
# Tree-level random intercept
# Lognormal likelihood (log-transform data)
duncan_lmer <- stan_lmer(log(width) ~ patch + (1 | tree), data = duncan, 
                         chains = getOption("mc.cores"), iter = 5000, warmup = 1000)

print(duncan_lmer,2)
summary(duncan_lmer)

## Diagnostic plots ##

yrep <- posterior_predict(duncan_lmer)
indx <- sample(nrow(yrep), 2000)

# Posterior predictive check: histograms
ppc_hist(duncan_lmer$y, yrep[indx[1:3],])

# Posterior predictive check: mean
ppc_stat_grouped(duncan_lmer$y, yrep[indx,], group = duncan_lmer$glmod$fr$patch, stat = mean)

# Posterior predictive check: SD
ppc_stat_grouped(duncan_lmer$y, yrep[indx,], group = duncan_lmer$glmod$fr$patch, stat = sd)

# Normal QQ plot of tree-level random intercept point estimates, grouped by patch
grp_intercept <- as.data.frame(duncan_lmer, regex_pars = "b\\[") %>% as.matrix()

colMedians(grp_intercept) %>% data.frame() %>% setNames("intercept") %>% 
  mutate(patch = factor(tapply(duncan$patch, duncan$tree, unique))) %>% 
  ggplot(aes(sample = intercept)) + stat_qq(size = 2) + geom_qq_line() +
  theme_bw() + facet_wrap(vars(patch), ncol = 2)


#===========================================================================
# DUNCAN RING-TO-PITH CORRECTIONS
#===========================================================================





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
                  width = as.vector(exp(yrep[indx,])))

duncan %>% ggplot(aes(x = tree, y = width)) +
  geom_violin(aes(x = tree, y = width, group = tree),
              data = dat, color = "darkgray", fill = "darkgray", alpha = 0.8) +
  geom_jitter(shape = 16, alpha = 0.5, size = 1, width = 0.1, height = 0) +
  scale_y_log10() + xlab("Tree") + ylab("Inner ring width (cm)") +
  theme_bw(base_size = 16) + theme(axis.text.x = element_blank()) + 
  facet_wrap(vars(patch), scales = "free_x")

if(save_plot) dev.off()


         