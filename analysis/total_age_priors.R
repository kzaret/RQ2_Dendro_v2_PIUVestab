## INFORMATIVE PRIORS ON TOTAL TREE AGE
##
## Combine model-based posteriors of coring-height correction and 
## Duncan rings-to-pith correction with prior on missing outer rings
## to arrive at informative priors on the age / recruitment year of
## each tree
## Use individual age priors to reconstruct prior on annual recruitment


#===========================================================================
# SETUP
#===========================================================================

library(matrixStats)
library(dplyr)
library(tidyr)
library(forcats)
library(rstanarm)
library(ggplot2)
library(bayesplot)
library(ggridges)
library(here)
if(.Platform$OS.type == "windows") options(device = windows)

#===========================================================================
# DATA
#===========================================================================

# Load saved posterior draws
if(file.exists(here("analysis","results","coring_height_correction.RData")))
  load(here("analysis","results","coring_height_correction.RData"))

if(file.exists(here("analysis","results","duncan_rings-to-pith.RData")))
  load(here("analysis","results","duncan_rings-to-pith.RData"))

# Merge coring-height and Duncan corrections and associated data 
# and compute total age and establishment year
# Note that iter is just an index; iterations across different models are 
# obviously uncoupled
# NOTE: currently one tree ("F07T01b") in Duncan data set not in cores data,
#       so will be ignored for now
tree_age <- left_join(add_rings_height, rings_to_pith, by = c("patch","tree","iter")) %>% 
  mutate(rings_to_pith = replace_na(rings_to_pith, 0),
         age = round(ring_count + add_rings_height + rings_to_pith),
         establishment_year = year - age - 1) # -1 b/c sampling is early in calendar year

# Prior on establishment year of each tree
# (for time-uncertain state-space model)
tree_year <- pivot_wider(tree_age, id_cols = establishment_year, names_from = tree, 
                         values_from = iter, values_fn = length, values_fill = 0) %>% 
  complete(establishment_year = min(establishment_year):max(establishment_year)) %>% 
  as.data.frame() %>% arrange(establishment_year) %>% replace(is.na(.), 0)

tree_year[,-1] <- sweep(tree_year[,-1], 2, colSums(tree_year[,-1]), "/")

# Prior on number of trees established per year
recruitment <- as.data.frame(xtabs(~ iter + establishment_year + patch, tree_age)) %>% 
  rename(year = establishment_year, N = Freq) %>% 
  mutate(year = as.numeric(as.character(year)), iter = as.numeric(as.character(iter)))


#===========================================================================
# FIGURES
#===========================================================================

# Credible intervals of posterior distribution of age estimates
save_plot <- TRUE
if(save_plot) {
  png(filename=here("analysis", "results", "tree_age_intervals.png"),
      width=7, height=7, units="in", res=300, type="cairo-png")
} else dev.new()

tree_age %>% pivot_wider(id_cols = c(tree,iter), names_from = tree, values_from = age) %>% 
  select(-iter) %>% mcmc_intervals_data(prob = 0.8, prob_outer = 0.95) %>% 
  mutate(patch = tree_age$patch[match(parameter, tree_age$tree)], 
         tree = fct_reorder(parameter, .x = m, .fun = identity)) %>% 
  ggplot(aes(x = tree, y = m)) +
  geom_linerange(aes(ymin = l, ymax = h), size = 1.5, color = "darkgray") +
  geom_linerange(aes(ymin = ll, ymax = hh), size = 0.5, color = "darkgray") +
  geom_point(pch = 16, size = 1, alpha = 0.5) +
  labs(x = "Tree", y = "Age") + theme_bw(base_size = 16) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        panel.grid = element_blank()) + 
  facet_wrap(vars(patch), scales = "free_x")

if(save_plot) dev.off()

# Joyplots of posterior distribution of rings-to-pith estimates
save_plot <- TRUE
if(save_plot) {
  png(filename=here("analysis", "results", "tree_age_joyplots.png"),
      width=7, height=7, units="in", res=300, type="cairo-png")
} else dev.new()

tree_age %>% mutate(tree = fct_reorder(.f = tree, .x = age, .fun = median)) %>% 
  ggplot(aes(x = age, y = tree, height = stat(density))) + 
  geom_density_ridges(color = "white", fill = "black") +
  scale_x_log10(expand = expansion(0,0)) + scale_y_discrete(expand = expansion(mult = c(0.02,0.08))) +
  labs(x = "Age", y ="Tree") + theme_bw(base_size = 16) + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(), 
        panel.background = element_rect(fill = "black", color = "black"), panel.border = element_blank(),
        strip.background = element_rect(color = "black", fill = "black"),
        strip.text = element_text(color = "white")) +
  facet_wrap(vars(patch), scales = "free")

if(save_plot) dev.off()

# Credible intervals of posterior distribution of annual recruitment
# Note: takes a while
save_plot <- FALSE
if(save_plot) {
  png(filename=here("analysis", "results", "recruitment_intervals.png"),
      width=8, height=6, units="in", res=300, type="cairo-png")
} else dev.new(width=8, height=6)

recruitment %>% 
  pivot_wider(id_cols = c(patch,year,iter), names_from = c(patch,year), values_from = N, values_fill = 0) %>%
  select(-iter) %>% # b/c pars arg to mcmc_intervals_data() does not work as advertised
  mcmc_intervals_data(point_est = "mean", prob = 0.8, prob_outer = 0.95) %>%
  separate(parameter, c("patch","year"), sep = "_") %>% mutate(year = as.numeric(year)) %>% 
  ggplot(aes(x = year, y = m)) +
  geom_ribbon(aes(ymin = ll, ymax = hh), fill = "black", alpha = 0.3) +
  geom_line(alpha = 0.5, lwd = 1) +
  scale_y_continuous(breaks = seq(0, max(hh))) +
  labs(x = "Year", y = "Recruitment") + theme_bw(base_size = 16) + 
  theme(panel.grid = element_blank(), plot.margin = unit(c(5.5, 15, 5.5, 5.5), "points"),
        panel.spacing = unit(15, "points")) +
  facet_wrap(vars(patch), scales = "free_y")

if(save_plot) dev.off()




