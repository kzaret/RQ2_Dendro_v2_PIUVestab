### Coring height corrections:  regressions of ring counts on height along stem for harvested cross section.  To be used to estimate rings not captured due to coring height of cores.

library(mlmRev)#for use in mlm tutorial (https://mc-stan.org/users/documentation/case-studies/tutorial_rstanarm.html)
library(lme4)
library(rstanarm)
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default())


### tutorial: 
data(Gcsemv, package = "mlmRev")
summary(Gcsemv)

# Make Male the reference category and rename variable
Gcsemv$female <- relevel(Gcsemv$gender, "M")


# Use only total score on coursework paper 
GCSE <- subset(x = Gcsemv, 
               select = c(school, student, female, course))

# Count unique schools and students
J <- length(unique(GCSE$school))
N <- nrow(GCSE)

# Check structure of data frame
str(GCSE)

#varying intercept and slope with single predictor
M3 <- lmer(formula = course ~ 1 + female + (1 + female | school), 
           data = GCSE, 
           REML = FALSE)
summary(M3) 

#using stan to do above
M3_stanlmer <- stan_lmer(formula = course ~ female + (1 + female | school), 
                         data = GCSE,
                         seed = 349)

prior_summary(object = M3_stanlmer)

M3_stanlmer


?stan_glmer


#############################################################################


# Data:  harvested PIUV from all patches; pot_outlier = 1 if had been removed based on linear regressions and/or abnormal ring counts in consecutive rounds.

#Re Bog forest outliers:  For linear regression, had removed BF0802 entirely [9 cross sections] due to nonlinear growth pattern that differs from the rest of the harvested trees.  It’s also one of the older BF trees. Removed cross-sections 06, 05 and 04 from BF0902 (bc ring counts of 92, 91, 90 and 90 in a row and notes re 04 being rotten and may want to leave out). Removed ff bc high on stem and no cross sections in between:  BF0303_110to112, BF0305_100to102. Removed BF0202 [5 cross sections]. Removed BF0301_110to112.

help(priors, package = 'rstanarm')  #The poisson family function defaults to using the log link

vignette("priors", package = 'rstanarm')

harv <- read.csv("C:/Users/Marsh Hawk Zaret/Documents/Big_D/Data_Analysis/RQ1v2_PIUVestab/Data/RC_Height_cross_sections.csv", header = TRUE)

colnames(harv)
head(harv)

harv_glmer <- stan_glmer(formula = AddRings ~ Height_RC + Height_RC:Patch + (-1 + Height_RC | Sapling), data = harv, family = poisson(link = "identity")) #model includes intercept for fixed effects, but no intercept for random effect of individual harvested tree.  Note:  could add dummy data such that each tree has a basal CS (height = 0) and add. rings = 0 for each basal CS.

prior_summary(harv_glmer)
