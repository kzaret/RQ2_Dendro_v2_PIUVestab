#---------------------------------------------------------------------------------
# SIMULATION TESTING FOR POISSON STATE-SPACE MODELS
#
# Simulate data from the basic univariate Poisson SS model and
# fit the pseudo-data with the same model, the equivalent 
# Poisson-multinomial model (Poisson for total count, multinomial
# for conditional counts), and a multinomial model (conditioned on
# total count).
#
# Simulate data from a novel model that adds observation error to
# the sample of times generated from the Poisson SS model. Fit
# the observed counts with the generating model or the standard
# multinomial model, and compare fit if the true times were known.
#---------------------------------------------------------------------------------

options(device = windows)

## @knitr setup
library(rstan)
library(shinystan)
library(yarrr)
library(matrixStats)
library(here)
options(mc.cores = parallel::detectCores() - 1)
## @knitr

if(file.exists(here("analysis","results","Poisson_SS.RData")))
  load(here("analysis","results","Poisson_SS.RData"))

#---------------------------------------------------------------------------------
# POISSON STATE-SPACE MODEL
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# Simulate univariate Poisson state-space model
# x[t] = x[t-1] + w[t], w[t] ~ N(0,sigma), x[0] ~ N(0,sigma)
# y[t] ~ Pois(exp(x[t]))
#---------------------------------------------------------------------------------

## @knitr Poisson_SS_sim
set.seed(34567)
sigma <- 0.3
TT <- 200

x <- vector("numeric",TT)
x[1] <- rnorm(1,0,sigma)
for(tt in 2:TT)
  x[tt] <- rnorm(1, x[tt-1], sigma)
y <- rpois(TT, exp(x))
## @knitr

#---------------------------------------------------------------------------------
# Fit standard Poisson state-space model
#---------------------------------------------------------------------------------

## @knitr fit_pois
fit_pois <- stan(file = here("analysis","Poisson_SS.stan"), 
                 data = list(T = TT, y = y), pars = c("sigma","x"),
                 chains = 3, iter = 2000, warmup = 1000,
                 control = list(adapt_delta = 0.99, max_treedepth = 12))
## @knitr print_fit_pois
print(fit_pois, pars = "sigma", probs = c(0.025,0.5,0.975))
## @knitr

#---------------------------------------------------------------------------------
# Fit Poisson-multinomial state-space model
# (do not condition on total)
#---------------------------------------------------------------------------------

## @knitr fit_pois_mn
fit_pois_mn <- stan(file = here("analysis","Poisson_multinomial_SS.stan"), 
                    data = list(T = TT, y = y), pars = c("sigma","x"),
                    init = function() list(sigma = runif(1,0.1,0.5)),
                    chains = 3, iter = 2000, warmup = 1000,
                    control = list(adapt_delta = 0.99, max_treedepth = 12))
## @knitr print_fit_pois_mn
print(fit_pois_mn, pars = "sigma", probs = c(0.025,0.5,0.975))
## @knitr

#---------------------------------------------------------------------------------
# Fit multinomial model to state-space Poisson counts
# (condition on total)
#---------------------------------------------------------------------------------

## @knitr fit_mn
fit_mn <- stan(file = here("analysis","multinomial_SS.stan"), 
               data = list(T = TT, y = y), pars = c("sigma","pi","lambda"),
               chains = 3, iter = 2000, warmup = 1000,
               control = list(adapt_delta = 0.99, max_treedepth = 12))
## @knitr print_fit_mn
print(fit_mn, pars = "sigma", probs = c(0.025,0.5,0.975))
## @knitr

#---------------------------------------------------------------------------------
# Plot data, states, and fits
#---------------------------------------------------------------------------------

dev.new(width = 7, height = 5)
## @knitr plot_pois_mn
par(mar = c(5.1,4.1,2,1))
# Poisson-multinomial state-space model
lambda <- exp(as.matrix(fit_pois_mn, "x"))
plot(1:TT, exp(x), type = "l", col = "dodgerblue", lwd = 3, 
     cex.axis = 1.2, cex.lab = 1.5, cex.main = 1.5, xlab = "time", ylab = "count", 
     ylim = range(colQuantiles(lambda, probs = 0.975), y), yaxt = "n")
rug(seq(0,TT,10)[seq(0,TT,10) %% 50 != 0], side = 1, ticksize = -0.03) 
axis(2, at = 0:par("usr")[4], las = 1, cex.axis = 1.2)
points(1:TT, y, type = "h", col = transparent("black", 0.3))
polygon(c(1:TT, TT:1), 
        c(colQuantiles(lambda, probs = 0.025), rev(colQuantiles(lambda, probs = 0.975))),
        col = transparent("dimgray", 0.5), border = NA)
lines(1:TT, colMedians(lambda), col = "dimgray", lwd = 3)
legend("topright", bty = "n", text.col = "white", cex = 1.2,
       legend = expression(lambda[italic(t)], italic(y)[italic(t)], widehat(lambda[italic(t)])),
       pch = c(NA,"I",NA), lwd = c(3,NA,15), 
       col = c(NA, transparent("black", 0.3), transparent("dimgray", 0.5)))
legend("topright", bty = "n", cex = 1.2, 
       legend = expression(lambda[italic(t)], italic(y)[italic(t)], widehat(lambda[italic(t)])), 
       lwd = c(3,NA,3), col = c(transparent("dodgerblue", 0.3), NA, "dimgray"))
## @knitr


#---------------------------------------------------------------------------------
# POISSON STATE-SPACE MODEL WITH OBSERVATION ERROR IN TIMES
#---------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# Simulate univariate Poisson state-space model with observation error in times
# x[tau] = x[tau-1] + w[tau], w[tau] ~ N(0,sigma), x[0] ~ N(0,sigma)
# y[t] ~ Pois(exp(x[t])) <=> tau[i] ~ Multinom(1,pi) for i = 1, ..., N
# t[i] ~ Multinom(1, gamma[i]), 
# where gamma[i] is the observation error distribution for time i.
# Example: 
# gamma[i,j] = dgeom(tau[j] - t[i], r) <=> t[i] ~ tau[i] + Geom(r)
#---------------------------------------------------------------------------------

## @knitr Poisson_tobs_SS_sim
set.seed(321)
TT <- 200
N <- 500   # total sample size
sigma <- 0.3
r <- 0.2   # probability parameter for geometric obs error in time

x <- vector("numeric",TT)
x[1] <- rnorm(1,0,sigma)
for(tt in 2:TT)
  x[tt] <- rnorm(1, x[tt-1], sigma)
pi <- exp(x)/sum(exp(x))
chi <- as.vector(rmultinom(1,N,pi))  # true counts
tau <- rep(1:TT, times = chi)        # true times
tt <- pmin(tau + rgeom(N,r), TT)     # observed times
tab <- table(tt)
y <- replace(rep(0,TT), as.numeric(names(tab)), tab) # observed counts
## @knitr

#---------------------------------------------------------------------------------
# Fit standard Poisson-multinomial model to true, unknown times  
#---------------------------------------------------------------------------------

## @knitr fit_tau
fit_tau <- stan(file = here("analysis","Poisson_multinomial_SS.stan"), 
              data = list(T = TT, y = chi), pars = c("sigma","x"),
              chains = 3, iter = 2000, warmup = 1000,
              control = list(adapt_delta = 0.99, max_treedepth = 12))
## @knitr print_fit_tau
print(fit_tau, pars = "sigma", probs = c(0.025,0.5,0.975))
## @knitr

#---------------------------------------------------------------------------------
# Fit standard Poisson-multinomial model to observed times  
#---------------------------------------------------------------------------------

## @knitr fit_t
fit_t <- stan(file = here("analysis","Poisson_multinomial_SS.stan"), 
              data = list(T = TT, y = y), pars = c("sigma","x"),
              chains = 3, iter = 2000, warmup = 1000,
              control = list(adapt_delta = 0.99, max_treedepth = 12))
## @knitr print_fit_t
print(fit_t, pars = "sigma", probs = c(0.025,0.5,0.975))
## @knitr

#---------------------------------------------------------------------------------
# Fit multinomial model with obs error in time 
#---------------------------------------------------------------------------------

## @knitr fit_tobs
fit_tobs <- stan(file = here("analysis","Poisson_multinomial_tobs_SS.stan"), 
                 data = list(T = TT, y = y, r = r), pars = c("sigma","x"),
                 chains = 3, iter = 2000, warmup = 1000,
                 control = list(adapt_delta = 0.99, max_treedepth = 12))
## @knitr print_fit_tobs
print(fit_tobs, pars = "sigma", probs = c(0.025,0.5,0.975))
## @knitr

#---------------------------------------------------------------------------------
# Plot data, states, and fits
#---------------------------------------------------------------------------------

dev.new(width = 10, height = 6)

## @knitr plot_tobs
par(mfrow = c(2,2), mar = c(3,2.5,2.5,1), oma = c(1.5,2,0,0))

lambda <- N*pi
lambda_t <- exp(as.matrix(fit_t, "x"))
lambda_tau <- exp(as.matrix(fit_tau, "x"))
lambda_tobs <- exp(as.matrix(fit_tobs, "x"))
ul <- max(colQuantiles(lambda_t, probs = 0.975), 
          colQuantiles(lambda_tau, probs = 0.975),
          colQuantiles(lambda_tobs, probs = 0.975))

# states and observations
plot(1:TT, lambda, type = "l", col = "dodgerblue", lwd = 3, 
     las = 1, cex.axis = 1.2, cex.lab = 1.5, cex.main = 1.5, ylim = c(0,ul), 
     xlab = "", ylab = "count", main = "States and observations", font.main = 1, xpd = NA)
rug(seq(0,TT,10)[seq(0,TT,10) %% 50 != 0], side = 1, ticksize = -0.02) 
points(1:TT, chi, pch = 1, col = "dodgerblue")
points(1:TT, y, type = "h", col = transparent("black", 0.3))
legend("topleft", bty = "n", text.col = "white", cex = 1.2, pt.cex = 1,
       legend = expression(lambda[italic(t)], chi[italic(t)], italic(y)[italic(t)]),
       pch = c(NA,NA,"I"), col = c(NA, NA, transparent("black", 0.3)))
legend("topleft", bty = "n", cex = 1.2, pt.cex = 1,
       legend = expression(lambda[italic(t)], chi[italic(t)], italic(y)[italic(t)]), 
       pch = c(NA,1,NA), lty = c(1,NA,NA), lwd = c(3,NA,NA), 
       col = c("dodgerblue", "dodgerblue", NA))

# fit to true, unknown times
plot(1:TT, lambda, type = "l", col = "dodgerblue", lwd = 3,
     las = 1, cex.axis = 1.2, cex.lab = 1.5, cex.main = 1.5, ylim = c(0,ul), 
     xlab = "", ylab = "", main = bquote("Poisson-multinomial fit to" ~ chi[italic(t)]), 
     font.main = 1)
rug(seq(0,TT,10)[seq(0,TT,10) %% 50 != 0], side = 1, ticksize = -0.02) 
polygon(c(1:TT, TT:1), 
        c(colQuantiles(lambda_tau, probs = 0.025), 
          rev(colQuantiles(lambda_tau, probs = 0.975))),
        col = transparent("dimgray", 0.5), border = NA)
lines(1:TT, colMedians(lambda_tau), col = "dimgray", lwd = 3)
points(1:TT, chi, pch = 1, col = "dodgerblue")
legend("topleft", bty = "n", legend = expression(widehat(italic(lambda)[italic(t)])), 
       text.col = "white", cex = 1.2, lwd = 15, col = transparent("dimgray", 0.5))
legend("topleft", bty = "n", legend = expression(widehat(italic(lambda)[italic(t)])), 
       cex = 1.2, lwd = 3, col = "dimgray")

# fit to observed times
plot(1:TT, lambda, type = "l", col = "dodgerblue", lwd = 3, 
     las = 1, cex.axis = 1.2, cex.lab = 1.5, cex.main = 1.5, 
     ylim = c(0,ul), xlab = "time", ylab = "count", 
     main = bquote("Poisson-multinomial fit to" ~ italic(y)[italic(t)]), 
     font.main = 1, xpd = NA)
rug(seq(0,TT,10)[seq(0,TT,10) %% 50 != 0], side = 1, ticksize = -0.02) 
points(1:TT, y, type = "h", col = transparent("black", 0.3))
polygon(c(1:TT, TT:1), 
        c(colQuantiles(lambda_t, probs = 0.025), 
          rev(colQuantiles(lambda_t, probs = 0.975))),
        col = transparent("dimgray", 0.5), border = NA)
lines(1:TT, colMedians(lambda_t), col = "dimgray", lwd = 3)

# fit to observed times, accounting for obs error
plot(1:TT, lambda, type = "l", col = "dodgerblue", lwd = 3, 
     las = 1, cex.axis = 1.2, cex.lab = 1.5, cex.main = 1.5, 
     font.main = 1, ylim = c(0,ul), xlab = "time", ylab = "", 
     main = bquote("Time-uncertain Poisson-multinomial fit to" ~ italic(y)[italic(t)]), xpd = NA)
rug(seq(0,TT,10)[seq(0,TT,10) %% 50 != 0], side = 1, ticksize = -0.02) 
points(1:TT, y, type = "h", col = transparent("black", 0.3))
polygon(c(1:TT, TT:1), 
        c(colQuantiles(lambda_tobs, probs = 0.025), 
          rev(colQuantiles(lambda_tobs, probs = 0.975))),
        col = transparent("dimgray", 0.5), border = NA)
lines(1:TT, colMedians(lambda_tobs), col = "dimgray", lwd = 3)
## @knitr

#---------------------------------------------------------------------------------
# Plot observation error distribution
#---------------------------------------------------------------------------------

dev.new(width = 7, height = 5)

## @knitr plot_geom_obs
tau_i <- 50                      # true time index
r <- 0.2                         # probability parameter for geometric obs error in time 
p_t_i <- dgeom(1:TT - tau_i, r)  # P(t_i | tau_i, r)

par(mar = c(5.1,5.1,2,1))
barplot(p_t_i, las = 1, cex.axis = 1.2, cex.lab = 1.5, 
        col = "darkgray", border = "white", space = 0,
        xlim = c(1,TT), xaxs = "i", xaxt = "n", ylim = range(p_t_i)*1.05,
        xlab = bquote(italic(t)[italic(i)]), ylab = bquote(gamma[italic(it)]))
axis(1, at = c(1, seq(50, TT, 50)), cex.axis = 1.2)
rug(seq(0,TT,10)[seq(0,TT,10) %% 50 != 0], side = 1, ticksize = -0.02) 
arrows(x0 = tau_i, y0 = -0.035, y1 = -0.025, col = "dodgerblue", length = 0.1, lwd = 2, xpd = NA)
text(tau_i, -0.042, labels = bquote(tau[italic(i)]), cex = 1.5, col = "dodgerblue", xpd = NA)
box()
## @knitr


#---------------------------------------------------------------------------------
# SAVE STANFIT OBJECTS
#---------------------------------------------------------------------------------

save(list = ls()[sapply(ls(), function(x) do.call(class, list(as.name(x)))) == "stanfit"], 
     file = here("analysis","results","Poisson_SS.RData"))


