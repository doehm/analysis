# libraries
library(tidyverse)
library(patchwork)
library(invgamma)

set.seed(1600)

# set par names
pars <- c("b0", "b1", "sd")

# model set up
nobs <- 50
b0 <- 2
b1 <- 4
sd <- 2

# make df
pars_true <- tibble(
  name = pars,
  pars_true = c(b0, b1, sd)
)

# simulate data
x <- rnorm(nobs, 10, sqrt(5))
y <- b0 + b1*x + rnorm(nobs, 0, sd)

# plot data
tibble(x, y) |>
  ggplot(aes(x, y)) +
  geom_point() +
  geom_abline(slope = b1, intercept = b0, col = "red")


# posterior function ------------------------------------------------------

posterior <- function(par){

  # take the input parameters of the intercept, slope and std dev
  b0 <- par[1]
  b1 <- par[2]
  sd <- par[3]

  # compute the expected value given the input parameters
  y_hat <- b0 + b1*x

  # compute the log likelihoods
  loglikelihoods <- dnorm(y, mean = y_hat, sd = sd, log = TRUE)

  # sum the log likelihoods and return
  sumll <- sum(loglikelihoods)

  # priors - non-formative
  b0_prior <- dnorm(b0, sd = 5, log = TRUE)
  b1_prior <- dnorm(b1, sd = 5, log = TRUE)
  sd2_prior <- dinvgamma(sd^2, shape = 1, log = TRUE) # inverse gamma prior on the variance

  # now return the sum of all components
  sum(sumll, b0_prior, b1_prior, sd2_prior)
}


# Metropolis-Hastings algorithm -------------------------------------------

metropHastingsMCMC <- function(theta0, proposal_sd, iter, burnin){

  # initialise the chain
  chain = matrix(NA, nrow=iter+1, ncol=length(theta0))
  chain[1,] = theta0
  acceptance <- array(0, c(iter, length(theta0)))

  # each iteraction take a draw from the proposal for each parameter
  # calculate teh acceptance probability
  # either accept or reject the proposal given the probability r
  for (i in 1:iter){
    theta_star <- rnorm(length(pars), mean = chain[i,], sd = proposal_sd)
    r <- exp(posterior(theta_star) - dnorm(theta_star, mean = chain[i,], sd = proposal_sd, log = TRUE) -
               posterior(chain[i,]) + dnorm(chain[i,], mean = theta_star, sd = proposal_sd, log = TRUE))

    # now that r is a vector of 3 values representing each parameter we need to accept/reject each individually
    for(k in 1:length(pars)){
      if (runif(1) < r[k]){
        chain[i+1,k] = theta_star[k]
        acceptance[i,k] <- 1
      }else{
        chain[i+1,k] = chain[i,k]
      }
    }
  }
  cat("Acceptance rate:", mean(acceptance[burnin:iter]))
  chain[burnin:iter,]
}


# run sim -----------------------------------------------------------------

theta0 <- c(1, 2, 1)
prop_sd <- c(0.5, 0.1, 0.25)
its <- 100000
burn <- 25000
sims <- metropHastingsMCMC(theta0, prop_sd, its, burn)


# plot --------------------------------------------------------------------

df_sims <- sims |>
  as_tibble() |>
  set_names(pars) |>
  mutate(iter = 1:n()) |>
  pivot_longer(c(b0, b1, sd))

df_mean <- df_sims |>
  group_by(name) |>
  summarise(value = mean(value))

g_hist <- df_sims |>
  ggplot(aes(value)) +
  geom_histogram(fill = "darkcyan", bins = 60) +
  geom_vline(aes(xintercept = pars_true), pars_true, col = "red", lty = 2) +
  geom_vline(aes(xintercept = value), df_mean, lty = 2) +
  facet_wrap(~name, nrow = length(pars), scales = "free_x")

g_chains <- df_sims |>
  ggplot(aes(iter, value)) +
  geom_line() +
  geom_hline(aes(yintercept = pars_true), pars_true, col = "red", lty = 2) +
  geom_hline(aes(yintercept = value), df_mean, lty = 2) +
  facet_wrap(~name, nrow = length(pars), scales = "free_y")

g_hist + g_chains

