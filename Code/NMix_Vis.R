## ----------------------------------------------------------
##
## Script name: MCMC Viz
##
## Script purpose:
##
## Author: Spencer R Keyser
##
## Date Created: 2026-10-14
##
## Email: srk252@cornell.edu
##
## Github: https://github.com/skeyser
##
## -----------------------------------------------------------
##
## Notes:
##
##
## -----------------------------------------------------------

## Defaults
options(scipen = 6, digits = 4)

## -----------------------------------------------------------

## Package Loading
library(dplyr)
library(ggplot2)
library(data.table)
library(here)
library(coda)
library(tidybayes)
library(MCMCvis)

## -----------------------------------------------------------

dat <- readRDS(here("Data/NMix_Input/PFW_CA_OR_LESGOL_Ragged.RDS"))
Sites <- dat$SiteYr
SiteInd <- dat$siteyr_ind
samp <- readRDS(here("Data/NMix_Output/PFW_CA_OR_LESGOL_Test.rds"))

## Convert the chains into a list
mcmc_list <- mcmc.list(
  mcmc(samp$chain1), 
  mcmc(samp$chain2), 
  mcmc(samp$chain3)
)

## Remove samples
rm(samp)

## Check convergence
## Check convergence
convergence_check <- function(x, 
                              params, 
                              probs = c(0.055, 0.5, 0.945),
                              ...){
  for(par in params){
    print(paste("Checking", par))
    mod_sum <- MCMCvis::MCMCsummary(x,
                                    params = par,
                                    round = 2,
                                    exact = F,
                                    probs = probs,
                                    ...)
    rh <- mod_sum$Rhat > 1.1
    check <- any(rh, na.rm = T)
    if(!check){
      print(paste("All values < Rhat 1.1 for", par))
    } else {
      print(paste("Some values >= Rhat 1.1 for", par))
      failRhat <- mod_sum[which(rh),]
      return(failRhat)
    }
  }
}

## Check parameter convergence
convergence_check(x = mcmc_list, params = "alpha_snow", ISB = F)
convergence_check(x = mcmc_list, params = "alpha_eff", ISB = F)
convergence_check(x = mcmc_list, params = "beta0")
convergence_check(x = mcmc_list, params = "beta1")
convergence_check(x = mcmc_list, params = "eps", ISB = F)
convergence_check(x = mcmc_list, params = "sd.eps.loc", ISB = F)

## Trace plots
MCMCvis::MCMCtrace(mcmc_list,
                  params = "alpha0",
                  n.eff = T,
                  Rhat = T,
                  pdf = F)

MCMCvis::MCMCtrace(mcmc_list,
                   params = "alpha_eff",
                   n.eff = T,
                   Rhat = T,
                   pdf = F)

MCMCvis::MCMCtrace(mcmc_list,
                   params = "alpha_snow",
                   n.eff = T,
                   Rhat = T,
                   pdf = F)

MCMCvis::MCMCtrace(mcmc_list,
                   params = "beta0",
                   n.eff = T,
                   Rhat = T,
                   pdf = F)

MCMCvis::MCMCtrace(mcmc_list,
                   params = "beta1",
                   n.eff = T,
                   Rhat = T,
                   pdf = F)

MCMCvis::MCMCtrace(mcmc_list,
                   params = "sd.eps.loc",
                   n.eff = T,
                   Rhat = T,
                   pdf = F)

## Pull samples into a df
n_chains <- 3
n_iter <- nrow(mcmc_list[[1]])
str(mcmc_list)
samples_df <- mcmc_list |> 
  as.matrix() |> 
  as_tibble() |> 
  select(contains(c("beta", "alpha", "eps", "N"))) |> 
  mutate(.chain = rep(1:n_chains, each = n_iter),
         .iteration = rep(1:n_iter, n_chains)) |> 
  pivot_longer(cols = -c(.chain, .iteration),
               names_to = "variable",
               values_to = "value") |> 
  group_by(variable) |> 
  summarise(MeanVal = mean(value),
            QL = quantile(value, probs = 0.055),
            QH = quantile(value, probs = 0.945))

samples_df |> 
  filter(str_detect(variable, "^N")) |> 
  mutate(SiteID = SiteYr) |>
  separate(col = SiteID, into = c("Site", "Year"), sep = "_") |>
  mutate(Year = factor(Year)) |> 
  group_by(Year) |> 
  summarise(C = median(MeanVal),
            QL = median(QL),
            QH = median(QH)) |>
  ungroup() |>  
  ggplot(aes(x = Year, y = C, group = 1)) + 
  geom_point() +
  geom_pointrange(aes(ymin = QL, ymax = QH)) +
  geom_line() +
  theme_bw()

MCMCplot(mcmc_list, params = "beta1")
MCMCplot(mcmc_list, params = "beta0")
MCMCplot(mcmc_list, params = "alpha_eff")
MCMCplot(mcmc_list, params = "alpha_snow")

plogis(samples_df$MeanVal[samples_df$variable == "alpha0"] + samples_df$MeanVal[samples_df$variable == "alpha_eff[4]"]) - plogis(samples_df$MeanVal[samples_df$variable == "alpha0"])


# Combine chains; rows are posterior draws
post <- do.call(rbind, mcmc_list)

# Values used for year in the fitted model, e.g. 0, 1, 2, 3, ...
year_vals <- sort(unique(dat$year))

# For every posterior draw and each year, calculate expected abundance
year_effect_draws <- rbindlist(
  lapply(year_vals, function(yr) {
    
    lambda <- exp(
      post[, "beta0"] +
        post[, "beta1"] * yr
    )
    
    rellam <- exp(
      post[,"beta1"] * yr
    )
    
    data.table(
      year = yr,
      lambda = lambda,
      rl = rellam
    )
  })
)

# Summarize posterior distribution at each year
year_effect_summary <- year_effect_draws[, .(
  median = median(lambda),
  mean = mean(lambda),
  lower = quantile(lambda, 0.055),
  upper = quantile(lambda, 0.945),
  rlmedian = median(rl),
  rlmean = mean(rl),
  rllower = quantile(rl, 0.055),
  rlupper = quantile(rl, 0.945)
), by = year]

first_year <- 1989

year_effect_summary[, calendar_year := first_year + year]

year_effect_summary

lesgo_tr <- ggplot(year_effect_summary, aes(x = calendar_year, y = median)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    fill = "darkgreen",
    alpha = 0.25
  ) +
  geom_line(
    color = "darkgreen",
    linewidth = 1
  ) +
  geom_point(
    color = "darkgreen",
    size = 2
  ) +
  labs(
    x = "Year",
    y = "E[Count/site]",
    title = "Marginal effect of year on expected Lesser Goldfinch abundance in CA & OR",
    subtitle = "line = posterior median; ribbon = 89% credible interval"
  ) +
  theme_bw()

lesgo_tr_rl <- ggplot(year_effect_summary, aes(x = calendar_year, y = rlmedian-1)) +
  geom_ribbon(
    aes(ymin = rllower-1, ymax = rlupper-1),
    fill = "darkgreen",
    alpha = 0.25
  ) +
  geom_line(
    color = "darkgreen",
    linewidth = 1
  ) +
  geom_point(
    color = "darkgreen",
    size = 2
  ) +
  labs(
    x = "Year",
    y = "Expected % Increase Relative to 1989",
    title = "Marginal effect of year on expected Lesser Goldfinch abundance in CA & OR",
    subtitle = "line = posterior median; ribbon = 89% credible interval"
  ) +
  scale_y_continuous(labels = scales::percent) + 
  theme_bw()


ggsave(lesgo_tr, filename = here("Figures/Prelim/LESGOL_CA_OR_Trend.png"),
       height = 8, width = 8, dpi = 600)

ggsave(lesgo_tr_rl, filename = here("Figures/Prelim/LESGOL_CA_OR_RelPctTrend.png"),
       height = 8, width = 8, dpi = 600)

