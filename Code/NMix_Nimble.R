## ----------------------------------------------------------
##
## Script name: PFW N-Mixture Model
##
## Script purpose:
##
## Author: Spencer R Keyser
##
## Date Created: 2026-10-08
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

## Load renv for cli runs
renv::load()

## Package Loading
library(dplyr)
library(ggplot2)
library(data.table)
library(nimble)
## -----------------------------------------------------------

## -----------------------------------------------------------
##
## Begin Section: Nimble models 
##
## -----------------------------------------------------------
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: N-Mixture in stacked design
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Simplest N-Mixture model w/ covariates stacked design (site x yr)
fw_nmix <- nimbleCode({
  
  ####### Priors ######
  ## Abundance
  beta0 ~ dnorm(0, 0.1)
  beta1 ~ dnorm(0, 0.1)
  
  ## Site-level raneff
  for(s in 1:nSite){
    eps.loc[s] ~ dnorm(0, tau.eps.loc)
  }
  sd.eps.loc ~ dunif(0, 2)
  tau.eps.loc <- pow(sd.eps.loc, -2)
  
  ## Detection
  alpha0 ~ dnorm(0, 0.1)
  
  ## Contrasts for the factor effects on detection
  alpha_eff[1] <- 0
  for(e in 2:nEff){
    alpha_eff[e] ~ dnorm(0, 0.1)
  }
  
  alpha_snow[1] <- 0
  for(g in 2:nSnow){
    alpha_snow[g] ~ dnorm(0, 0.1)
  }
  
  ## Diffuse normal priors on tmin and prcp (continuous)
  alpha_tmin ~ dnorm(0, 0.1)
  alpha_prcp ~ dnorm(0, 0.1)
  
  ###### Likelihood ######
  ## Abundance
  for(i in 1:nSiteYr){
    N[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta0 + beta1 * year[i] + eps.loc[site_ind[i]]
  }
  for(v in 1:nObs){
    ## Observation
    y[v] ~ dbin(p_obs[v], N[siteyr_ind[v]])
    
    logit(p_obs[v]) <- alpha0 +
      alpha_eff[eff[v]] +
      alpha_snow[snow[v]] +
      alpha_tmin * tmin_sc[v] + 
      alpha_prcp * prcp_sc[v]
  }
}
)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: NMix with full site x year grid
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Simplest N-Mixture model w/ covariates stacked design (site x yr)
fw_nmix_t <- nimbleCode({
  
  ####### Priors ######
  ## Abundance
  beta0 ~ dnorm(0, 0.1)
  beta1 ~ dnorm(0, 0.1)
  
  ## Site-level raneff
  for(s in 1:nSite){
    eps.loc[s] ~ dnorm(0, tau.eps.loc)
  }
  sd.eps.loc ~ dunif(0, 2)
  tau.eps.loc <- pow(sd.eps.loc, -2)
  
  ## Detection
  alpha0 ~ dnorm(0, 0.1)
  
  alpha_eff[1] <- 0
  for(e in 2:nEff){
    alpha_eff[e] ~ dnorm(0, 0.1)
  }
  
  alpha_snow[1] <- 0
  for(g in 2:nSnow){
    alpha_snow[g] ~ dnorm(0, 0.1)
  }
  
  alpha_tmin ~ dnorm(0, 0.1)
  alpha_prcp ~ dnorm(0, 0.1)
  
  ###### Likelihood ######
  ## Abundance
  for(i in 1:nSite){
    for(t in 1:nYear){
      N[i,t] ~ dpois(lambda[i,t])
      log(lambda[i,t]) <- beta0 + beta1 * year[t] + eps.loc[site_ind[i]]
    }
  }
  
  ## Observation
  for(v in 1:nObs){
    
    y[v] ~ dbin(p_obs[v], N[siteyr_ind[v]])
    
    logit(p_obs[v]) <- alpha0 +
      alpha_eff[eff[v]] +
      alpha_snow[snow[v]] 
      alpha_tmin * tmin[v] + 
      alpha_prcp * prcp[v] 
  }
  
  ## Derived parameters
  for(t in 1:nYear){
    bigN[t] <- sum(N[1:nSite,t])
  }
  
}
)

## -----------------------------------------------------------
##
## Begin Section: Data load and model settings
##
## -----------------------------------------------------------

## Load the data
#dat <- readRDS(here("Data/NMix_Input/PFW_WI_DEJU_Ragged.RDS"))
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript PFW_NMix.R <input_data.rds> <output_samples.rds>",
    call. = FALSE
  )
}

input_file <- args[1]
output_file <- args[2]

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file, call. = FALSE)
}

dat <- readRDS(input_file)

if (!is.list(dat)) {
  stop("Input RDS must contain a list of NIMBLE data/constants.", call. = FALSE)
}

message("Loaded data from: ", input_file)
str(dat)

## Inits
Nst <- apply(dat$y_wide, 1, function(x) max(x, na.rm = T))
inits <- function() {
  list(
    N = Nst,
    alpha0 = rnorm(1, 0, 1),
    alpha_tmin = rnorm(1, 0, 1),
    alpha_prcp = rnorm(1, 0, 1),
    beta0 = rnorm(1, 0, 1),
    beta1 = rnorm(1, 0, 1),
    
    eps.loc = rnorm(dat$nSite, 0, 0.1),
    sd.eps.loc = runif(1, 0.05, 1)
    
    #eps.obs = rnorm(dat$nSiteYr, 0, 0.1),
    #sd.eps.obs = runif(1, 0.05, 1)
  )
}

## MCMC Settings
nburnin = 10000
niter = 40000
nthin = 20
nchain = 3

## Specify constants
constants <- list(
  nSite = dat$nSite,
  nSiteYr = dat$nSiteYr,
  nObs = dat$nObs,
  nEff = dat$nEff,
  nSnow = dat$nSnow,
  year = dat$year,
  site_ind = dat$site_ind,
  siteyr_ind = dat$siteyr_ind,
  eff = dat$eff,
  snow = dat$snow,
  tmin_sc = dat$tmin_sc,
  prcp_sc = dat$prcp_sc,
  NMax = 150
)

## Build Model
FWcmod <- nimbleModel(
  code = fw_nmix,
  data = list(y = dat$y),
  constants = constants,
  inits = inits()
)

## Monitors
pmon <- c(
  "lambda",
  "N",
  "p_obs",
  "alpha0",
  "alpha_eff",
  "alpha_snow",
  "alpha_tmin",
  "alpha_prcp",
  "beta0",
  "beta1",
  "eps.loc",
  #"eps.obs",
  "sd.eps.loc"
  #"sd.eps.obs"
)

## Compilation
Cmodel <- compileNimble(FWcmod)

## Configure
conf <- configureMCMC(FWcmod, monitors = pmon, enableWAIC = T)

## Build it
Rmcmc <- buildMCMC(conf)

## Compile
Cmcmc <- compileNimble(Rmcmc, project = Cmodel)

## Run the samples
samples <- runMCMC(Cmcmc,
                   niter = niter,
                   nburnin = nburnin,
                   thin = nthin,
                   nchains = nchain,
                   samplesAsCodaMCMC = FALSE,
                   summary = FALSE,
                   WAIC = TRUE)

## Write the samples
saveRDS(samples, output_file)
message("Saved samples to: ", output_file)
