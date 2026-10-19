## ----------------------------------------------------------
##
## Script name: N-Mixture Model Data Prep
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

## Package Loading
library(dplyr)
library(ggplot2)
library(data.table)
library(sf)
library(tidyverse)
library(here)
## -----------------------------------------------------------

## Load data
pfw <- read_csv(here("Data/PFW_AllYr_TestSpp.csv"))

## Check out the data quickly
glimpse(pfw)

## Let's start with a manageable area to begin
## Wisconsin has about 105k submitted lists
pfw <- pfw |> 
  filter(SUBNATIONAL1_CODE %in% c("US-CA", "US-OR")) |> 
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326, remove = F)

## Let's visualize
us <- st_read("/Users/srk252/Documents/GIS/Vector_Data/states.shp")
us <- us |> filter(STATE_NAME %in% c("California", "Oregon")) |> st_transform(crs = 4326)

ggplot() + 
  geom_sf(data = pfw) + 
  geom_sf(data = us, color = "blue", fill = NA)

## Let's set a sampling time period where we want to assume "closure"
## Prime winter months DEC - FEB?
pfw <- pfw |>
  filter(Month %in% c(12, 1, 2)) |> 
  arrange(LOC_ID, Month, PROJ_PERIOD_ID) |> 
  group_by(LOC_ID, PROJ_PERIOD_ID) |> 
  mutate(NSurv = n()) |> 
  ungroup()
  
length(unique(pfw$LOC_ID))

## We should get a sense for the average gap in surveys
## Per the PFW protocol we should have a gap of atleast
## 5 days between surveys

## Below we can bin the data into specific SSPs

## Set the PFW years
syears <- 1989:2024

## Start date based on filter above
start_dates <- lubridate::ymd(paste0(1988:2023, "-", 12, "-", 01))

## Find the leap years
lyrs <- leap_year(1989:2024)

## Get the end dates respecting the leap year
end_dates <- vector(mode = "list", length = length(lyrs))
for(i in 1:length(lyrs)){
  end_dates[[i]] <- lubridate::ymd(ifelse(lyrs[i], paste0(syears[i], "-", 02, "-", 29), paste0(syears[i], "-", 02, "-", 28)))
}
end_dates <- do.call(c, end_dates)

## Find the max number of dates between the hypothetical start and end date
maxdyr <- end_dates - start_dates
ssps <- as.numeric(maxdyr)/7

## Make a function to split the dates into equal bins
tspl_f <- function(x, breaks = 7){
  tspl <- split(seq(start_dates[x], end_dates[x]), ceiling(seq_along(seq(start_dates[x], end_dates[x]))/breaks))
  tspl <- lapply(tspl, range)
  tspl <- data.frame(SSP = names(tspl), sdate = do.call(c, lapply(tspl, min)), edate = do.call(c, lapply(tspl, max)))
  return(tspl)
}

## Split the dates up
tspl <- lapply(seq_along(syears), tspl_f)
tspl <- do.call(rbind, tspl)

## Create effort and timing information and associate potential SSPs for binning if desired
pfw <- pfw |> 
  mutate(Month = ifelse(str_detect(Month, "^\\d{1}$"), as.character(paste0("0", as.character(Month))), as.character(Month))) |> 
  mutate(Day = ifelse(str_detect(Day, "^\\d{1}$"), as.character(paste0("0", as.character(Day))), as.character(Day))) |> 
  mutate(Date = lubridate::ymd(paste0(Year, "-", Month, "-", Day))) |> 
  group_by(LOC_ID, PROJ_PERIOD_ID) |> 
  arrange(LOC_ID, PROJ_PERIOD_ID, Date) |> 
  mutate(SInterval = max(Date) - min(Date)) |> 
  mutate(LagInt = Date - lag(Date)) |> 
  mutate(
    SSP = map_int(
      Date, ~ as.integer(tspl$SSP[.x %within% interval(tspl$sdate, tspl$edate)])
    )
  )

## Take a look
hist(as.numeric(pfw$SInterval))
min(as.numeric(pfw$LagInt), na.rm = T)
hist(pfw$NSurv)
hist(pfw$EFFORT_HRS_ATLEAST)
hist(log(pfw$lesgol))
hist(log(pfw$pinsis))
hist(log(pfw$evegro))

## Good split across bins
table(pfw$SSP)

## Let's see how the distribution of count sites within a bin are
pfw |> 
  group_by(LOC_ID, PROJ_PERIOD_ID, SSP) |> 
  mutate(BinCheck = n()) |> 
  ggplot() + 
  geom_histogram(aes(x = BinCheck))

## Check the distribution on the number of people doing at least two checks (for detection prob)
## Some folk submitted 40+ checklists?!
pfw |> 
  ggplot() + 
  geom_histogram(aes(x = NSurv)) + 
  geom_vline(aes(xintercept = 2), color = "red") +
  geom_vline(aes(xintercept = mean(pfw$NSurv)), color = "blue") +
  theme_bw()

ggplot(data = pfw) + 
  geom_sf(aes(color = NSurv)) + 
  geom_sf(data = us, color = "blue", fill = NA) + 
  facet_wrap(~PROJ_PERIOD_ID)

ggplot(data = pfw |> mutate(lesgol_pa = ifelse(lesgol > 0, 1, 0))) + 
  geom_sf(aes(color = lesgol_pa)) +
  scale_color_viridis_c() +
  geom_sf(data = us, color = "blue", fill = NA) + 
  facet_wrap(~PROJ_PERIOD_ID)


## Let's take format the data from SpAbundance to start with only one covariate for year
pfw <- pfw |> 
  ungroup() |> 
  group_by(PROJ_PERIOD_ID) |> 
  mutate(NSiteYr = n()) |> 
  mutate(ProjYr = as.numeric(str_extract(PROJ_PERIOD_ID, "\\d+")))

## Participant count
ggplot(data = pfw |> 
         group_by(PROJ_PERIOD_ID) |> 
         summarise(NSiteYr = max(NSiteYr),
                   across(amegfi:purfin, ~sum(.x, na.rm = T)),
                   ProjYr = max(ProjYr)) |> 
         pivot_longer(cols = amegfi:purfin,
                      names_to = "Species",
                      values_to = "Max_Count"), 
       aes(x = ProjYr, y = Max_Count/NSiteYr)) +
  geom_point() + 
  geom_line() + 
  theme_bw() + 
  ylab("Total Obs./Site") + 
  xlab("Year") +
  facet_wrap(~Species, scales = "free")


## Export the data for WI
pfw <- pfw |> 
  relocate(amegfi:purfin, .after = colnames(pfw)[length(colnames(pfw))]) |> 
  st_drop_geometry()

write_csv(pfw, file = "Data/PFW_CA_OR_Test.csv")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Format the data for Nimble
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pfw <- read_csv(here("Data/PFW_CA_OR_Test.csv"))

## For the first pass we will do a "stacked" model where years
## are long wise but we will use site as a raneff to deal with 
## repeated measures
glimpse(pfw)
max(pfw$SSP)
hist(pfw$NSurv)
quantile(pfw$NSurv)

## For a first pass we can drop the sites that went over the protocol
pfwf <- pfw |> 
  #filter(is.na(LagInt) | LagInt <= 15) |> 
  group_by(LOC_ID, PROJ_PERIOD_ID, SSP) |> 
  slice_sample(n = 1) |> 
  group_by(LOC_ID, PROJ_PERIOD_ID) |> 
  mutate(LagInt = as.numeric(Date - lag(Date))) |> 
  arrange(LOC_ID, Date) |> 
  ungroup()

hist(pfwf$LagInt)
abline(v = 7)
## Shape into the M x J matrix
pfww <- pfwf |> 
  mutate(SiteYr = paste0(LOC_ID, "_", ProjYr)) |> 
  select(SiteYr, SSP, lesgol) |>
  tidyr::pivot_wider(names_from = SSP,
                     values_from = lesgol,
                     names_glue = "Rep_{SSP}") |> 
  relocate(SiteYr, paste0("Rep_", 1:13)) |> 
  mutate(Complete = rowSums(!is.na(across(starts_with("Rep_"))))) |> 
  filter(Complete >= 3) |> 
  select(-Complete)

## Site Year combos
gsites <- unique(pfww$SiteYr)

## Snow cover
snow_cov <- pfwf |> 
  ungroup() |>
  mutate(SiteYr = paste0(LOC_ID, "_", ProjYr)) |>
  filter(SiteYr %in% gsites) |> 
  select(SiteYr, SSP, SNOW_DEP_ATLEAST) |>
  tidyr::pivot_wider(names_from = SSP,
                     values_from = SNOW_DEP_ATLEAST,
                     names_glue = "Snow_{SSP}",
                     values_fill = NA) |> 
  relocate(SiteYr, paste0("Snow_", 1:13))|>
  tibble::column_to_rownames("SiteYr")

## Climate (2-day survey average)
clim_cov <- fread(here("Data/PFW_All_SurveyData_Daymet.csv"))
clim_cov <- as.data.frame(clim_cov)
## Select the vars needed
clim_cov <- clim_cov |> select(SUB_ID,
                               LOC_ID,
                               LATITUDE,
                               LONGITUDE,
                               Date,
                               tmin,
                               tmax,
                               prcp) |> 
  mutate(Date = as.Date(Date))

# Merge
clim_cov <- pfwf |> 
  left_join(clim_cov)

clim_cov <- clim_cov |> 
  ungroup() |>
  mutate(SiteYr = paste0(LOC_ID, "_", ProjYr)) |>
  filter(SiteYr %in% gsites) |> 
  select(SiteYr, SSP, tmin, tmax, prcp) |>
  tidyr::pivot_wider(names_from = SSP,
                     values_from = c(tmin, tmax, prcp),
                     names_glue = "{.value}_{SSP}",
                     values_fill = NA,
                     values_fn = mean)

## Take the vars of interest
tmin <- clim_cov |> 
  select(SiteYr, starts_with("tmin")) |> 
  tibble::column_to_rownames("SiteYr")

tmax <- clim_cov |> 
  select(SiteYr, starts_with("tmax")) |> 
  tibble::column_to_rownames("SiteYr")

prcp <- clim_cov |> 
  select(SiteYr, starts_with("prcp")) |> 
  tibble::column_to_rownames("SiteYr")


## Effort
eff_cov <- pfwf |> 
  ungroup() |>
  mutate(SiteYr = paste0(LOC_ID, "_", ProjYr)) |>
  filter(SiteYr %in% gsites) |> 
  select(SiteYr, SSP, EFFORT_HRS_ATLEAST) |>
  tidyr::pivot_wider(names_from = SSP,
                     values_from = EFFORT_HRS_ATLEAST,
                     names_glue = "EFF_{SSP}",
                     values_fill = NA) |> 
  relocate(SiteYr, paste0("EFF_", 1:13)) |> 
  tibble::column_to_rownames("SiteYr")

## Year data
year_cov <- pfwf |>
  mutate(SiteYr = paste0(LOC_ID, "_", ProjYr)) |>
  filter(SiteYr %in% gsites) |> 
  select(SiteYr, LOC_ID, ProjYr) |> 
  #mutate(Year = ProjYr - min(ProjYr)) |> 
  distinct(SiteYr, LOC_ID, ProjYr)


## Get aligned matrices for filtering data
## by incompleteness across the three sets
## y, eff, and snow

## y-matrix
y_wide <- as.matrix(pfww |> select(starts_with("Rep_")))
siteyr_ids <- pfww$SiteYr
rownames(y_wide) <- siteyr_ids

## Ensure effort and snow have rows in the exact same SiteYr order
## and columns in the exact same replicate order as y_wide.
eff_names  <- sub("^Rep_", "EFF_", colnames(y_wide))
snow_names <- sub("^Rep_", "Snow_", colnames(y_wide))
tmin_names <- sub("^Rep_", "tmin_", colnames(y_wide))
tmax_names <- sub("^Rep_", "tmax_", colnames(y_wide))
prcp_names <- sub("^Rep_", "prcp_", colnames(y_wide))

stopifnot(all(siteyr_ids %in% rownames(eff_cov)))
stopifnot(all(siteyr_ids %in% rownames(snow_cov)))
stopifnot(all(eff_names %in% colnames(eff_cov)))
stopifnot(all(snow_names %in% colnames(snow_cov)))
stopifnot(all(tmin_names %in% colnames(tmin)))
stopifnot(all(tmax_names %in% colnames(tmax)))
stopifnot(all(prcp_names %in% colnames(prcp)))

eff_wide <- as.matrix(
  eff_cov[match(siteyr_ids, rownames(eff_cov)), eff_names, drop = FALSE]
)

snow_wide <- as.matrix(
  snow_cov[match(siteyr_ids, rownames(snow_cov)), snow_names, drop = FALSE]
)

tmin_wide <- as.matrix(
  tmin[match(siteyr_ids, rownames(tmin)), tmin_names, drop = FALSE]
)

tmax_wide <- as.matrix(
  tmax[match(siteyr_ids, rownames(tmax)), tmax_names, drop = FALSE]
)

prcp_wide <- as.matrix(
  prcp[match(siteyr_ids, rownames(prcp)), prcp_names, drop = FALSE]
)

# Give all matrices identical dimnames for easier checking
rownames(eff_wide) <- siteyr_ids
rownames(snow_wide) <- siteyr_ids
colnames(eff_wide) <- colnames(y_wide)
colnames(snow_wide) <- colnames(y_wide)

stopifnot(identical(dim(y_wide), dim(eff_wide)))
stopifnot(identical(dim(y_wide), dim(snow_wide)))


## Indexing based on complete obs for three vars
obs_available <- !is.na(y_wide) &
  !is.na(eff_wide) &
  !is.na(snow_wide) &
  !is.na(tmin_wide) &
  !is.na(tmax_wide) &
  !is.na(prcp_wide)

valid_idx <- which(obs_available, arr.ind = TRUE)

## House keeping variables for Nimble
nObs <- nrow(valid_idx)
SiteYr <- unique(rownames(valid_idx))
nSiteYr <- length(SiteYr)
nSite <- length(unique(gsub("_\\d{4}", "", SiteYr)))
site_id_by_yr <- gsub("_\\d{4}", "", SiteYr)
site_id_by_obs <- gsub("_\\d{4}", "", names(valid_idx[,1]))
siteyr_by_obs <- names(valid_idx[,1])
Site <- unique(gsub("_\\d{4}", "", SiteYr))

## Indices for the nested vars
## Site by year
site_ind <- match(site_id_by_yr, Site)
## Site by rep/obs
site_obs_ind <- match(site_id_by_obs, Site)
## Site yr ind
siteyr_ind <- match(siteyr_by_obs, SiteYr)

obs_site <- valid_idx[, 1]
obs_rep  <- valid_idx[, 2]

## Extract aligned vectors directly using the same matrix indices
y_long    <- y_wide[valid_idx]
eff_long  <- eff_wide[valid_idx]
snow_long <- snow_wide[valid_idx]
tmin_long  <- tmin_wide[valid_idx]
tmax_long <- tmax_wide[valid_idx]
prcp_long <- prcp_wide[valid_idx]

## Checks
stopifnot(length(y_long) == nObs)
stopifnot(length(eff_long) == nObs)
stopifnot(length(snow_long) == nObs)
stopifnot(length(tmin_long) == nObs)
stopifnot(length(tmax_long) == nObs)
stopifnot(length(prcp_long) == nObs)

stopifnot(!anyNA(y_long))
stopifnot(!anyNA(eff_long))
stopifnot(!anyNA(snow_long))
stopifnot(!anyNA(tmin_long))
stopifnot(!anyNA(tmax_long))
stopifnot(!anyNA(prcp_long))

## Trim the y_wide to match
y_wide <- y_wide[rownames(y_wide) %in% SiteYr,]

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Abundance covariates
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
year_cov <- year_cov |> 
  mutate(SiteYr = paste0(LOC_ID, "_", ProjYr)) |> 
  filter(SiteYr %in% unique(names(obs_site))) |> 
  arrange(match(SiteYr, unique(names(obs_site)))) |>
  mutate(Year = ProjYr - min(ProjYr)) |> 
  pull(Year)

# ## ***********************************************************
# ##
# ## Section Notes: We have a series of sites where 
# ## we are missing both count and cov data bc of non-surveyed
# ## data.
# ## Option 1: Ragged array
# ## Option 2: Imputation from prior
# ##
# ## ***********************************************************  
# 
# ## 40% of sites x reps are missing...
# ## 44% of sites don't have effort covariates
# ## 48% of sites don't have snow cover data
# sum(is.na(pfww[,2:ncol(pfww)])) / prod(dim(pfww[,2:ncol(pfww)]))
# sum(is.na(eff_cov)) / prod(dim(eff_cov))
# sum(is.na(snow_cov)) / prod(dim(snow_cov))
# 
# ## Move into a ragged array
# pfwa <- as.data.frame(pfww)
# rownames(pfwa) <- pfwa$SiteYr
# pfwa <- pfwa[,2:ncol(pfwa)]
# valid_idx <- which(!is.na(pfwa), arr.ind = TRUE)
# nObs <- nrow(valid_idx)
# nSite <- length(unique(gsub("_\\d{4}", "", rownames(valid_idx)))) 
# nSiteYr <- length(unique(rownames(valid_idx)))                 
# validSamp <- unique(rownames(valid_idx))
# 
# ## Indexing
# obs_site <- as.integer(valid_idx[,1])
# obs_rep <- as.integer(valid_idx[,2])
# 
# site_ids <- sort(unique(obs_site))
# site_map <- setNames(seq_along(site_ids), site_ids)
# obs_site <- unname(site_map[as.character(obs_site)])
# max(obs_site)
# nSiteYr == length(unique(obs_site))
# 
# y_wide <- pfwa |> select(starts_with("Rep"))
# 
# ## Make a long count data
# y_long <- matrix(NA, nObs)
# eff_long <- matrix(NA, nObs)
# snow_long <- matrix(NA, nObs)
# for(v in 1:nObs){
#   i <- obs_site[v]
#   j <- obs_rep[v]
#   
#   ## Count
#   y_long[v] <- y_wide[i,j]
#   
#   ## Effort
#   eff_long[v] <- eff_cov[i,j]
#   
#   ## Snow
#   snow_long[v] <- snow_cov[i,j]
# }
# 
# sum(is.na(y_long))
# sum(is.na(eff_long))
# sum(is.na(snow_long))

## Nimble Ragged Data
nim.rag <- list(
  y = y_long,
  y_wide = y_wide,
  eff = as.numeric(factor(eff_long)),
  snow = as.numeric(factor(snow_long)),
  tmin = tmin_long,
  tmin_sc = as.vector(scale(tmin_long)),
  tmax = tmax_long,
  tmax_sc = as.vector(scale(tmax_long)),
  prcp = prcp_long,
  prcp_sc = as.vector(scale(prcp_long)),
  nEff = length(unique(eff_long)),
  nSnow = length(unique(snow_long)),
  nSiteYr = nSiteYr,
  SiteYr = SiteYr,
  nSite = nSite,
  Site = Site,
  nObs = nObs,
  site_obs = obs_site,
  site_ind = site_ind,
  site_obs_ind = site_obs_ind,
  siteyr_ind = siteyr_ind,
  rep_obs = obs_rep,
  year = year_cov
)

str(nim.rag)

## Write the df
saveRDS(nim.rag, file = here("Data/NMix_Input/PFW_CA_OR_LESGOL_Ragged.RDS"))
