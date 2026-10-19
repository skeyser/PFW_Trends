## ----------------------------------------------------------
##
## Script name: PFW Environmental Extraction
##
## Script purpose:
##
## Author: Spencer R Keyser
##
## Date Created: 2026-10-12
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
library(terra)
library(readr)
library(here)
## -----------------------------------------------------------

## Load in the PFW Survey info
pfws <- read_csv(here("Data/PFW_All_SurveyData.csv"))
glimpse(pfws)

## Process climate data for the 2-day survey period
## Start with the average climate conditions of those days
## Month and Day will define the start so we can add a day
## for the 2-day period

## This might need to be done CLI
## Make a simple text file for lat long
dmt_coords <- pfws |> 
  select(LOC_ID, LATITUDE, LONGITUDE) |> 
  distinct(LOC_ID, .keep_all = TRUE) |> 
  mutate(LOC_ID = paste0(LOC_ID, ".csv"))

## Write as a text file
write_delim(dmt_coords, file = "/Users/srk252/Documents/GIS/Env_Data/Daymet/dmt_coords.txt", delim = ",", col_names = FALSE)

## Need to re-run process so let's save the trouble by pulling out the completed runs
runs <- list.files("/Users/srk252/Documents/GIS/Env_Data/Daymet/", pattern = "*.csv")

dmt_coords_sub <- dmt_coords |> 
  filter(!LOC_ID %in% runs)

## Check we got them all
nrow(dmt_coords) - nrow(dmt_coords_sub) == length(runs)

## save the new data
write_delim(dmt_coords_sub, file = "/Users/srk252/Documents/GIS/Env_Data/Daymet/dmt_coords_sub.txt", delim = ",", col_names = FALSE)


## Missed in the download
bad_files <- readRDS(here("Data/Daymet_Bad_Extract.rds"))

dmt_coords_bad <- dmt_coords |> 
  filter(LOC_ID %in% bad_files)

## Quick map to see if these are coastal sites
library(sf)
us <- st_read("/Users/srk252/Documents/GIS/Vector_Data/Westernhemi.shp")

dmt_coords_bad |>
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326) |>
  ggplot() +
  geom_sf(data = us) +
  geom_sf()


## Check we got them all
nrow(dmt_coords) - nrow(dmt_coords_sub) == length(runs)

## save the new data
write_delim(dmt_coords_bad, file = "/Users/srk252/Documents/GIS/Env_Data/Daymet/dmt_coords_sub.txt", delim = ",", col_names = FALSE)
