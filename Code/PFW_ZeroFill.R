## ----------------------------------------------------------
##
## Script name: Project FeederWatch Data Prep
##
## Script purpose:
##
## Author: Spencer R Keyser
##
## Date Created: 2026-10-07
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
library(tidyverse)
## -----------------------------------------------------------

## ***********************************************************
##
## Section Notes: Data fields
## LOC_ID = LOCATION ID OF CHECKLIST
## LAT/LONG = COORDS (GCS)
## SUBNATIONAL1_CODE = GEOPOLITICAL ASSOC.
## SUB_ID = ID FOR THE 2-DAY SAMPLING EVENT
## OBS_ID = UNIQUE FOR EACH ROW OF DATA SUBMISSION
## Month/Day/Year = Time of survey
## PORJ_PERIOD_ID = PFW SEASON ID
## SPECIES_CODE = eBird Species ID
## alt_full_spp_code = 
## HOW_MANY = COUNTS
## VALID = FLAG FOR VALID (1) OR INVALID (0)
## REVIEWED = 
## DAY1_AM/DAY1_PM/DAY2_AM/DAY2_PM
## EFFORT_HRS_ATLEAST = COMBINED EFFORT PER SURVEY
## SNOW_DEP_ATLEAST = SNOW AT FW SITE
##
## 
## ***********************************************************

source(here("Code/PFW_Funs.R"))

## Load PFW files
pfw.files <- list.files("/Users/srk252/Documents/PFW_Raw/", pattern = "Public", full.names = T)

#pfw <- fread("/Users/srk252/Documents/PFW_Raw/PFW_all_2006_2010_May2024_Public.csv")
pfw <- lapply(pfw.files, fread)
pfw <- lapply(pfw, as.data.frame)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: PFW Survey Information
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
glimpse(pfw)
vars <- c("LOC_ID", "LATITUDE", "LONGITUDE", "SUBNATIONAL1_CODE",
          "SUB_ID", "Month", "Day", "Year", "PROJ_PERIOD_ID")

pfw_spat_get <- function(x, vars){
  o <- x |> 
    select(all_of(vars)) |> 
    distinct(SUB_ID, .keep_all = TRUE)
}

pfw_surv <- lapply(pfw, function(x) pfw_spat_get(x=x, vars=vars))
pfw_surv <- do.call(rbind, pfw_surv)
#write_csv(pfw_surv, file = here("Data/PFW_All_SurveyData.csv"))
rm(pfw_surv)
## -----------------------------------------------------------
##
## End Section: PFW Survey Information
##
## -----------------------------------------------------------


## All species helper
uniq_sp <- function(x){
  unique(x$SPECIES_CODE)
}

usp <- lapply(pfw, uniq_sp)
usp <- do.call(c, usp)
usp <- unique(usp)

## Species translation table
spt <- read.csv("/Users/srk252/Documents/PFW_Raw/PFW_spp_translation_table_May2024.csv")
spt |> filter(american_english_name %in% c("Evening Grosbeak", "Purple Finch", "Carolina Wren", "Lesser Goldfinch"))

## Find the number of unique species
unique(pfw$SPECIES_CODE)
unique(pfw$alt_full_spp_code)

## Apply the zero fill process
## Discrepancy between mine and Wes's function
## on how to collate species based on the Review Flag
## To truly get a perfect square matrix we need to provide
## all detected species across each PFW sampling era


pfw <- lapply(pfw, function(x) ZeroFillPFW(x,
                                           SpeciesCodes = c("evegro",
                                                            "purfin",
                                                            "houfin",
                                                            "daejun",
                                                            "amegfi",
                                                            "blujay",
                                                            "dowwoo",
                                                            "pinsis",
                                                            "lesgol",
                                                            "carwre"),
                                           rollup = T))

pfw <- lapply(pfw, function(x) x |> filter(VALID == 1))

## Bind all the years together for a large df
pfw <- do.call(rbind, pfw)
glimpse(pfw)

## Save file for model testing
readr::write_csv(pfw, file = here("Data/PFW_AllYr_TestSpp.csv"))
