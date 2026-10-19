## ----------------------------------------------------------
##
## Script name: Land Cover Download Script
##
## Script purpose:
##
## Author: Spencer R Keyser
##
## Date Created: 2026-09-28
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
## -----------------------------------------------------------

## -----------------------------------------------------------
##
## Begin Section: Load Data
##
## -----------------------------------------------------------

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Subsection: Land Cover Data Extraction
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Project FeederWatch years
years <- 1986:2025

dwnld_nlcd_fun <- function(years, 
                           overwrite = F,
                           droot = "/Users/srk252/Documents/GIS/Env_Data/NLCD/"){
  options(timeout = 800)
  for(yr in years){
    print(paste0("Downloading NLCD for year:", yr))
    urlhead <- "https://www.mrlc.gov/downloads/sciweb1/shared/mrlc/data-bundles/"
    urltail <- paste0("Annual_NLCD_LndCov_", yr, "_CU_C1V2.zip")
    url <- paste0(urlhead, urltail)
    dest <- paste0(droot, gsub(urlhead, "", url))
    
    if(overwrite){
      download.file(url = url, destfile = dest)
    } else {
      if(any(list.files(droot) == gsub("*.txt|*.tif|*.zip", "", urltail))){
        message("File name detected in destination dir. Skipping download. Set overwrite = T to download file anyways.")
        next 
      } else {
        print("File/Dir not detected, downloading...")
        download.file(url = url, destfile = dest)
      }
    }
  }
  options(timeout = 60)
}

dwnld_cdl_fun <- function(years, 
                           overwrite = F,
                           droot = "/Users/srk252/Documents/GIS/Env_Data/CDL/"){
  options(timeout = 600)
  for(yr in years){
    print(paste0("Downloading CDL for year:", yr))
    urlhead <- "https://www.nass.usda.gov/Research_and_Science/Cropland/Release/datasets/"
    urltail <- paste0(yr, "_30m_cdls.zip")
    url <- paste0(urlhead, urltail)
    dest <- paste0(droot, gsub(urlhead, "", url))
    
    if(overwrite){
      tryCatch({
        download.file(url = url, destfile = dest)
        },
        error = function(e){
          message("An error occurred")
          message(paste("Error:", e$message))
        })
    } else {
      if(any(list.files(droot) == gsub("*.txt|*.tif|*.zip", "", urltail))){
        message("File name detected in destination dir. Skipping download. Set overwrite = T to download file anyways.")
        next 
      } else {
        print("File/Dir not detected, downloading...")
        tryCatch({
          download.file(url = url, destfile = dest)
        },
        error = function(e){
          message("An error occurred")
          message(paste("Error:", e$message))
        })
      }
    }
  }
  options(timeout = 60)
}

## Download NLCD files
dwnld_nlcd_fun(years)
#dwnld_cdl_fun(years)

