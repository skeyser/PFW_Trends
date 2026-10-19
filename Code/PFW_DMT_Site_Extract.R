# **********************************************************
# Script Info ------------------------------------------------
#
# Script name: DMT site-level extractions
#
# Script purpose: Summarise daymet data for PFW sites at the site
# level
#
# Author: Spencer R Keyser
#
# Date Created: 2026-10-19
#
# Email: srk252@cornell.edu
#
# Github: https://github.com/skeyser
#
#
#
# Notes:
#
#
#

# Defaults
options(scipen = 6, digits = 4)

renv::load()

# Package Loading -------------------------------------------
library(tidyverse)
library(ggplot2)
library(data.table)
library(here)

# ***********************************************************
#
# Data Loading ####
#
# ***********************************************************

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript PFW_DMT_Site_Extract.R <pfw.csv> <output dir>",
    call. = FALSE
  )
}

pfws_input <- args[1]
# Read in the survey data
#pfws_input <- "Data/PFW_All_SurveyData.csv"
dmt_path <- "/Users/srk252/Documents/GIS/Env_Data/Daymet/"
output_file <- args[2]

if (!file.exists(pfws_input)) {
  stop("Input file does not exist: ", pfws_input, call. = FALSE)
}

pfws <- as.data.frame(fread(pfws_input))
pfws <- pfws |> 
  mutate(Month = ifelse(str_detect(Month, "^\\d{1}$"), as.character(paste0("0", as.character(Month))), as.character(Month))) |> 
  mutate(Day = ifelse(str_detect(Day, "^\\d{1}$"), as.character(paste0("0", as.character(Day))), as.character(Day))) |> 
  mutate(Date = lubridate::ymd(paste0(Year, "-", Month, "-", Day)))


message("Loaded data from: ", pfws_input)

# Find the Daymet files
# dmt.files <- list.files("/Users/srk252/Documents/GIS/Env_Data/Daymet",
#                         pattern = "^L[0-9]+\\.csv", full.names = T)

dmt.files <- list.files(dmt_path,
                        pattern = "^L[0-9]+\\.csv", 
                        full.names = T)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## ---- Daymet Link Site Level
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dmt_link_site <- function(pfw,
                          dmt,
                          start_date = "12-01-1980",
                          end_date = "02-29-2025"){
  
  dmt_check <- function(dmt){
    # Check whether the file has anything beyond the 6 metadata lines
    file_lines <- readr::read_lines(dmt, progress = FALSE)
    
    if (length(file_lines) <= 6 || !any(nzchar(trimws(file_lines[-seq_len(6)])))) {
      message("Site ", gsub(".csv", "", basename(dmt)), " has an empty Daymet file...skipping...")
      return(NULL)
    }
    dmt_tmp <- tryCatch(
      readr::read_csv(
        dmt,
        skip = 6,
        show_col_types = FALSE,
        progress = FALSE
      ),
      error = function(e) {
        message(
          "Could not read Daymet file for site ", basename(dmt),
          ": ", conditionMessage(e),
          "...skipping..."
        )
        NULL
      }
    )
    dmt_tmp$LOC_ID <- gsub(".csv", "", basename(dmt))
    return(dmt_tmp)
  }
  
  # Run the function
  dmt_list <- lapply(dmt, dmt_check)
  # Bind the dfs
  dmt_df <- do.call(rbind, dmt_list)
  
  # Remove the dmt_list variables and clean
  message("Clean up working memory.")
  rm(dmt_list)
  gc()
  
  # Check the dates
  if(!is.Date(start_date) | !is.Date(end_date)){
    start_date <- lubridate::mdy(start_date)
    if(leap_year(mdy(paste0("01-01", "-", str_extract(end_date, "\\d{4}"))))){
      message("The end_date falls on a leap year. Consider adjusting your date
              if your end date is the last day in February.")
      end_date <- lubridate::mdy(end_date) 
    } else {
      message("The end_date does not fall on a leap year. If you supplied an out-of-range date, adjustments will be made. Proceed with caution and verify.")
      if(str_detect(end_date, "02[[:punct:]]29|29[[:punct:]]-02")){
        message("Leap year date set for non-leap year, adjusting from 02-29 to 02-28...verify this is correct.")
        end_date <- str_replace(end_date, "29", "28")
      }
      end_date <- lubridate::mdy(end_date)
    }
  }
  
  message("Start filter and merge.")
  
  # Filter and summarise
  dmt_df <- dmt_df |>
    mutate(Date = as.Date(parse_date_time(paste(year, yday), orders = "yj"))) |> 
    rename_with(~ ifelse(str_detect(.x, "prcp"), "prcp", .x)) |> 
    rename_with(~ ifelse(str_detect(.x, "tmax"), "tmax", .x)) |> 
    rename_with(~ ifelse(str_detect(.x, "tmin"), "tmin", .x)) |> 
    filter(Date %within% (start_date %--% end_date)) |> 
    mutate(Winter_Year = as.character(ifelse(month(Date) >= 11, year(Date) + 1, year(Date)))) |> # variable for grouping winter months
    filter(month(Date) >= month(start_date) | month(Date) <= month(end_date)) |> # this constrains the months to be after or before the months that define start and ends (e.g., 1-12 would take all month but 4-12 would take Apr-Dec)
    filter(day(Date) >= day(start_date) | day(Date) <= day(end_date)) |> # day constraints
    group_by(LOC_ID, Winter_Year) |> 
    summarise(prcp = sum(prcp),
              tmin = mean(tmin),
              tmax = mean(tmax))
  
  # Merge with PFW data
  pfw_out <- pfw |> 
    select(-SUB_ID, -Month, -Day, -Year, -Date) |> 
    distinct() |> 
    mutate(Winter_Year = gsub("PFW_", "", PROJ_PERIOD_ID)) |> 
    left_join(dmt_df)
  
  return(pfw_out) 
  
}

# Run function ----
pfw_dmt_site <- dmt_link_site(pfw = pfws,
                              dmt = dmt.files,
                              start_date = "12-01-1987",
                              end_date = "02-29-2024")

# Write the data ----
write_csv(pfw_dmt_site, file = args[2])
message("Script finished. Files written to: ", output_file)
