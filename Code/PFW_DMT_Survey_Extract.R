## Script Info ------------------------------------------------
##
## Script name: PFW Daymet extractions 2-day interval
##
## Script purpose:
##
## Author: Spencer R Keyser
##
## Date Created: 2026-10-15
##
## Email: srk252@cornell.edu
##
## Github: https://github.com/skeyser
##
##
##
## Notes:
##
##
##

## Defaults
options(scipen = 6, digits = 4)

renv::load()

## Package Loading -------------------------------------------
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
    "Usage: Rscript PFW_DMT_Extract.R <pfw.csv> <output dir>",
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

# ***********************************************************
#
# PFW x Daymet Function ####
#
# ***********************************************************
## ***********************************************************
##
## Section Notes: PFW site-level information
## Daymet data is location specific and includes 
## all years and all days
##
## ***********************************************************
dmt_link_survey <- function(pfw,
                     dmt,
                     days = 2,
                     save_progress = F,
                     save_interval = 10000){
  usites <- unique(pfw$LOC_ID)
  
  s <- 0
  clim_out <- vector(mode = "list", length = nrow(pfw))
  
  message("Starting site loop.")
  last_time <- Sys.time()
  for(site in usites){
    s <- s + 1
    if (s %% 1000 == 0) {
      current_time <- Sys.time()
      
      # Diff time every 1000th iteration
      time_diff <- difftime(current_time, last_time, units = "mins")
      
      # Message for tracking
      cat(sprintf("\rProcessing Site: %d / %d | Last 1,000 steps: %.2f mins", 
                  s, length(usites), time_diff))
      utils::flush.console()
      
      # Update tracker
      last_time <- current_time
    }
    # Filter and structure date for PFW survey at location site
    pfw_tmp <- pfw |> 
      filter(LOC_ID == site)
    
    # Find the matching file
    matched_file <- dmt[str_equal(gsub(".csv", "", basename(dmt)), site)]
    
    # Checks
    if (length(matched_file) == 0) {
      message("No Daymet file found for site ", site, "...skipping...")
      next
    }
    
    # Checks
    if (length(matched_file) > 1) {
      message("Multiple Daymet files found for site ", site, "...skipping...")
      next
    }
    
    # Check whether the file has anything beyond the 6 metadata lines
    file_lines <- readr::read_lines(matched_file, progress = FALSE)
    
    if (length(file_lines) <= 6 || !any(nzchar(trimws(file_lines[-seq_len(6)])))) {
      message("Site ", site, " has an empty Daymet file...skipping...")
      next
    }
    dmt_tmp <- tryCatch(
      readr::read_csv(
        matched_file,
        skip = 6,
        show_col_types = FALSE,
        progress = FALSE
      ),
      error = function(e) {
        message(
          "Could not read Daymet file for site ", site,
          ": ", conditionMessage(e),
          "...skipping..."
        )
        NULL
      }
    )
    
    if (is.null(dmt_tmp) || ncol(dmt_tmp) == 0 || nrow(dmt_tmp) == 0) {
      message("Site ", site, " has no usable Daymet data...skipping...")
      next
    }
    
    dmt_tmp <- dmt_tmp |> 
      mutate(Date = as.Date(parse_date_time(paste(year, yday), orders = "yj"))) |> 
      rename_with(~ ifelse(str_detect(.x, "prcp"), "prcp", .x)) |> 
      rename_with(~ ifelse(str_detect(.x, "tmax"), "tmax", .x)) |> 
      rename_with(~ ifelse(str_detect(.x, "tmin"), "tmin", .x))
    
    # We need to extract and summarise the mean clim for a survey on interval days
    date_pairs <- data.frame(sdate = pfw_tmp$Date, edate = (pfw_tmp$Date + (days - 1)))
    
    climpair_out <- vector(mode = "list", length = nrow(date_pairs))
    
    # Loop for summarizing the data
    for(i in 1:nrow(date_pairs)){
      dt.tmp <- date_pairs[i,]
      clim_pair <- dmt_tmp |> 
        filter(Date %within% lubridate::interval(dt.tmp$sdate, dt.tmp$edate)) |> 
        summarise(year = min(year),
                  yday = min(yday),
                  prcp = mean(prcp),
                  tmax = mean(tmax),
                  tmin = mean(tmin),
                  Date = min(Date)) |> 
        select(Date, prcp, tmin, tmax)
      climpair_out[[i]] <- clim_pair
    }
    
    climpair_out <- do.call(rbind, climpair_out)
    climpair_out$LOC_ID <- site
    pfw_tmp <- left_join(pfw_tmp, climpair_out, relationship = "many-to-many", by = c("LOC_ID", "Date"))
    
    clim_out[[s]] <- climpair_out
    
    ## Option to save progress intermittently
    if(isTRUE(save_progress)){
      if(s %% save_interval == 0){
        clim_tmp_out <- do.call(rbind, clim_out)
        write.csv(clim_tmp_out, file = args[2])
        cat("Saved interim files to: ", output_file)
      }
    }
    
  }
  
  clim_out <- do.call(rbind, clim_out)
  return(clim_out)
}


# Run the function ----
pfw_dmt_survey <- dmt_link(pfw = pfws,
                    dmt = dmt.files,
                    days = 2,
                    save_progress = F,
                    save_interval = 10000)


# Merge w/ PFW ----
pfw_dmt_m <- left_join(pfws, pfw_dmt, by = c("LOC_ID", "Date"))

# Write the data ----
write_csv(pfw_dmt_m, file = args[2])
message("Script finished. Files written to: ", output_file)


# ## Code below to isolating some problematic Daymet files
# check_dmt <- function(i, files){
#   file_path <- files[i]
#   
#   tryCatch({
#     d <- read_csv(file_path, skip = 6, n_max = 1, show_col_types = F, progress = F)
#     
#     if(ncol(d) < 5){
# 
#       message(paste0("Error at Index [", i, "]: ", file_path)) 
#       stop("File has less than 5 columns.")
#     }
#     
#     return(TRUE) 
#     
#   }, error = function(e) {
# 
#     return(paste("Failed:", e$message))
#   })
# }
# 
# results <- lapply(seq_along(dmt.files), check_dmt, files = dmt.files)
# 
# bad_files <- basename(dmt.files[which(unlist(results != "TRUE"))])
# 
# saveRDS(bad_files, file = here("Data/Daymet_Bad_Extract.rds"))
