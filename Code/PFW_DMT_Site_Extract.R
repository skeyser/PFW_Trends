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
    "Usage: Rscript PFW_DMT_Site_Extract.R <pfw.csv> <output dir> [chunk_size] [start_date] [end_date]",
    call. = FALSE
  )
}

pfws_input <- args[1]
# Read in the survey data
#pfws_input <- "Data/PFW_All_SurveyData.csv"
dmt_path <- "/Users/srk252/Documents/GIS/Env_Data/Daymet/"
output_file <- args[2]

# Options for chunk_size, start_date and end_date passed via CLI
if(length(args) < 3) print("No chunk size given...proceeding file by file.")
cs <- if(length(args) >= 3) as.integer(args[3]) else NULL
sdate <- if(length(args) >= 4) as.character(args[4]) else "12-01-1987"
edate <- if(length(args) >= 5) as.character(args[5]) else "02-29-2024"

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
                          start_date,
                          end_date,
                          chunk_size,
                          temp_dir = tempdir()){
  
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
  
  # ***********************************************************
  #
  # #### Process one daymet file at a time function
  #
  # ***********************************************************
  process_one_daymet_file <- function(f) {
    
    loc_id <- sub("\\.csv$", "", basename(f))
    
    # Check the first line to confirm data exists
    first_data_line <- tryCatch(
      readr::read_lines(
        f,
        skip = 6,
        n_max = 1,
        progress = FALSE
      ),
      error = function(e) character(0)
    )
    
    if (length(first_data_line) == 0 || !nzchar(trimws(first_data_line))) {
      message("Site ", loc_id, " has an empty Daymet file; skipping.")
      return(NULL)
    }
    
    # Read in the file with tryCatch
    d <- tryCatch(
      readr::read_csv(
        f,
        skip = 6,
        show_col_types = FALSE,
        progress = FALSE
      ),
      error = function(e) {
        message(
          "Could not read Daymet file for site ", loc_id,
          ": ", conditionMessage(e)
        )
        NULL
      }
    )
    
    if (is.null(d) || nrow(d) == 0) {
      return(NULL)
    }
    
    # Summarize immediately file by file
    d |>
      mutate(Date = as.Date(parse_date_time(paste(year, yday), orders = "yj"))) |> 
      rename_with(~ ifelse(str_detect(.x, "prcp"), "prcp", .x)) |> 
      rename_with(~ ifelse(str_detect(.x, "tmax"), "tmax", .x)) |> 
      rename_with(~ ifelse(str_detect(.x, "tmin"), "tmin", .x)) |> 
      dplyr::mutate(
        LOC_ID = loc_id,
        Date = lubridate::make_date(year, 1, 1) +
          lubridate::days(yday - 1)
      ) |>
      dplyr::filter(
        Date >= start_date,
        Date <= end_date
      ) |>
      dplyr::mutate(
        Winter_Year = as.character(
          ifelse(
            lubridate::month(Date) >= 11,
            lubridate::year(Date) + 1,
            lubridate::year(Date)
          )
        )
      ) |>
      dplyr::filter(
        lubridate::month(Date) >= lubridate::month(start_date) |
          lubridate::month(Date) <= lubridate::month(end_date)
      ) |>
      dplyr::filter(
        lubridate::day(Date) >= lubridate::day(start_date) |
          lubridate::day(Date) <= lubridate::day(end_date)
      ) |>
      dplyr::group_by(LOC_ID, Winter_Year) |>
      dplyr::summarise(
        prcp = sum(prcp, na.rm = TRUE),
        tmin = mean(tmin, na.rm = TRUE),
        tmax = mean(tmax, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  # ***********************************************************
  #
  # #### Check chunk_size
  #
  # ***********************************************************
  if(!is.null(chunk_size)){
    # Write chunk-level summaries to disk
    temp_summary_file <- file.path(temp_dir, "daymet_site_winter_summary.csv")
    
    if (file.exists(temp_summary_file)) {
      file.remove(temp_summary_file)
    }
    
    n_chunks <- ceiling(length(dmt) / chunk_size)
    
    message(
      "Processing ", length(dmt), " files in ",
      n_chunks, " chunks of up to ", chunk_size, " files."
    )
    
    for (i in seq_len(n_chunks)) {
      
      start_idx <- (i - 1) * chunk_size + 1
      end_idx <- min(i * chunk_size, length(dmt))
      
      file_chunk <- dmt[start_idx:end_idx]
      
      message(
        "Processing chunk ", i, " of ", n_chunks,
        " (files ", start_idx, " through ", end_idx, ")"
      )
      
      # Each element is now only a small summary table, not daily observations.
      chunk_results <- lapply(file_chunk, process_one_daymet_file)
      
      chunk_summary <- data.table::rbindlist(
        chunk_results,
        use.names = TRUE,
        fill = TRUE
      )
      
      if (nrow(chunk_summary) > 0) {
        file_exists <- file.exists(temp_summary_file)
        
        data.table::fwrite(
          chunk_summary,
          temp_summary_file,
          append = file_exists,
          col.names = !file_exists
        )
      }
      
      rm(chunk_results, chunk_summary)
      gc(FALSE)
    }
  } else {
    message("Chunk size not detected...proceeding for one by one file processing.")
    # Do the full one by one run
    dmt_summary <- lapply(dmt, process_one_daymet_file)
    dmt_summary <- tibble::as_tibble(dmt_summary)
  }
  
  if(!is.null(chunk_size) & file.exists(temp_summary_file)){
    # Read the temp summary file if run in chunks
    dmt_summary <- data.table::fread(temp_summary_file) |>
      tibble::as_tibble()
    
    # Clean up after loading the data
    unlink(temp_summary_file)
    
  }
  
  message("Merge DMT summaries w/ PFW data.")
  
  # Merge with PFW data
  pfw_out <- pfw |>
    select(-SUB_ID, -Month, -Day, -Year, -Date) |>
    distinct() |>
    mutate(Winter_Year = gsub("PFW_", "", PROJ_PERIOD_ID)) |>
    left_join(dmt_summary |> mutate(Winter_Year = as.character(Winter_Year)))
  
  return(pfw_out)
  
}

# Run function ----
pfw_dmt_site <- dmt_link_site(pfw = pfws,
                              dmt = dmt.files,
                              start_date = sdate,
                              end_date = edate,
                              chunk_size = cs)

# Write the data ----
write_csv(pfw_dmt_site, file = args[2])
message("Script finished. Files written to: ", output_file)
