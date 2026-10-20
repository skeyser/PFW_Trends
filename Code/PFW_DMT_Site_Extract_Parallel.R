# **********************************************************
# Script Info ------------------------------------------------
#
# Script name:
#
# Script purpose:
#
# Author: Spencer R Keyser
#
# Date Created: 2026-10-20
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

# Package Loading -------------------------------------------
renv::load()

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(lubridate)
  library(stringr)
  library(future)
  library(future.apply)
})

# ***********************************************************
#
# #### Command line arguments
#
# ***********************************************************

# Example run-- Rscript PFW_DMT_Site_Extract_Parallel.R PFW_All_SurveyData.csv [100] [4]

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    paste(
      "Usage:",
      "Rscript PFW_DMT_Site_Extract.R",
      "<pfw.csv>",
      "<output.csv>",
      "[chunk_size]",
      "[n_workers]"
    ),
    call. = FALSE
  )
}

pfws_input <- args[1]
output_file <- args[2]

# Files per chunk
chunk_size <- if (length(args) >= 3) as.integer(args[3]) else 100

# Default will be 4 of 1 less than available workers unless specified in CLI
default_workers <- max(1, min(4, future::availableCores() - 1))
n_workers <- if (length(args) >= 4) as.integer(args[4]) else default_workers

if (is.na(chunk_size) || chunk_size < 1) {
  stop("chunk_size must be a positive integer.", call. = FALSE)
}

if (is.na(n_workers) || n_workers < 1) {
  stop("n_workers must be a positive integer.", call. = FALSE)
}

# ***********************************************************
#
# #### Get the paths for the script
#
# ***********************************************************
dmt_path <- "/Users/srk252/Documents/GIS/Env_Data/Daymet/"

if (!file.exists(pfws_input)) {
  stop("Input file does not exist: ", pfws_input, call. = FALSE)
}

if (!dir.exists(dmt_path)) {
  stop("Daymet directory does not exist: ", dmt_path, call. = FALSE)
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)


# ***********************************************************
#
# #### Daymet files
#
# ***********************************************************
dmt.files <- list.files(
  dmt_path,
  pattern = "^L[0-9]+\\.csv$",
  full.names = TRUE
)

##############################
######### Test ###############
# dmt.files <- dmt.files[1:100]

if (length(dmt.files) == 0) {
  stop(
    "No Daymet files found in: ",
    dmt_path,
    call. = FALSE
  )
}

message("Found ", length(dmt.files), " Daymet files.")

# ***********************************************************
#
# #### Dates
#
# ***********************************************************
start_date <- lubridate::mdy("12-01-1987")
end_date <- lubridate::mdy("02-29-2024")

if (is.na(start_date) || is.na(end_date)) {
  stop("Could not parse start_date or end_date.", call. = FALSE)
}

# ***********************************************************
#
# #### Single daymet file processing function
#
# ***********************************************************
process_one_daymet_file <- function(f, start_date, end_date) {
  
  loc_id <- sub("\\.csv$", "", basename(f))
  
  # Check the file has data
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
  
  # Read the daymet file
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
  
  # Process the file with the filters and summaries for the variables of interest
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
    ) |>
    mutate(PROJ_PERIOD_ID = paste0("PFW_", as.character(Winter_Year)),
           Winter_Year = as.character(Winter_Year)) |> 
    relocate(PROJ_PERIOD_ID, .before = prcp) |> 
    relocate(Winter_Year, .before = PROJ_PERIOD_ID)
}


# ***********************************************************
#
# #### Daymet processing by chunk with temp file writing
#
# ***********************************************************
process_daymet_chunk <- function(file_chunk,
                                 chunk_id,
                                 n_chunks,
                                 start_date,
                                 end_date,
                                 chunk_output_dir) {
  
  # Add a check on data.table to ensure limited parallel operations
  # internally
  data.table::setDTthreads(1)
  
  # Process the files per chunk
  chunk_results <- lapply(
    file_chunk,
    process_one_daymet_file,
    start_date = start_date,
    end_date = end_date
  )
  
  # Get the summaries per chunk
  chunk_summary <- data.table::rbindlist(
    chunk_results,
    use.names = TRUE,
    fill = TRUE
  )
  
  # File path per chunk id
  chunk_output_file <- file.path(
    chunk_output_dir,
    sprintf("daymet_summary_chunk_%05d.csv", chunk_id)
  )
  
  # Always write a valid CSV, even if no usable records occur in the chunk.
  if (nrow(chunk_summary) > 0) {
    data.table::fwrite(
      chunk_summary,
      chunk_output_file
    )
  } else {
    data.table::fwrite(
      data.table::data.table(
        LOC_ID = character(),
        Winter_Year = character(),
        PROJ_PERIOD_ID = character(),
        prcp = numeric(),
        tmin = numeric(),
        tmax = numeric()
      ),
      chunk_output_file
    )
  }
  
  rm(chunk_results, chunk_summary)
  gc(FALSE)
  
  message(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | Completed chunk ", chunk_id,
    " of ", n_chunks,
    " (", length(file_chunk), " files)"
  )
  
  return(chunk_output_file)
}

# ***********************************************************
#
# #### Split the Daymet files into appropriate chunks
#
# ***********************************************************
file_chunks <- split(
  dmt.files,
  ceiling(seq_along(dmt.files) / chunk_size)
)

n_chunks <- length(file_chunks)

message("Chunk size: ", chunk_size, " files")
message("Number of chunks: ", n_chunks)
message("Parallel workers: ", n_workers)

# Create a temp directory with files per chunk
chunk_output_dir <- file.path(
  tempdir(),
  paste0("daymet_chunk_summaries_", Sys.getpid())
)

message("Temporary chunk directory: ", chunk_output_dir)

dir.create(
  chunk_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Run the function in parallel using future.apply
message("Starting parallel Daymet processing...")

future::plan(
  future::multisession,
  workers = n_workers
)

chunk_files <- tryCatch(
  
  future.apply::future_lapply(
    seq_along(file_chunks),
    FUN = function(i) {
      
      process_daymet_chunk(
        file_chunk = file_chunks[[i]],
        chunk_id = i,
        n_chunks = n_chunks,
        start_date = start_date,
        end_date = end_date,
        chunk_output_dir = chunk_output_dir
      )
    },
    
    # Roughly two futures per worker
    future.scheduling = 2,
    
    future.seed = TRUE,
    
    future.stdout = TRUE,
    future.conditions = "condition"
    
  ),
  
  finally = {
    future::plan(future::sequential)
  }
)

message("Finished processing all Daymet chunks.")

# ***********************************************************
#
# #### Combine files from temp dir into output
#
# ***********************************************************
message("Combining chunk summaries...")

chunk_files <- unlist(chunk_files)

chunk_files <- chunk_files[
  file.exists(chunk_files) &
    file.info(chunk_files)$size > 0
]

if (length(chunk_files) == 0) {
  stop(
    "No chunk summary files were created.",
    call. = FALSE
  )
}

dmt_summary <- data.table::rbindlist(
  lapply(chunk_files, data.table::fread),
  use.names = TRUE,
  fill = TRUE
)

if (nrow(dmt_summary) == 0) {
  stop(
    "No usable Daymet observations were produced after filtering.",
    call. = FALSE
  )
}

message(
  "Combined Daymet summary rows: ",
  nrow(dmt_summary)
)

# ***********************************************************
#
# #### Load PFW data
#
# ***********************************************************
message("Reading PFW data...")

pfws <- data.table::fread(pfws_input) |>
  as.data.frame()

required_pfw_cols <- c("LOC_ID", "PROJ_PERIOD_ID", "Year", "Month", "Day")

missing_pfw_cols <- setdiff(required_pfw_cols, names(pfws))

if (length(missing_pfw_cols) > 0) {
  stop(
    "PFW input is missing required column(s): ",
    paste(missing_pfw_cols, collapse = ", "),
    call. = FALSE
  )
}

pfws <- pfws |>
  mutate(
    Month = stringr::str_pad(as.character(Month), width = 2, pad = "0"),
    Day = stringr::str_pad(as.character(Day), width = 2, pad = "0"),
    Date = lubridate::ymd(paste(Year, Month, Day, sep = "-"))
  )

message("Loaded PFW data from: ", pfws_input)
message("PFW rows: ", nrow(pfws))
message("Unique PFW sites: ", dplyr::n_distinct(pfws$LOC_ID))

# ***********************************************************
#
# #### Join with PFW Site date for final output
#
# ***********************************************************
message("Joining Daymet summaries to PFW data...")

# Currently formatted to return winter summaries for each
# year x site...to pair with only sites surveyed remove
# PROJ_PERIOD_ID from select and add to join by argument
pfw_dmt_site <- pfws |>
  #select(-any_of(c("SUB_ID", "Month", "Day", "Year", "Date", "PROJ_PERIOD_ID"))) |>
  select(LOC_ID, LATITUDE, LONGITUDE, SUBNATIONAL1_CODE) |> 
  distinct() |>
  left_join(
    as.data.frame(dmt_summary),
    by = c("LOC_ID")
  )

# ***********************************************************
#
# #### Write the output to file
#
# ***********************************************************
data.table::fwrite(
  pfw_dmt_site,
  output_file
)

message("Final output written to: ", output_file)

# ***********************************************************
#
# #### Clean everything up
#
# ***********************************************************
unlink(chunk_output_dir, recursive = TRUE)

message("Script finished successfully.")
