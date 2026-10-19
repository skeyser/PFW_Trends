## ----------------------------------------------------------
##
## Script name: PFW Functions
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

## -----------------------------------------------------------

## Zero fill w/ function from Wes H.
ZeroFillPFW_Wes <- function(InputData, SpeciesCodes, rollup = TRUE){
  
  #The first three steps collect information that
  # will be necessary for taxonomic roll-up and zero-filling. First, we create
  # a table that contains only the information describing the observation 
  # periods (and not anything specific to any single birds species). Second,
  # we need to count the number of bird species for which zero-filling has been
  # requested, in order to know how many times we will need to loop through the
  # zero-filling process. Third, we need to check the species codes for which 
  # zero-filling was requested by comparing these to a list of all species codes 
  # found in InputData.
  
  #Create the table of sampling event information from all observation periods
  PFW_SED <- InputData %>%
    select(-OBS_ID, -SPECIES_CODE, -alt_full_spp_code, -HOW_MANY, -PLUS_CODE, 
           -VALID, -REVIEWED) %>%
    distinct()
  
  #Count the number of species codes
  NSpecies <- length(SpeciesCodes)
  
  #Identify any species codes that are not found in the data
  InvalidSppCodes <- InputData %>%
    select(SPECIES_CODE) %>%
    distinct() %>%
    anti_join(y = ., x = tibble(SPECIES_CODE = SpeciesCodes), by = "SPECIES_CODE")
  #For each non-existing species code, print an error message
  if (nrow(InvalidSppCodes) >= 1){
    for (i in 1:nrow(InvalidSppCodes)){
      cat(paste("THE SPECIES CODE \'", 
                InvalidSppCodes$SPECIES_CODE[i], 
                "\' DOES NOT EXIST\n",
                sep = ""))
      cat(" IN THE DATA TABLE, AND ZERO-FILLING IS NOT POSSIBLE.\n\n")
    }
    #break out of the function by returning a null result
    cat("THIS FUNCTION IS STOPPING BECAUSE ZERO-FILLING OF ALL\n")
    cat(" REQUESTED SPECIES CODES IS NOT POSSIBLE.\n\n")
    return(NULL)
  }
  
  #Loop through the species codes in turn, possibly doing a taxonomic roll-up,
  # and always zero-filling each species' data.
  for (i in 1:NSpecies){
    #
    #Check if this species code is "nobird" and if so skip over this 
    # species code.
    if (SpeciesCodes[i] == "nobird"){
      cat("THE SPECIES CODE \'nobird\' INDICATES AN ABSENCE OF BIRDS,\n")
      cat(" AND IS NOT AN ACTUAL BIRD SPECIES. NO ZERO-FILLING WILL\n")
      cat("BE DONE FOR THE SPECIES_CODE \'nobird\'.\n\n")
      next
    }
    #Place the SPECIES_CODE for the ith species into Spp_to_Zero_Fill.
    Spp_to_Zero_Fill <- SpeciesCodes[i]
    #
    #If a taxonomic roll-up is to be done, then do so for the focal species.
    if (rollup){
      InputData <- InputData %>%
        mutate(SPECIES_CODE = case_when(
          (is.na(alt_full_spp_code)
           ~ SPECIES_CODE),
          ((!is.na(alt_full_spp_code) &
              alt_full_spp_code == Spp_to_Zero_Fill)
           ~ alt_full_spp_code),
          ((!is.na(alt_full_spp_code) &
              alt_full_spp_code != Spp_to_Zero_Fill)
           ~ SPECIES_CODE))
        )
    }
    
    #Create template for zero-count records for the focal species.
    FakeSppSpecificFields <- as_tibble_row(list(OBS_ID = "OBSnull",
                                                SPECIES_CODE = Spp_to_Zero_Fill,
                                                HOW_MANY = 0,
                                                PLUS_CODE = NA,
                                                VALID = 1,
                                                REVIEWED = 0))
    
    #Zero-fill the data for the focal species.
    ZeroFilledFocalSpecies <- InputData %>%
      filter(SPECIES_CODE == Spp_to_Zero_Fill) %>%
      bind_rows(., bind_cols(PFW_SED, FakeSppSpecificFields)) %>%
      group_by(across(c(names(PFW_SED), "SPECIES_CODE"))) %>%
      summarize(OBS_ID = min(OBS_ID), HOW_MANY = sum(HOW_MANY), 
                PLUS_CODE = max(PLUS_CODE),
                VALID = min(VALID), REVIEWED = max(REVIEWED)) %>%
      ungroup()
    ZeroFilledFocalSpecies$PLUS_CODE <- as.logical(ZeroFilledFocalSpecies$PLUS_CODE)
    
    #Put the zero-filled data into a table that accumulates data from all of
    # the species named in the vector SpeciesCodes. Either create the tibble
    # to be returned if it does not already exist, or append rows to an
    # existing object.
    if ((exists("OutputData")))
      OutputData <- bind_rows(OutputData, ZeroFilledFocalSpecies)
    else
      OutputData <- ZeroFilledFocalSpecies
  }
  #After all species' data are zero-filled, return the zero-filled data to the 
  # user. Unless the user of the function directs the output into some object,
  # the output will be printed to the console and not be saved for subsequent use.
  #Also, if a taxonomic roll-up has been done, print a message to remind users
  # this was done.
  if (rollup){
    cat("\nNOTE: ANY RECORDS OF BIRDS REPORTED TO THE LEVEL OF A SUBSPECIES OR\n")
    cat(" OTHER RECOGNIZABLE FORM HAVE BEEN TREATED AS MEMBERS OF THE FULL\n")
    cat(" SPECIES FOR THE PURPOSE OF ZERO-FILLING, AND 'species_code' VALUES\n")
    cat(" FOR THESE RECOGNIZABLE FORMS HAVE BEEN REPLACE BY THE CODES FOR THE\n")
    cat(" FULL SPECIES IN THE OUTPUT.\n\n")
  }
  return(OutputData)
}

## Second version of Zero-filling but putting the data in a wide format
## and allowing for "All_Sp" selection
ZeroFillPFW <- function(DAT,
                        SpeciesCodes = "All_Sp",
                        rollup = FALSE) {
  
  require(dplyr)
  require(tidyr)
  
  ## Helper functions
  
  ## Take first non-NA
  first_non_missing <- function(x) {
    x[which(!is.na(x))[1]]
  }
  
  ## Max Count for non-NA
  max_or_na <- function(x) {
    if (all(is.na(x))) {
      x[NA_integer_]
    } else {
      max(x, na.rm = TRUE)
    }
  }
  
  ## Sum count for non-NA
  sum_or_na <- function(x) {
    if (all(is.na(x))) {
      x[NA_integer_]
    } else {
      sum(x, na.rm = TRUE)
    }
  }
  
  ## Columns we expect at a minimum from PFW
  required_cols <- c(
    "SUB_ID", "SPECIES_CODE", "HOW_MANY",
    "VALID", "REVIEWED", "alt_full_spp_code"
  )
  
  ## Check the missing ones
  missing_cols <- setdiff(required_cols, names(DAT))
  
  ## Have to have those
  if (length(missing_cols) > 0) {
    stop(
      "DAT is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  ## If Tax rollup then we replace the SPECIES_CODE with the
  ## alt_full_spp_code
  if (rollup) {
    DAT <- DAT %>%
      mutate(
        SPECIES_CODE = coalesce(alt_full_spp_code, SPECIES_CODE)
      )
  }
  
  ## Find the species we want to ZF, either All_Sp or some provided vec
  if (length(SpeciesCodes) == 1 && SpeciesCodes == "All_Sp") {
    
    SpeciesCodes <- DAT %>%
      distinct(SPECIES_CODE) %>%
      pull(SPECIES_CODE) %>%
      sort()
    
    message("Zero-filling all species found in DAT...")
    
  } else {
    
    SpeciesCodes <- sort(unique(SpeciesCodes))
    
    message("Zero-filling supplied species codes...")
  }
  
  ## "nobird" is not a species count column
  SpeciesCodes <- setdiff(SpeciesCodes, "nobird")
  
  ## Get OBS-level column names
  ## Check against the DAT 
  observation_cols <- c(
    "OBS_ID",
    "SPECIES_CODE",
    "alt_full_spp_code",
    "HOW_MANY",
    "VALID",
    "REVIEWED"
  )
  
  metadata_cols <- setdiff(names(DAT), observation_cols)
  
  ## SUB_ID is the event ID, so it should not be summarized inside across()
  metadata_cols_no_subid <- setdiff(metadata_cols, "SUB_ID")
  
  ## Get one value per sampling event
  ## Max for the REVIEW and VALID Data
  event_metadata <- DAT %>%
    group_by(SUB_ID) %>%
    summarise(
      across(
        all_of(metadata_cols_no_subid),
        first_non_missing
      ),
      VALID = max_or_na(VALID),
      REVIEWED = max_or_na(REVIEWED),
      .groups = "drop"
    )
  
  ## Find species that were reviewed and store them
  reviewed_species <- DAT %>%
    filter(
      VALID == 1,
      REVIEWED == 1,
      !is.na(HOW_MANY),
      SPECIES_CODE %in% SpeciesCodes
    ) %>%
    distinct(SUB_ID, SPECIES_CODE) %>%
    group_by(SUB_ID) %>%
    summarise(
      Rspec = paste(sort(SPECIES_CODE), collapse = "|"),
      .groups = "drop"
    )
  
  ## Get the max count at the sampling event per species
  species_counts <- DAT %>%
    filter(SPECIES_CODE %in% SpeciesCodes) %>%
    group_by(SUB_ID, SPECIES_CODE) %>%
    summarise(
      HOW_MANY = max_or_na(HOW_MANY),
      .groups = "drop"
    )
  
  ## Expand the grid so we have Species x Sampling matrix
  complete_counts <- expand_grid(
    SUB_ID = event_metadata$SUB_ID,
    SPECIES_CODE = SpeciesCodes
  ) %>%
    left_join(
      species_counts,
      by = c("SUB_ID", "SPECIES_CODE")
    ) %>%
    mutate(
      ## Zero-fill
      HOW_MANY = coalesce(HOW_MANY, 0L)
    )
  
  ## Put this is a wide format
  species_wide <- complete_counts %>%
    pivot_wider(
      names_from = SPECIES_CODE,
      values_from = HOW_MANY,
      values_fill = 0L
    )
  
  ## Bundle the data
  ZF_out <- event_metadata %>%
    left_join(reviewed_species, by = "SUB_ID") %>%
    left_join(species_wide, by = "SUB_ID") %>%
    relocate(
      any_of(c(
        "VALID",
        "REVIEWED",
        "PLUS_CODE",
        "Rspec"
      )),
      .before = SpeciesCodes[1]
    )
  
  if (rollup) {
    message(
      "Taxonomic rollup applied: alt_full_spp_code was used where available."
    )
  }
  
  ZF_out
}
