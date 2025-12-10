#!/usr/bin/env Rscript

# nohup Rscript 05_Remove_fragile_sites.R CH12F3 3 5 0 UM_TKO_CIT_Duv_CLEAN --Dep --dens > nohup_05Remove_fragile_sites_Dep_dens.out 2>&1 &
# nohup Rscript 05_Remove_fragile_sites.R CH12F3 3 5 0 UM_TKO_CIT_Duv_CLEAN > nohup_05Remove_fragile_sites_DepIndep.out 2>&1 &

library(GenomicRanges)
library(readr)
library(dplyr)
library(stringr)

##########################
#       FUNCTIONS        #
##########################

# Parse command-line flags for Dep/Indep and mut/dens options
parse_flags <- function(flags) {

  RUN_DEP <- TRUE
  RUN_INDEP <- TRUE
  RUN_MUT <- TRUE
  RUN_DENS <- TRUE

  TYPE_SET <- FALSE
  MODE_SET <- FALSE
  ALL_SEEN <- FALSE

  if (length(flags) > 0) {
    flags <- flags[flags != ""]

    for (f in flags) {
      if (f == "--ALL") {
        ALL_SEEN <- TRUE
      }
      if (f == "--Dep") {
        RUN_DEP <- TRUE
        RUN_INDEP <- FALSE
        TYPE_SET <- TRUE
      }
      if (f == "--Indep") {
        RUN_DEP <- FALSE
        RUN_INDEP <- TRUE
        TYPE_SET <- TRUE
      }
      if (f == "--mut") {
        RUN_MUT <- TRUE
        RUN_DENS <- FALSE
        MODE_SET <- TRUE
      }
      if (f == "--dens") {
        RUN_MUT <- FALSE
        RUN_DENS <- TRUE
        MODE_SET <- TRUE
      }
    }

    # Apply ALL only to non-set categories
    if (ALL_SEEN && !TYPE_SET) {
      RUN_DEP <- TRUE
      RUN_INDEP <- TRUE
    }
    if (ALL_SEEN && !MODE_SET) {
      RUN_MUT <- TRUE
      RUN_DENS <- TRUE
    }
  }

  cat("[INFO] Flags parsed:\n")
  cat("       RUN_DEP =", RUN_DEP, "\n")
  cat("       RUN_INDEP =", RUN_INDEP, "\n")
  cat("       RUN_MUT =", RUN_MUT, "\n")
  cat("       RUN_DENS =", RUN_DENS, "\n")

  return(list(
    RUN_DEP=RUN_DEP,
    RUN_INDEP=RUN_INDEP,
    RUN_MUT=RUN_MUT,
    RUN_DENS=RUN_DENS
  ))
}

# Function to log messages (optional for better debugging)
log_message <- function(type, message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("%s [%s] %s\n", timestamp, toupper(type), message))
}

# Reads a BED file into a GRanges object
df_to_granges <- function(df) {
  # Check if required columns exist
  required_cols <- c("chrom", "start", "end")
  if (!all(required_cols %in% colnames(df))) {
    stop("Dataframe must contain columns: chrom, start, end")
  }

  # Convert start and end to numeric if needed
  df$start <- as.numeric(df$start)
  df$end <- as.numeric(df$end)
  if (any(is.na(df$start)) || any(is.na(df$end))) {
    stop("Start or end columns contain non-numeric values after coercion")
  }

  GRanges(seqnames = df$chrom, ranges = IRanges(start = df$start, end = df$end))
}


# Core filtering logic, moved to a function for the loop
filter_hotspots <- function(cell_line, threshold_num, min_count, enlargement, sample, type, dens_mut) {
    
    log_message("INFO", sprintf("Starting filtering for Type=%s, Mode=%s", type, dens_mut))

# Determine the directory name based on 'type' (Dep/Indep)
    type_dir_name <- if (type == "Dep") {
        "Dependent_hotspots"
    } else if (type == "Indep") {
        "Independent_hotspots"
    } else {
        stop("[ERROR] Invalid type specified.")
    }

    base_dir <- file.path("scratch/DepIndep_dataset",
                          paste0("Threshold_", threshold_num),
                          paste0("Enlargement_", enlargement),
                          type_dir_name,
                          sample,
                          paste0("Ranking_", type, "_hotspots"),
                          dens_mut)

    input_filename <- sprintf("%d_%s_%s_%s_Ranked_hotspots_%s_%d_%d.tsv",
                              enlargement, sample, type, cell_line, dens_mut, threshold_num, min_count)
    input_filepath <- file.path(base_dir, input_filename)
    blacklist_filepath <- "/scratch/reference_data/mm9-blacklist.bed"

    if (!file.exists(input_filepath)) {
      log_message("WARNING", paste("Input file does not exist, skipping:", input_filepath))
      return(NULL)
    }
    log_message("INFO", paste("Reading input TSV:", input_filepath))
    input_df <- read_tsv(input_filepath, col_types = cols())
    
    log_message("INFO", paste("Reading blacklist BED:", blacklist_filepath))
    blacklist_df <- read_tsv( blacklist_filepath, col_names = FALSE, col_types = cols(), col_select = 1:3)
    colnames(blacklist_df) <- c("chrom", "start", "end")


    log_message("INFO", paste("Generating GRanges objects:"))
    input_gr <- df_to_granges(input_df)
    blacklist_gr <- df_to_granges(blacklist_df)

    log_message("INFO", "Identifying overlapping regions...")
    hits <- findOverlaps(input_gr, blacklist_gr)
    to_remove <- unique(queryHits(hits)) # Unique indices of input_df to remove
    log_message("INFO", paste(length(to_remove), "overlapping regions to remove"))

    # filtered and deleted 
    filtered_df <- input_df[-to_remove, ]
    removed_df <- input_df[to_remove, ]

    # Save filtered rows
    output_dir <- file.path(base_dir, "Filtered_hotspots")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      log_message("INFO", paste("Created output directory:", output_dir))
    }
    output_filepath <- file.path(output_dir,
                                 sprintf("%d_%s_%s_%s_Ranked_hotspots_%s_%d%d_blackfilt.tsv",
                                         enlargement, sample, type, cell_line, dens_mut, threshold_num, min_count))

    log_message("INFO", paste("Writing filtered output to:", output_filepath))
    write_tsv(filtered_df, output_filepath)

    # Saving deleted rows (blacklist)
    removed_filepath <- file.path(output_dir,
                                  sprintf("%d_%s_%s_%s_Deleted_hotspots_%s_%d%d_blacklist.tsv",
                                          enlargement, sample, type, cell_line, dens_mut, threshold_num, min_count))
    log_message("INFO", paste("Writing overlapping (blacklisted) regions to:", removed_filepath))
    if (nrow(removed_df) > 0) {
      write_tsv(removed_df, removed_filepath)
      log_message("INFO", paste("Blacklisted regions saved to:", removed_filepath))
    } else {
      log_message("INFO", "No overlapping regions found, no blacklist file written")
    }
    
    log_message("INFO", sprintf("\nFiltering for Type=%s, Mode=%s completed.", type, dens_mut))
}


##########################
#          MAIN          #
##########################


main <- function(args) {
    # Expected Positional Arguments (5 total): cell_line, threshold_num, min_count, enlargement, sample
    # Remaining arguments are optional flags (e.g., --Dep, --dens)

    if (length(args) < 5) {
        stop("Error: At least 5 positional parameters required: cell_line, threshold_num, min_count, enlargement, sample. Optional flags follow.")
    }

    # Positional Arguments
    cell_line <- args[1]
    threshold_num <- as.integer(args[2])
    min_count <- as.integer(args[3])
    enlargement <- as.integer(args[4])
    sample <- args[5]

    # Optional Flags
    extra_flags <- if (length(args) >= 6) {
        args[6:length(args)]
    } else {
        character(0) #  Empty vector if no flags
    }

    flags <- parse_flags(extra_flags)

    # Build combinations list based on parsed flags
    combinations <- list()
    
    if (flags$RUN_DEP) {
        if (flags$RUN_DENS) combinations <- append(combinations, list(list(type="Dep", dens_mut="density")))
        if (flags$RUN_MUT) combinations <- append(combinations, list(list(type="Dep", dens_mut="nmutations")))
    }
    if (flags$RUN_INDEP) {
        if (flags$RUN_DENS) combinations <- append(combinations, list(list(type="Indep", dens_mut="density")))
        if (flags$RUN_MUT) combinations <- append(combinations, list(list(type="Indep", dens_mut="nmutations")))
    }
    
    if (length(combinations) == 0) {
        stop("[WARNING] No combinations selected to run. Check your flags.")
    }
    
    log_message("INFO", sprintf("Total combinations to process: %d", length(combinations)))

    # Loop over each combination and process
    for (combo in combinations) {
        type_i <- combo$type
        dens_mut_i <- combo$dens_mut
        
        filter_hotspots(cell_line, threshold_num, min_count, enlargement, sample, type_i, dens_mut_i)
    }

    log_message("INFO", "All filtering processes completed successfully.")
}

args <- commandArgs(trailingOnly = TRUE)
main(args)