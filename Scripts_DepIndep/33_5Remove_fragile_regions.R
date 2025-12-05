#!/usr/bin/env Rscript

# nohup Rscript 33_5Remove_fragile_regions.R 3 5 0 Dep UM_TKO_CIT_Duv_CLEAN nmutations CH12F3 > nohup_5Remove_fragile_regionsUMTKODuvnmut.out 2>&1 &
# nohup Rscript 33_5Remove_fragile_regions.R 3 5 0 Dep UM_TKO_CIT_Duv_CLEAN density CH12F3 > nohup_5Remove_fragile_regionsUMTKODuvdens.out 2>&1 &

library(GenomicRanges)
library(readr)
library(dplyr)
library(stringr)

##########################
#       FUNCTIONS        #
##########################

# Function to log messages (optional for better debugging)
log_message <- function(type, message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("%s [%s] %s\n", timestamp, toupper(type), message))
}

# Reads a BED file into a GRanges object
# Create a GRanges object from a dataframe with columns chrom, start, end
# Assumes 'df' is already read and contains at least these three columns
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


##########################
#          MAIN          #
##########################


main <- function(args) {
  if (length(args) < 7) {
    stop("Error: At least 7 parameters required: threshold_num, min_count, enlargement, type (Dep or Indep), sample, dens_mut, cell_line.")
  }

  threshold_num <- as.integer(args[1])
  min_count <- as.integer(args[2])
  enlargement <- as.integer(args[3])
  type <- args[4]
  sample <- args[5]
  dens_mut <- args[6]
  cell_line <- args[7]

  base_dir <- file.path("./Dep_Indep_dataset",
                        paste0("Threshold_", threshold_num),
                        paste0("Enlargement_", enlargement),
                        sample,
                        paste0("Ranking_", type, "_hotspots"),
                        dens_mut)

  input_filename <- sprintf("%d_%s_%s_%s_Ranked_hotspots_%s_%d%d.tsv",
                            enlargement, sample, type, cell_line, dens_mut, threshold_num, min_count)
  input_filepath <- file.path(base_dir, input_filename)
  blacklist_filepath <- "./Controls/mm9-blacklist.bed"

  if (!file.exists(input_filepath)) {
    stop(paste("Input file does not exist:", input_filepath))
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
  removed_df  <- input_df[to_remove, ]

  # Save filtered rows
  output_filepath <- file.path(base_dir,
                               sprintf("%d_%s_%s_%s_Ranked_hotspots_%s_%d%d_blackfilt.tsv",
                                       enlargement, sample, type, cell_line, dens_mut, threshold_num, min_count))

  log_message("INFO", paste("Writing filtered output to:", output_filepath))
  write_tsv(filtered_df, output_filepath)

  # Saving deleted rows (blacklist)
  removed_filepath <- file.path(base_dir,
                                sprintf("%d_%s_%s_%s_Deleted_hotspots_%s_%d%d_blacklist.tsv",
                                        enlargement, sample, type, cell_line, dens_mut, threshold_num, min_count))
  log_message("INFO", paste("Writing overlapping (blacklisted) regions to:", removed_filepath))
  if (nrow(removed_df) > 0) {
  write_tsv(removed_df, removed_filepath)
  log_message("INFO", paste("Blacklisted regions saved to:", removed_filepath))
} else {
  log_message("INFO", "No overlapping regions found, no blacklist file written")
}


  log_message("INFO", "Process completed successfully")


}

args <- commandArgs(trailingOnly = TRUE)
main(args)

