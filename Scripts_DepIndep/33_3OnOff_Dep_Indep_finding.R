#!/usr/bin/env Rscript

# nohup Rscript 33_3OnOff_Dep_Indep_finding.R 3 5 0 Dep UM_TKO_CIT_Duv_CLEAN density CH12F3 > nohup_350Dep_UM_TKO_CIT_Duv_density.out 2>&1 &
# nohup Rscript 33_3OnOff_Dep_Indep_finding.R 3 5 0 Dep UM_TKO_CIT_Duv_CLEAN nmutations CH12F3 > nohup_350Dep_UM_TKO_CIT_Duv_nmutations.out 2>&1 &

# nohup Rscript 33_3OnOff_Dep_Indep_finding.R 3 5 0 Indep UM_TKO_CIT_Duv_CLEAN density CH12F3 > nohup_350Indep_UM_TKO_CIT_Duv_density.out 2>&1 &
# nohup Rscript 33_3OnOff_Dep_Indep_finding.R 3 5 0 Indep UM_TKO_CIT_Duv_CLEAN nmutations CH12F3 > nohup_350Indep_UM_TKO_CIT_Duv_nmutations.out 2>&1 &

library(GenomicRanges)
library(readr)

# Function to check file existence with error message
check_file_exists <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(sprintf("Error: file '%s' does not exist!", file_path))
  } else {
    cat(sprintf("File exists: %s\n", file_path))
  }
}

# Function to read BED-like file (chrom, start, end)
read_regions <- function(filepath) {
  check_file_exists(filepath)
  regions <- read_tsv(filepath, col_names = FALSE, col_types = cols(), col_select = 1:3)
  colnames(regions) <- c("chrom", "start", "end")
  return(regions)
}

# Function to create GRanges from dataframe with chrom, start, end
create_granges <- function(df) {
  return(GRanges(seqnames = df$chrom,
                 ranges = IRanges(start = df$start, end = df$end)))
}

# Function to find overlapping rows between input and target regions
get_overlapping_subset <- function(input_df, target_granges) {
  input_gr <- create_granges(input_df)
  overlaps <- findOverlaps(input_gr, target_granges)
  hits <- queryHits(overlaps)
  subset_df <- input_df[hits, ]
  return(subset_df)
}

# Main script function 
main <- function(args) {
  # Check correct number of arguments
  if (length(args) < 7) {
    stop("Error: At least 7 parameters required: threshold_num, min_count, enlargement, type (Dep or Indep), sample, dens_mut (density or nmutations), cell_line.")
  }
  
  # Parse arguments
  threshold_num <- as.integer(args[1])
  min_count <- as.integer(args[2])
  enlargement <- as.integer(args[3])
  type <- args[4]
  sample <- args[5]
  dens_mut <- args[6]
  cell_line <- args[7]
  
  # Build file paths
   base_dir <- file.path("./Dep_Indep_dataset",
                      paste0("Threshold_", threshold_num),
                      paste0("Enlargement_", enlargement),
                      sample,
                      paste0("Ranking_", type, "_hotspots"),
                      dens_mut)
  input_filename <- sprintf("%d_%s_%s_%s_Ranked_hotspots_%s_%d%d.tsv",
                            enlargement, sample, type, cell_line, dens_mut, threshold_num, min_count)
  input_filepath <- file.path(base_dir, input_filename)
  
  # Target region files
  on_target_file <- "./Controls/ON_target_CH12F3.bed"
  off_target_file <- "./Controls/OFF_target_CH12F3.bed"
  
  # Check files
  check_file_exists(input_filepath)
  check_file_exists(on_target_file)
  check_file_exists(off_target_file)
  
  # Read input and target regions
  input_data <- read_tsv(input_filepath, col_names = TRUE)
  on_target_regions <- read_regions(on_target_file)
  off_target_regions <- read_regions(off_target_file)
  
  # Create GRanges for target regions
  gr_on_target <- create_granges(on_target_regions)
  gr_off_target <- create_granges(off_target_regions)
  
  # Find overlaps and subset input data
  input_on_target <- get_overlapping_subset(input_data, gr_on_target)
  input_off_target <- get_overlapping_subset(input_data, gr_off_target)
  
  # Define output filenames
  output_on <- file.path(base_dir, sprintf("%d_%s_%s_%s_Ranked_hotspots_%s_ON_target_%d%d.tsv",
                                          enlargement, sample, type, cell_line, dens_mut, threshold_num, min_count))
  output_off <- file.path(base_dir, sprintf("%d_%s_%s_%s_Ranked_hotspots_%s_OFF_target_%d%d.tsv",
                                           enlargement, sample, type, cell_line, dens_mut, threshold_num, min_count))
  
  # Write outputs
  write_tsv(input_on_target, output_on)
  write_tsv(input_off_target, output_off)
  
  # Print summary
  cat("Output files saved:\n")
  cat(output_on, "\n")
  cat(output_off, "\n")
  cat("Processing completed successfully.\n")
}

# Run main with command line arguments
args <- commandArgs(trailingOnly = TRUE)
main(args)
