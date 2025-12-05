#!/usr/bin/env Rscript

# nohup Rscript 05_OnOff_DepIndep_finding.R 3 5 0 CH12F3 UM_TKO_CIT_Duv_CLEAN density Dep > nohup_350Dep_UM_TKO_CIT_Duv_density.out 2>&1 &
# nohup Rscript 05_OnOff_DepIndep_finding.R 3 5 0 CH12F3 UM_TKO_CIT_Duv_CLEAN nmutations Dep > nohup_350Dep_UM_TKO_CIT_Duv_nmutations.out 2>&1 &

# nohup Rscript 05_OnOff_DepIndep_finding.R 3 5 0 CH12F3 UM_TKO_CIT_Duv_CLEAN density Indep > nohup_350Indep_UM_TKO_CIT_Duv_density.out 2>&1 &
# nohup Rscript 05_OnOff_DepIndep_finding.R 3 5 0 CH12F3 UM_TKO_CIT_Duv_CLEAN nmutations Indep > nohup_350Indep_UM_TKO_CIT_Duv_nmutations.out 2>&1 &

library(GenomicRanges)
library(readr)

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


# Function to check file existence with error message
check_file_exists <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(sprintf("[ERROR] file '%s' does not exist!", file_path))
  } else {
    cat(sprintf("[INFO] File exists: %s\n", file_path))
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
  if (length(args) < 5) {
    stop("[ERROR] At least 5 parameters required: cell_line, threshold_num, min_count, enlargement, sample, type (Dep or Indep), dens_mut (density or nmutations).")
  }
  
  # Parse arguments
  cell_line <- args[1]
  threshold_num <- as.integer(args[2])
  min_count <- as.integer(args[3])
  enlargement <- as.integer(args[4])
  sample <- args[5]


  # Extra arguments are flags
  extra_flags <- args[8:length(args)]
  flags <- parse_flags(extra_flags)

  # Build combinations to run based on flags
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
  

  # Loop over each combination and process
  for (combo in combinations) {
    type_i <- combo$type
    dens_mut_i <- combo$dens_mut
    
    cat(sprintf("[INFO] Processing: type=%s, dens_mut=%s\n", type_i, dens_mut_i))
    # Build file paths
    base_dir <- file.path("/scratch/DepIndep_dataset",
                      paste0("Threshold_", threshold_num),
                      paste0("Enlargement_", enlargement),
                      paste0(type_i, "endent_hotspots"),
                      sample,
                      paste0("Ranking_", type_i, "_hotspots"),
                      dens_mut_i)

    input_filename <- sprintf("%d_%s_%s_%s_Ranked_hotspots_%s_%d%d.tsv",
                          enlargement, sample, type_i, cell_line, dens_mut_i, threshold_num, min_count)
    input_filepath <- file.path(base_dir, input_filename)
    
    # Target region files
    target_dir <- "/scratch/reference_data"
    on_target_file <- file.path(target_dir, "ON_target_CH12F3.bed")
    off_target_file <- file.path(target_dir, "OFF_target_CH12F3.bed")
    
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
    output_dir <- file.path(base_dir, "OnOff_target_ranking")
    if (!dir.exists(output_dir)) {
      cat("[INFO] The output directory doesn't exist, creating it...\n")
      dir.create(output_dir, recursive = TRUE)
    }
    output_on <- file.path(output_dir, sprintf("%d_%s_%s_%s_Ranked_hotspots_%s_ON_target_%d%d.tsv",
                                           enlargement, sample, type_i, cell_line, dens_mut_i, threshold_num, min_count))

    output_off <- file.path(output_dir, sprintf("%d_%s_%s_%s_Ranked_hotspots_%s_OFF_target_%d%d.tsv",
                                            enlargement, sample, type_i, cell_line, dens_mut_i, threshold_num, min_count))
    
    # Write outputs
    write_tsv(input_on_target, output_on)
    write_tsv(input_off_target, output_off)
    
    # Print summary
    cat("[INFO] Output files saved:\n")
    cat(output_on, "\n")
    cat(output_off, "\n")
    cat("[INFO] Processing completed successfully.\n")
  }
  cat("[INFO] All requested combinations processed successfully.\n")
}

# Run main with command line arguments
args <- commandArgs(trailingOnly = TRUE)
main(args)
