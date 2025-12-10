#!/usr/bin/env Rscript

# nohup Rscript 06_OnOff_DepIndep_finding.R CH12F3 3 5 0 UM_TKO_CIT_Duv_CLEAN --Dep --dens > nohup_350Dep_UM_TKO_CIT_Duv_density.out 2>&1 &
# Note: This script processes BOTH the original Ranked_hotspots and the _blackfilt.tsv files.

library(GenomicRanges)
library(readr)
library(dplyr) 

##########################
#       FUNCTIONS        #
##########################

# Function to log messages (Re-introduced for consistency)
log_message <- function(type, message) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("%s [%s] %s\n", timestamp, toupper(type), message))
}

# Parse command-line flags (Original logic maintained)
parse_flags <- function(flags) {

    RUN_DEP <- TRUE; RUN_INDEP <- TRUE; RUN_MUT <- TRUE; RUN_DENS <- TRUE
    TYPE_SET <- FALSE; MODE_SET <- FALSE; ALL_SEEN <- FALSE
    
    if (length(flags) > 0) {
        flags <- flags[flags != ""]
        for (f in flags) {
            if (f == "--ALL") { ALL_SEEN <- TRUE }
            if (f == "--Dep") { RUN_DEP <- TRUE; RUN_INDEP <- FALSE; TYPE_SET <- TRUE }
            if (f == "--Indep") { RUN_DEP <- FALSE; RUN_INDEP <- TRUE; TYPE_SET <- TRUE }
            if (f == "--mut") { RUN_MUT <- TRUE; RUN_DENS <- FALSE; MODE_SET <- TRUE }
            if (f == "--dens") { RUN_MUT <- FALSE; RUN_DENS <- TRUE; MODE_SET <- TRUE }
        }
        if (ALL_SEEN && !TYPE_SET) { RUN_DEP <- TRUE; RUN_INDEP <- TRUE }
        if (ALL_SEEN && !MODE_SET) { RUN_MUT <- TRUE; RUN_DENS <- TRUE }
    }
    
    log_message("INFO", "Flags parsed:")
    log_message("INFO", paste("RUN_DEP =", RUN_DEP))
    log_message("INFO", paste("RUN_INDEP =", RUN_INDEP))
    log_message("INFO", paste("RUN_MUT =", RUN_MUT))
    log_message("INFO", paste("RUN_DENS =", RUN_DENS))

    return(list(RUN_DEP=RUN_DEP, RUN_INDEP=RUN_INDEP, RUN_MUT=RUN_MUT, RUN_DENS=RUN_DENS))
}


# Function to check file existence with error message
check_file_exists <- function(file_path) {
 if (!file.exists(file_path)) {
  stop(sprintf("[ERROR] File '%s' does not exist!", file_path))
 } else {
  log_message("INFO", paste("File exists:", file_path))
 }
}

# Function to read BED-like file (chrom, start, end)
read_regions <- function(filepath) {
 check_file_exists(filepath)
 # Target regions do not have header
 regions <- readr::read_tsv(filepath, col_names = FALSE, col_types = readr::cols(), col_select = 1:3)
 colnames(regions) <- c("chrom", "start", "end")
 return(regions)
}

# Function to create GRanges from dataframe with chrom, start, end
create_granges <- function(df) {
 return(GenomicRanges::GRanges(seqnames = df$chrom,
  ranges = IRanges::IRanges(start = df$start, end = df$end)))
}

# Function to find overlapping rows between input and target regions
get_overlapping_subset <- function(input_df, target_granges) {
 input_gr <- create_granges(input_df)
 overlaps <- GenomicRanges::findOverlaps(input_gr, target_granges)
 hits <- S4Vectors::queryHits(overlaps)
 subset_df <- input_df[hits, ]
 return(subset_df)
}


# CORE LOGIC: Contains the main processing logic
process_hotspots <- function(cell_line, threshold_num, min_count, enlargement, sample, type, dens_mut, input_suffix) {
    
    log_message("INFO", sprintf("Starting processing: Type=%s, Mode=%s, Input Suffix=%s", type, dens_mut, input_suffix))

    # Determine the directory name based on 'type' (Dep/Indep)
    type_dir_name <- if (type == "Dep") {
        "Dependent_hotspots"
    } else if (type == "Indep") {
        "Independent_hotspots"
    } else {
        stop("[ERROR] Invalid type specified.")
    }

    # Build base_dir up to ranking folder
    base_dir_ranking <- file.path("/scratch/DepIndep_dataset",
                            paste0("Threshold_", threshold_num),
                            paste0("Enlargement_", enlargement),
                            type_dir_name,
                            sample,
                            paste0("Ranking_", type, "_hotspots"),
                            dens_mut)

    # Define the input path based on suffix
    if (input_suffix == "_blackfilt") {
        # Filtered Hotspot Case (found in a 'Filtered_hotspots' subdirectory)
        input_dir <- file.path(base_dir_ranking, "Filtered_hotspots")
        input_filename <- sprintf("%d_%s_%s_%s_Ranked_hotspots_%s_%d%d%s.tsv",
                                  enlargement, sample, type, cell_line, dens_mut, threshold_num, min_count, input_suffix)
    } else {
        # Original Hotspot Case (found in the ranking folder)
        input_dir <- base_dir_ranking
        input_filename <- sprintf("%d_%s_%s_%s_Ranked_hotspots_%s_%d%d.tsv",
                                  enlargement, sample, type, cell_line, dens_mut, threshold_num, min_count)
    }
    
    input_filepath <- file.path(input_dir, input_filename)
    
    # Target region files (Are the same for both inputs)
    target_dir <- "/scratch/reference_data"
    on_target_file <- file.path(target_dir, "ON_target_CH12F3.bed")
    off_target_file <- file.path(target_dir, "OFF_target_CH12F3.bed")
    
    # Check required input files
    check_file_exists(input_filepath)
    check_file_exists(on_target_file)
    check_file_exists(off_target_file)
    
    # Read input data (Assume ranked hotspots have header)
    input_data <- readr::read_tsv(input_filepath, col_names = TRUE, col_types = readr::cols())
    on_target_regions <- read_regions(on_target_file)
    off_target_regions <- read_regions(off_target_file)
    
    # Create GRanges for target regions
    gr_on_target <- create_granges(on_target_regions)
    gr_off_target <- create_granges(off_target_regions)
    
    # Find overlaps and subset input data
    log_message("INFO", "Finding ON-target overlaps...")
    input_on_target <- get_overlapping_subset(input_data, gr_on_target)
    log_message("INFO", "Finding OFF-target overlaps...")
    input_off_target <- get_overlapping_subset(input_data, gr_off_target)
    
    # Define output directory and filenames
    # The output goes into an 'OnOff_target_ranking' subdirectory created INSIDE the input directory
    output_dir <- file.path(input_dir, "OnOff_target_ranking") 
    
    if (!dir.exists(output_dir)) {
        log_message("INFO", paste("The output directory doesn't exist, creating it:", output_dir))
        dir.create(output_dir, recursive = TRUE)
    }
    
    # Output filenames include the suffix to distinguish between original and filtered runs
    output_on <- file.path(output_dir, sprintf("%d_%s_%s_%s_Ranked_hotspots_%s_ON_target_%d%d%s.tsv",
                                               enlargement, sample, type, cell_line, dens_mut, threshold_num, min_count, input_suffix))

    output_off <- file.path(output_dir, sprintf("%d_%s_%s_%s_Ranked_hotspots_%s_OFF_target_%d%d%s.tsv",
                                                enlargement, sample, type, cell_line, dens_mut, threshold_num, min_count, input_suffix))
    
    # Write outputs
    readr::write_tsv(input_on_target, output_on)
    readr::write_tsv(input_off_target, output_off)
    
    # Print summary
    log_message("INFO", "\nOutput files saved:")
    log_message("INFO", paste("ON-Target:", output_on))
    log_message("INFO", paste("OFF-Target:", output_off))
    log_message("INFO", sprintf("\nProcessing of %s completed successfully.", input_filepath))
}


##########################
#          MAIN          #
##########################
# Main script function 
main <- function(args) {
    # Expected Positional Arguments: cell_line, threshold_num, min_count, enlargement, sample
    # Flags follow: --Dep/--Indep, --dens/--mut
    if (length(args) < 5) {
        stop("[ERROR] At least 5 positional parameters required: cell_line, threshold_num, min_count, enlargement, sample. Flags follow.")
    }
    
    # Parse positional arguments
    cell_line <- args[1]
    threshold_num <- as.integer(args[2])
    min_count <- as.integer(args[3])
    enlargement <- as.integer(args[4])
    sample <- args[5]

    # Flags start from args[6]
    extra_flags <- if (length(args) >= 6) {
        args[6:length(args)]
    } else {
        character(0)
    }
    flags <- parse_flags(extra_flags)

    # Build combinations based on flags (original logic maintained)
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

    # Loop over each combination and process TWICE (Original and Filtered)
    for (combo in combinations) {
        type_i <- combo$type
        dens_mut_i <- combo$dens_mut
        
        # 1. Processing of Original Hotspots (input_suffix = "")
        process_hotspots(cell_line, threshold_num, min_count, enlargement, sample, type_i, dens_mut_i, input_suffix="")

        # 2. Processing of Blacklist Filtered Hotspots (input_suffix = "_blackfilt")
        process_hotspots(cell_line, threshold_num, min_count, enlargement, sample, type_i, dens_mut_i, input_suffix="_blackfilt")
    }

    log_message("INFO", "All requested combinations processed successfully.")
}

# Run main with command line arguments
args <- commandArgs(trailingOnly = TRUE)
main(args)