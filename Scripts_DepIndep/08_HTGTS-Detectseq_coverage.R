#!/usr/bin/env Rscript
# nohup Rscript 08_HTGTS-Detectseq_coverage.R CH12F3 3 5 0 UM_TKO_CIT_Duv --Dep --mut > nohup_08HTGTS-Detectseq_coverage.out 2>&1 &
# =========================================================================
# 08_HTGTS_Detectseq_Coverage_Unified.R
# This script unifies the overlap analysis between Detect-seq hotspots 
# and HTGTS target regions (ON/OFF). It performs analysis on all 
# combinations (Full, Top 200, Filtered) based on flags.
# =========================================================================

# Libraries required for genomic intervals and data handling
library(GenomicRanges)
library(readr)
library(dplyr)

########################
#       Functions      #
########################

# [INFO] Function to check if a file exists
check_file <- function(path, label) {
  if (!file.exists(path)) {
    cat("[WARNING]", label, "file", path, "doesn't exist! Skipping related analysis.\n")
    return(FALSE)
  } else {
    cat("[INFO]", label, "file found:", path, "\n")
    return(TRUE)
  }
}

# [INFO] Function to convert a data frame with chrom/start/end to GRanges
make_GRanges <- function(df) {
  GRanges(seqnames = df$chrom, ranges = IRanges(start = df$start, end = df$end))
}

# [INFO] Core Unified Analysis Function with dynamic output naming
analyze_coverage <- function(gr_detect, gr_ONtarget, gr_OFFtarget, outputDetect_dir, outputHTGTS_dir,
                             enlargement, sample, type, cell_line, mode, 
                             threshold, mincount, is_200, is_filtered) {
  
  # Construct dynamic base prefix: enlargement_sample_type_cellline_mode
  base_prefix <- paste(enlargement, sample, type, cell_line, mode, sep = "_")
  
  # Define tags for naming
  rank_tag <- if (is_200) "Summary200" else "Summary"
  filt_tag <- if (is_filtered) "_blackfilt" else ""
  
  cat(sprintf("[INFO] Starting analysis for: %s (%s)\n", base_prefix, rank_tag))
  
  # ----------------------------------------------------
  # A) Detect-seq vs HTGTS (How many Detect hotspots hit targets?)
  # ----------------------------------------------------
  overlaps_ON <- findOverlaps(gr_detect, gr_ONtarget)
  overlaps_OFF <- findOverlaps(gr_detect, gr_OFFtarget)
  
  total_detect <- length(gr_detect)
  covered_on <- length(unique(queryHits(overlaps_ON)))
  covered_off <- length(unique(queryHits(overlaps_OFF)))
  not_known <- total_detect - length(unique(c(queryHits(overlaps_ON), queryHits(overlaps_OFF))))
  
  summary_detect <- data.frame(
    Target = c("ON_targets_covered", "OFF_targets_covered", "Not_known_targets"),
    Count = c(covered_on, covered_off, not_known)
  )
  
  # Format: enlargement_sample_type_cellline_mode_TargetCoverageSummary(200)_Detect_thresholdmincount(_blackfilt).tsv
  out_name_detect <- paste0(base_prefix, "_TargetCoverage", rank_tag, "_Detect_", 
                            threshold, mincount, filt_tag, ".tsv")
  write_tsv(summary_detect, file.path(outputDetect_dir, out_name_detect))

  cat(sprintf("[INFO] Analysis finished. 1 table generated in: %s\n", outputDetect_dir))

  # ----------------------------------------------------
  # B) HTGTS vs Detect-seq (Sensitivity: 2 Outputs)
  # ----------------------------------------------------
  overlaps_ON_rev <- findOverlaps(gr_ONtarget, gr_detect)
  overlaps_OFF_rev <- findOverlaps(gr_OFFtarget, gr_detect)

  # Counts for ON targets
  total_on <- length(gr_ONtarget)
  covered_on_h <- length(unique(queryHits(overlaps_ON_rev)))
  not_covered_on <- total_on - covered_on_h

  # Counts for OFF targets
  total_off <- length(gr_OFFtarget)
  covered_off_h <- length(unique(queryHits(overlaps_OFF_rev)))
  not_covered_off <- total_off - covered_off_h

  # B.1) Detailed Summary
  summary_htgts_detailed <- data.frame(
    Target = c("ON_targets_covered", "ON_targets_not_covered", 
               "OFF_targets_covered", "OFF_targets_not_covered"),
    Count = c(covered_on_h, not_covered_on, covered_off_h, not_covered_off)
  )
  
  out_name_htgts_detailed <- paste0(base_prefix, "_TargetCoverage", rank_tag, "Detailed_HTGTS_", 
                                    threshold, mincount, filt_tag, ".tsv")
 write_tsv(summary_htgts_detailed, file.path(outputHTGTS_dir, out_name_htgts_detailed))

  # B.2) Global Summary
  not_covered_tot <- not_covered_on + not_covered_off
  summary_htgts_global <- data.frame(
    Target = c("ON_targets_covered", "OFF_targets_covered", "Targets_not_covered"),
    Count = c(covered_on_h, covered_off_h, not_covered_tot)
  )
  
  out_name_htgts_global <- paste0(base_prefix, "_TargetCoverage", rank_tag, "_HTGTS_", 
                                 threshold, mincount, filt_tag, ".tsv")
 write_tsv(summary_htgts_global, file.path(outputHTGTS_dir, out_name_htgts_global))
  
  cat(sprintf("[INFO] Analysis finished. 2 tables generated in: %s\n", outputHTGTS_dir))
}

########################
#          Main        #
########################

args <- commandArgs(trailingOnly = TRUE)

# [ERROR] Check for minimum positional arguments
if (length(args) < 5) {
  stop("[ERROR] At least 5 parameters required: cell_line, threshold_num, min_count, enlargement, sample.")
}

# Assign positional arguments
cell_line <- args[1]
threshold_num <- as.integer(args[2])
min_count <- as.integer(args[3])
enlargement <- as.integer(args[4])
sample <- args[5]

# [INFO] Parse optional flags for Dep/Indep and mode

parse_flags <- function(flags) {
  RUN_DEP <- FALSE; RUN_INDEP <- FALSE; RUN_MUT <- FALSE; RUN_DENS <- FALSE
  TYPE_SET <- FALSE; MODE_SET <- FALSE
  
  if (length(flags) > 0) {
    for (f in flags[flags != ""]) {
      if (f == "--Dep") { RUN_DEP <- TRUE; TYPE_SET <- TRUE }
      if (f == "--Indep") { RUN_INDEP <- TRUE; TYPE_SET <- TRUE }
      if (f == "--mut") { RUN_MUT <- TRUE; MODE_SET <- TRUE }
      if (f == "--dens") { RUN_DENS <- TRUE; MODE_SET <- TRUE }
    }
  }
  # Default behavior: run all if not specified
  if (!TYPE_SET) { RUN_DEP <- TRUE; RUN_INDEP <- TRUE }
  if (!MODE_SET) { RUN_MUT <- TRUE; RUN_DENS <- TRUE }
  cat("[INFO] Flags parsed:\n")
  cat("       RUN_DEP =", RUN_DEP, "\n")
  cat("       RUN_INDEP =", RUN_INDEP, "\n")
  cat("       RUN_MUT =", RUN_MUT, "\n")
  cat("       RUN_DENS =", RUN_DENS, "\n")

  return(list(RUN_DEP=RUN_DEP, RUN_INDEP=RUN_INDEP, RUN_MUT=RUN_MUT, RUN_DENS=RUN_DENS))
}

flags <- parse_flags(if (length(args) >= 6) args[6:length(args)] else character(0))

# Build combinations for the loop
combinations <- list()
if (flags$RUN_DEP) combinations <- append(combinations, list(list(type = "Dep", folder = "Dependent_hotspots")))
if (flags$RUN_INDEP) combinations <- append(combinations, list(list(type = "Indep", folder = "Independent_hotspots")))
modes <- c(); if (flags$RUN_DENS) modes <- c(modes, "density"); if (flags$RUN_MUT) modes <- c(modes, "nmutations")

# [INFO] Load HTGTS Reference Files
reference_dir <- "scratch/reference_data"
On_path <- file.path(reference_dir, paste0("ON_target_", cell_line, ".bed"))
Off_path <- file.path(reference_dir, paste0("OFF_target_", cell_line, ".bed"))

if (!check_file(On_path, "ON-ref") || !check_file(Off_path, "OFF-ref")) {
  stop("[ERROR] Critical HTGTS reference files missing.")
}

gr_ONtarget <- read_tsv(On_path, col_names = FALSE, col_types = cols(), col_select = 1:3) %>% 
  rename(chrom=1, start=2, end=3) %>% make_GRanges()
gr_OFFtarget <- read_tsv(Off_path, col_names = FALSE, col_types = cols(), col_select = 1:3) %>% 
  rename(chrom=1, start=2, end=3) %>% make_GRanges()

# Track execution summary
files_processed <- c(); files_skipped <- c()

# [INFO] MAIN LOOP: Iterate over all specified combinations
for (comb in combinations) {
  for (m in modes) {
    cat(sprintf("\n[INFO] =====================================================\n"))
    cat(sprintf("[INFO] STARTING COMBINATION: %s | %s\n", comb$type, m))
    
    # Define root path for current combination
    path <- "/scratch/DepIndep_dataset"
    path_root <- file.path(path, paste0("Threshold_", threshold_num), 
                           paste0("Enlargement_", enlargement), comb$folder, sample, 
                           paste0("Ranking_", comb$type, "_hotspots"), m)
    
    # Process both Full_hotspots and Filtered_hotspots
    set_configs <- list(
      list(subfolder = "Full_hotspots", filtered = FALSE),
      list(subfolder = "Filtered_hotspots", filtered = TRUE)
    )
    
    for (set in set_configs) {
      # Process both standard and Top 200 files
      file_variants <- list(
        list(tag = "Ranked_hotspots", is_200 = FALSE),
        list(tag = "Ranked200hotspots", is_200 = TRUE)
      )
      
      for (variant in file_variants) {
        suffix <- if (set$filtered) "_blackfilt.tsv" else ".tsv"
        
        # Build input file name based on requirements
        input_name <- paste0(enlargement, "_", sample, "_", comb$type, "_", cell_line, "_", 
                             variant$tag, "_", m, "_", threshold_num, min_count, suffix)
        
        input_path <- file.path(path_root, set$subfolder, input_name)
        outputDetect_dir <- file.path(path_root, set$subfolder, "Pie_chart", "Detect-seqvsHTGTS")
        outputHTGTS_dir <- file.path(path_root, set$subfolder, "Pie_chart", "HTGTSvsDetect-seq")
        
        if (check_file(input_path, variant$tag)) {
          input_data <- read_tsv(input_path, col_names = TRUE, show_col_types = FALSE)
          if (nrow(input_data) > 0) {
            if (!dir.exists(outputDetect_dir)) dir.create(outputDetect_dir, recursive = TRUE)
            if (!dir.exists(outputHTGTS_dir)) dir.create(outputHTGTS_dir, recursive = TRUE)
            
            # Execute analysis
            analyze_coverage(gr_detect = make_GRanges(input_data), 
                             gr_ONtarget = gr_ONtarget, gr_OFFtarget = gr_OFFtarget, 
                             outputDetect_dir = outputDetect_dir, outputHTGTS_dir = outputHTGTS_dir, 
                             enlargement = enlargement,
                             sample = sample, type = comb$type, cell_line = cell_line, 
                             mode = m, threshold = threshold_num, mincount = min_count, 
                             is_200 = variant$is_200, is_filtered = set$filtered)
            
            files_processed <- c(files_processed, input_name)
          } else {
            cat(sprintf("[WARNING] File %s is empty. Skipping.\n", input_name))
            files_skipped <- c(files_skipped, paste(input_name, "(EMPTY)"))
          }
        } else {
          files_skipped <- c(files_skipped, paste(input_name, "(NOT FOUND)"))
        }
      }
    }
  }
}

# [INFO] Final Summary Report
cat("\n=====================================================\n")
cat("                ANALYSIS SUMMARY REPORT              \n")
cat("=====================================================\n")
cat(sprintf("[INFO] Successfully processed: %d\n", length(files_processed)))
if (length(files_skipped) > 0) {
  cat("[DETAILS] Skipped files:\n")
  for (f in files_skipped) cat(paste0("  - ", f, "\n"))
}
cat("=====================================================\n")
cat("[INFO] All processing tasks completed successfully!\n")