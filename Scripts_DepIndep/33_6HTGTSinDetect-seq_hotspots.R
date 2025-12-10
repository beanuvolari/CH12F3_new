#!/usr/bin/env Rscript

# nohup Rscript 33_6HTGTSinDetect-seq_hotspots.R 3 5 0 Dep UM_TKO_CIT_Duv_CLEAN nmutations CH12F3 > nohup_350Dep_UMTKODuv_HTGTSnmut_blackfilt.out 2>&1 & 
# nohup Rscript 33_6HTGTSinDetect-seq_hotspots.R 3 5 0 Indep UM_TKO_CIT_Duv_CLEAN nmutations CH12F3 > nohup_350Indep_UMTKODuv_HTGTSnmut.out 2>&1 &

# ==========================
# HTGTS vs Detect-seq script
# ==========================
# How many HTGTS hotspots are covered by dtect-seq hotspots?
# Thi script measures te coverage of HTGTS hotspots by Detect-seq hotspots

library(GenomicRanges)
library(readr)

########################
#      Functions       #
########################

# Check if file exists
check_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(paste("[ERROR]", label, "file", path, "doesn't exist!"))
  } else {
    cat("[INFO]", label, "file found:", path, "\n")
  }
}

# Convert data.frame with chrom/start/end to GRanges
make_GRanges <- function(df) {
  GRanges(seqnames = df$chrom, ranges = IRanges(start = df$start, end = df$end))
}

# Save regions not covered
save_not_covered <- function(gr_target, overlaps, output_path) {
  not_cov_idx <- setdiff(seq_along(gr_target), unique(queryHits(overlaps)))
  not_cov_df  <- as.data.frame(gr_target[not_cov_idx])
  write_tsv(not_cov_df, output_path)
  cat("[INFO] Not-covered regions saved in:", output_path, "\n")
}

# Compute and save coverage summary
summarize_coverage <- function(total_on, covered_on, total_off, covered_off, file1, file2) {
  
  not_covered_on  <- total_on - covered_on
  not_covered_off <- total_off - covered_off
  
  # Detailed (ON + OFF split)
  pie_data1 <- data.frame(
    Target = c("ON_target_covered", "ON_target_not_covered", 
               "OFF_target_covered", "OFF_target_not_covered"),
    Count = c(covered_on, not_covered_on, covered_off, not_covered_off)
  )
  write_tsv(pie_data1, file1)
  cat("[INFO] Coverage summary saved in:", file1, "\n")
  
  # Global (ON+OFF vs Not covered)
  total <- total_on + total_off
  not_covered <- total - (covered_on + covered_off)
  
  pie_data2 <- data.frame(
    Target = c("ON_target_covered", "OFF_target_covered", "Not_covered_targets"),
    Count = c(covered_on, covered_off, not_covered)
  )
  write_tsv(pie_data2, file2)
  cat("[INFO] Coverage summary saved in:", file2, "\n")
}

# -------- Core function --------

process_input_file_HTGTS <- function(input_path, output_dir, output_prefix, 
                                     gr_ONtarget, gr_OFFtarget, summary_file, summary_2file) {
  
  check_file(input_path, "Input")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("[INFO] Output directory created:", output_dir, "\n")
  }
  
  # Read input
  input_data <- read_tsv(input_path, col_names = TRUE)
  gr_input   <- make_GRanges(input_data)
  
  # Overlaps
  overlaps_ON  <- findOverlaps(gr_ONtarget, gr_input)
  overlaps_OFF <- findOverlaps(gr_OFFtarget, gr_input)
  
  # Subset input
  on_hits_input  <- unique(subjectHits(overlaps_ON))
  off_hits_input <- unique(subjectHits(overlaps_OFF))
  
  input_ON_target  <- input_data[on_hits_input, ]
  input_OFF_target <- input_data[off_hits_input, ]
  
  # Save ON/OFF hits
  write_tsv(input_ON_target,  file.path(output_dir, paste0(output_prefix, "_ON.tsv")))
  write_tsv(input_OFF_target, file.path(output_dir, paste0(output_prefix, "_OFF.tsv")))
  
  # Coverage stats
  total_on  <- length(gr_ONtarget)
  covered_on <- length(unique(queryHits(overlaps_ON)))
  total_off  <- length(gr_OFFtarget)
  covered_off <- length(unique(queryHits(overlaps_OFF)))
  
  cat("[INFO] Coverage:\n",
      "ON-target: ", covered_on, "of", total_on, "\n",
      "OFF-target:", covered_off, "of", total_off, "\n")
  
  # Save not-covered regions
  save_not_covered(gr_ONtarget, overlaps_ON, 
                   file.path(output_dir, paste0(output_prefix, "_ON_notcovered.tsv")))
  save_not_covered(gr_OFFtarget, overlaps_OFF, 
                   file.path(output_dir, paste0(output_prefix, "_OFF_notcovered.tsv")))
  
  # Save summaries
  summarize_coverage(total_on, covered_on, total_off, covered_off, 
                     file.path(output_dir, summary_file),
                     file.path(output_dir, summary_2file))
  
  cat("[INFO] Processing finished for", input_path, "\n\n")
}

########################
#         Main         #
########################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop("[INFO] Usage: Rscript script.R threshold_num min_count enlargement type sample dens_mut cell_line")
}

threshold_num <- as.integer(args[1])
min_count     <- as.integer(args[2])
enlargement   <- as.integer(args[3])
type          <- args[4]
sample        <- args[5]
dens_mut      <- args[6]
cell_line     <- args[7]

# Input paths
file_dir <- file.path("Dep_Indep_dataset", 
                      paste0("Threshold_", threshold_num),
                      paste0("Enlargement_", enlargement), sample, 
                      paste0("Ranking_", type, "_hotspots"), dens_mut)

file_name     <- sprintf("%s_%s_%s_%s_Ranked_hotspots_%s_%d%d_blackfilt.tsv",
                         enlargement, sample, type, cell_line, dens_mut, threshold_num, min_count)
file200_name  <- sprintf("%s_%s_%s_%s_Ranked200hotspots_%s_%d%d_blackfilt.tsv",
                         enlargement, sample, type, cell_line, dens_mut, threshold_num, min_count)

input_file    <- file.path(file_dir, file_name)
input_file200 <- file.path(file_dir, file200_name)
check_file(input_file, "Input full hotspots")
check_file(input_file200, "Input top200 hotspots")

# Check ON/OFF references
On_target_path  <- file.path("Controls", "ON_target_CH12F3.bed")
Off_target_path <- file.path("Controls", "OFF_target_CH12F3.bed")
check_file(On_target_path, "ON-target")
check_file(Off_target_path, "OFF-target")

# Read ON/OFF target regions
On_target  <- read_tsv(On_target_path,  col_names = FALSE, col_types = cols(), col_select = 1:3)
Off_target <- read_tsv(Off_target_path, col_names = FALSE, col_types = cols(), col_select = 1:3)
colnames(On_target)  <- colnames(Off_target) <- c("chrom", "start", "end")

gr_ONtarget  <- make_GRanges(On_target)
gr_OFFtarget <- make_GRanges(Off_target)

# Output directory
output_dir_full  <- file.path(file_dir, "HTGTSvsDetect-seq", "FULL")
output_dir_200   <- file.path(file_dir, "HTGTSvsDetect-seq", "first_200")

# Output prefixes
output_prefix    <- sprintf("%s_%s_%s_%s_Rankedhot_%s_HTGTStargetFULL_%d%d_blackfilt", enlargement, sample, type, cell_line, dens_mut, threshold_num, min_count)
output_prefix200 <- sprintf("%s_%s_%s_%s_Rankedhot200_%s_HTGTStarget_%d%d_blackfilt", enlargement, sample, type, cell_line, dens_mut, threshold_num, min_count)

# Run analyses
process_input_file_HTGTS(input_file,    output_dir_full, output_prefix,    gr_ONtarget, gr_OFFtarget, 
                         sprintf("%s_%s_%s_%s_TargetCoverageSummaryFULL_HTGTS_%d%d_blackfilt.tsv", enlargement, sample, type, dens_mut, threshold_num, min_count),
                         sprintf("%s_%s_%s_%s_2TargetCoverageSummaryFULL_HTGTS_%d%d_blackfilt.tsv", enlargement, sample, type, dens_mut, threshold_num, min_count))

process_input_file_HTGTS(input_file200, output_dir_200, output_prefix200, gr_ONtarget, gr_OFFtarget, 
                         sprintf("%s_%s_%s_%s_TargetCoverageSummary200_HTGTS_%d%d_blackfilt.tsv", enlargement, sample, type, dens_mut, threshold_num, min_count),
                         sprintf("%s_%s_%s_%s_2TargetCoverageSummary200_HTGTS_%d%d_blackfilt.tsv", enlargement, sample, type, dens_mut, threshold_num, min_count))

cat("[INFO] All processing completed successfully!\n")
