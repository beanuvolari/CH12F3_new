#!/usr/bin/env Rscript

# nohup Rscript 33_7Detect-seqinHTGTS_hotspots.R 3 5 0 Dep UM_TKO_CIT_Duv_CLEAN nmutations CH12F3 > nohup_350Dep_UMTKODuv_Detectnmut_blackfilt.out 2>&1 & 
# nohup Rscript 33_7Detect-seqinHTGTS_hotspots.R 3 5 0 Indep UM_TKO_CIT_Duv_CLEAN nmutations CH12F3  > nohup_350Indep_UMTKODuv_Detectnmut.out 2>&1 &

# ==========================
# Detect-seq vs HTGTS script
# ==========================
# Question: How many Detect-seq hotspots fall in HTGTS hotspots (ON and OFF targets)?

library(GenomicRanges)
library(readr)

########################
#      Functions       #
########################

check_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(paste("[ERROR]", label, "file", path, "does not exist!"))
  } else {
    cat("[INFO]", label, "file found:", path, "\n")
  }
}

make_GRanges <- function(df) {
  GRanges(seqnames = df$chrom, ranges = IRanges(start = df$start, end = df$end))
}

# -------- Core function --------

process_detect_file <- function(input_path, output_dir, output_prefix, 
                                gr_ONtarget, gr_OFFtarget, summary_file) {
  
  check_file(input_path, "Input")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("[INFO] Output directory created:", output_dir, "\n")
  }
  
  # Load Detect-seq hotspots
  input_data <- read_tsv(input_path, col_names = TRUE)
  gr_input   <- make_GRanges(input_data)
  
  # Overlaps
  overlaps_ON  <- findOverlaps(gr_input, gr_ONtarget)
  overlaps_OFF <- findOverlaps(gr_input, gr_OFFtarget)
  
  # Subsets
  input_ON_target  <- input_data[subjectHits(overlaps_ON), ]
  input_OFF_target <- input_data[subjectHits(overlaps_OFF), ]
  
  # Save ON/OFF hits
  write_tsv(input_ON_target,  file.path(output_dir, paste0(output_prefix, "_ON.tsv")))
  write_tsv(input_OFF_target, file.path(output_dir, paste0(output_prefix, "_OFF.tsv")))
  
  # Coverage stats
  total       <- length(gr_input)
  covered_on  <- length(unique(queryHits(overlaps_ON)))
  covered_off <- length(unique(queryHits(overlaps_OFF)))
  not_known   <- total - (covered_on + covered_off)
  
  cat("[INFO] Coverage stats:\n",
      "  ON-target:  ", covered_on, "out of", total, "\n",
      "  OFF-target: ", covered_off, "out of", total, "\n",
      "  Not-known:  ", not_known, "out of", total, "\n")
  
  # Save summary
  pie_data <- data.frame(
    Target = c("ON_target_covered", "OFF_target_covered", "Not_known_targets"),
    Count  = c(covered_on, covered_off, not_known)
  )
  write_tsv(pie_data, file.path(output_dir, summary_file))
  
  cat("[INFO] Summary saved:", file.path(output_dir, summary_file), "\n\n")
}
process_detect_file_all_input <- function(input_path, output_dir, output_prefix, 
                                         gr_ONtarget, gr_OFFtarget, summary_file) {
  
  check_file(input_path, "Input")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("[INFO] Output directory created:", output_dir, "\n")
  }
  
  input_data <- read_tsv(input_path, col_names = TRUE)
  gr_input   <- make_GRanges(input_data)
  
  overlaps_ON  <- findOverlaps(gr_input, gr_ONtarget)
  overlaps_OFF <- findOverlaps(gr_input, gr_OFFtarget)
  
  # Subset input rows overlapped (qui NO unique, perché conto tutte le righe)
  input_ON_target  <- input_data[queryHits(overlaps_ON), ]
  input_OFF_target <- input_data[queryHits(overlaps_OFF), ]
  
  write_tsv(input_ON_target,  file.path(output_dir, paste0(output_prefix, "_ON.tsv")))
  write_tsv(input_OFF_target, file.path(output_dir, paste0(output_prefix, "_OFF.tsv")))
  
  # Coverage stats: qui conti tutte le righe input sovrapposte (anche duplicate)
  covered_on  <- length(queryHits(overlaps_ON))
  covered_off <- length(queryHits(overlaps_OFF))
  total       <- length(gr_input)
  not_known   <- total - (covered_on + covered_off)
  
  cat("[INFO] Coverage stats (all overlapping input rows):\n",
      "  ON-target:  ", covered_on, "out of", total, "\n",
      "  OFF-target: ", covered_off, "out of", total, "\n",
      "  Not-known:  ", not_known, "out of", total, "\n")
  
  pie_data <- data.frame(
    Target = c("ON_target_covered", "OFF_target_covered", "Not_known_targets"),
    Count  = c(covered_on, covered_off, not_known)
  )
  write_tsv(pie_data, file.path(output_dir, summary_file))
  
  cat("[INFO] Summary saved:", file.path(output_dir, summary_file), "\n\n")
}


########################
#         Main         #
########################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop("Usage: Rscript script.R threshold_num min_count enlargement type sample dens_mut cell_line")
}

threshold_num <- as.integer(args[1])
min_count     <- as.integer(args[2])
enlargement   <- as.integer(args[3])
type          <- args[4]
sample        <- args[5]
dens_mut      <- args[6]
cell_line     <- args[7]

# Input files
file_dir <- file.path("Dep_Indep_dataset", 
                      paste0("Threshold_", threshold_num),
                      paste0("Enlargement_", enlargement), sample, 
                      paste0("Ranking_", type, "_hotspots"), dens_mut)

detect_full_file <- file.path(file_dir, sprintf("%s_%s_%s_%s_Ranked_hotspots_%s_%d%d.tsv",
                                                enlargement, sample, type, cell_line, dens_mut, threshold_num, min_count))
detect_top200_file <- file.path(file_dir, sprintf("%s_%s_%s_%s_Ranked200hotspots_%s_%d%d.tsv",
                                                  enlargement, sample, type, cell_line, dens_mut, threshold_num, min_count))


detect_full_blackfilt <- file.path(file_dir, sprintf("%s_%s_%s_%s_Ranked_hotspots_%s_%d%d_blackfilt.tsv",
                                                enlargement, sample, type, cell_line, dens_mut, threshold_num, min_count))
detect_top200_blackfilt <- file.path(file_dir, sprintf("%s_%s_%s_%s_Ranked200hotspots_%s_%d%d_blackfilt.tsv",
                                                  enlargement, sample, type, cell_line, dens_mut, threshold_num, min_count))

check_file(detect_full_file, "Detect-seq full hotspots")
check_file(detect_top200_file, "Detect-seq top200 hotspots")

check_file(detect_full_blackfilt, "Detect-seq full hotspots filtered")
check_file(detect_top200_blackfilt, "Detect-seq top200 hotspots filtered")

# Reference files
On_target_path  <- file.path("Controls", "ON_target_CH12F3.bed")
Off_target_path <- file.path("Controls", "OFF_target_CH12F3.bed")
check_file(On_target_path, "ON-target")
check_file(Off_target_path, "OFF-target")

# Load references
On_target  <- read_tsv(On_target_path,  col_names = FALSE, col_types = cols(), col_select = 1:3)
Off_target <- read_tsv(Off_target_path, col_names = FALSE, col_types = cols(), col_select = 1:3)
colnames(On_target) <- colnames(Off_target) <- c("chrom", "start", "end")
gr_ONtarget  <- make_GRanges(On_target)
gr_OFFtarget <- make_GRanges(Off_target)

# Output dir
output_dir_full <- file.path(file_dir, "Detect-seqvsHTGTS", "FULL")
output_dir_200 <- file.path(file_dir, "Detect-seqvsHTGTS", "first_200")

# Run analyses
process_detect_file_all_input(detect_full_file, output_dir_full, 
                    sprintf("%s_%s_%s_Rankedhot_%s_DetecttargetFULL_%d%d", 
                            enlargement, sample, type, dens_mut, threshold_num, min_count),
                    gr_ONtarget, gr_OFFtarget,
                    sprintf("%s_%s_%s_%s_TargetCoverageSummaryFULL_Detect_%d%d.tsv", 
                            enlargement, sample, type, dens_mut, threshold_num, min_count))

process_detect_file_all_input(detect_top200_file, output_dir_200, 
                    sprintf("%s_%s_%s_Rankedhot200_%s_Detecttarget_%d%d", 
                            enlargement, sample, type, dens_mut, threshold_num, min_count),
                    gr_ONtarget, gr_OFFtarget,
                    sprintf("%s_%s_%s_%s_TargetCoverageSummary200_Detect_%d%d.tsv", 
                            enlargement, sample, type, dens_mut, threshold_num, min_count))


process_detect_file_all_input(detect_full_blackfilt, output_dir_full, 
                    sprintf("%s_%s_%s_Rankedhot_%s_DetecttargetFULL_%d%d_blackfilt", 
                            enlargement, sample, type, dens_mut, threshold_num, min_count),
                    gr_ONtarget, gr_OFFtarget,
                    sprintf("%s_%s_%s_%s_TargetCoverageSummaryFULL_Detect_%d%d_blackfilt.tsv", 
                            enlargement, sample, type, dens_mut, threshold_num, min_count))

process_detect_file_all_input(detect_top200_blackfilt, output_dir_200, 
                    sprintf("%s_%s_%s_Rankedhot200_%s_Detecttarget_%d%d_blackfilt", 
                            enlargement, sample, type, dens_mut, threshold_num, min_count),
                    gr_ONtarget, gr_OFFtarget,
                    sprintf("%s_%s_%s_%s_TargetCoverageSummary200_Detect_%d%d_blackfilt.tsv", 
                            enlargement, sample, type, dens_mut, threshold_num, min_count))




cat("[INFO] All processing completed successfully!\n")
