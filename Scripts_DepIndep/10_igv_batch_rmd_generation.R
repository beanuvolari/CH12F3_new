#!/usr/bin/env Rscript

################################################################################
# Script: 10_igv_batch_rmd_generation.R
# Description: Generates IGV resources into task-specific subdirectories
################################################################################

# [INFO] Setting up library paths
local_lib <- "/30tb/home/nuvobea/R/x86_64-pc-linux-gnu-library/4.4"
.libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(readr)
})

##########################
#       FUNCTIONS        #
##########################

log_message <- function(level, msg) {
  cat(paste0("[", level, "] ", msg, "\n"))
}

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
    return(list(RUN_DEP=RUN_DEP, RUN_INDEP=RUN_INDEP, RUN_MUT=RUN_MUT, RUN_DENS=RUN_DENS))
}

create_bed_reference <- function(data, output_dir, keyword) {
  bed_path <- file.path(output_dir, paste0("reference_", keyword, ".bed"))
  write_tsv(data[, 1:3], bed_path, col_names = FALSE)
  
  if (file.exists(bed_path) && file.size(bed_path) > 0) {
    log_message("INFO", paste(basename(bed_path), "successfully created"))
  }
}

create_igv_resources <- function(input_file, task_output_dir, keyword) {
  if (!file.exists(input_file)) {
    log_message("WARNING", paste("Input file missing, skipping:", input_file))
    return(NULL)
  }
  
  # [INFO] Creating task-specific directory
  if (!dir.exists(task_output_dir)) dir.create(task_output_dir, recursive = TRUE)
  
  rmd_name    <- paste0(keyword, "_report.rmd")
  batch_path  <- file.path(task_output_dir, paste0(keyword, ".batch"))
  rmd_path    <- file.path(task_output_dir, rmd_name)
  helper_path <- file.path(task_output_dir, "render_report.R") # Simplified name as it's inside its own folder
  
  # [INFO] Pre-cleanup
  files_to_clean <- c(batch_path, rmd_path, helper_path)
  for (f_old in files_to_clean) {
    if (file.exists(f_old)) file.remove(f_old)
  }

  hotspots <- unique(read_tsv(input_file, col_select = 1:3, show_col_types = FALSE))
  colnames(hotspots) <- c("chr", "start", "end")
  
  create_bed_reference(hotspots, task_output_dir, keyword)

  # [INFO] Initializing batch commands with snapshotDirectory header
  batch_commands <- c(paste0('snapshotDirectory "', task_output_dir, '"'))
  
 
  rmd_content <- c("---", paste0("title: 'IGV Report - ", keyword, "'"), "output: html_document", "---", "")
  
  for (i in seq_len(nrow(hotspots))) {
    padding <- (hotspots$end[i] - hotspots$start[i]) * 2
    batch_commands <- c(batch_commands, 
                        paste0("goto ", hotspots$chr[i], ":", (hotspots$start[i] - padding), "-", (hotspots$end[i] + padding)),
                        paste0("snapshot ", i, ".png"))
    rmd_content <- c(rmd_content, paste0("## ", i, " (", hotspots$chr[i], ")"), paste0("![](./", i, ".png)"), "")
  }
  
  writeLines(batch_commands, batch_path)
  writeLines(rmd_content, rmd_path)

  # Helper R script for rendering
  helper_content <- c(
    paste0('log_message <- function(level, msg) cat(paste0("[", level, "] ", msg, "\\n"))'),
    paste0('log_message("INFO", "Rendering HTML report for ', keyword, '...")'),
    paste0('rmarkdown::render("./', rmd_name, '", quiet = TRUE)'),
    'log_message("INFO", "Report rendered successfully. Cleaning up PNG files...")',
    'png_files <- list.files(pattern = "\\\\.png$")',
    'if (length(png_files) > 0) file.remove(png_files)',
    'log_message("INFO", "Cleanup complete.")'
  )
  writeLines(helper_content, helper_path)
  
  log_message("INFO", paste("Resources for", keyword, "saved in:", task_output_dir))
}

##########################
#          MAIN          #
##########################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("[ERROR] Usage: Rscript 10_igv_batch_rmd_generation.R <cell_line> <threshold> <min_count> <enlargement> <sample_name> [flags...]")
}

cell_line   <- args[1]
threshold   <- args[2]
min_count   <- args[3]
enlargement <- args[4]
sample_name <- args[5]
flags_in    <- if(length(args) >= 6) args[6:length(args)] else character(0)

config <- parse_flags(flags_in)
suffix <- paste0(threshold, min_count)

types_to_run <- c()
if (config$RUN_DEP) types_to_run <- c(types_to_run, "Dep")
if (config$RUN_INDEP) types_to_run <- c(types_to_run, "Indep")

modes_to_run <- c()
if (config$RUN_MUT) modes_to_run <- c(modes_to_run, "mut")
if (config$RUN_DENS) modes_to_run <- c(modes_to_run, "dens")

for (t_short in types_to_run) {
  for (m_short in modes_to_run) {
    
    mode_folder <- if(m_short == "mut") "nmutations" else "density"
    
    base_path <- file.path("/scratch/DepIndep_dataset", 
                           paste0("Threshold_", threshold), 
                           paste0("Enlargement_", enlargement),
                           paste0(t_short, "endent_hotspots"), 
                           sample_name,
                           paste0("Ranking_", t_short, "_hotspots"), 
                           mode_folder, "Filtered_hotspots")
    
    #
    output_parent_dir <- file.path(base_path, "IGV_Browser_files")
    
    tasks <- list(
      list(path = file.path(base_path, "Filtered_and_blacklist_top200", 
                            paste0("0_", sample_name, "_", t_short, "_", cell_line, "_Ranked_hotspots_", mode_folder, "_", suffix, "_200blackfilt.tsv")), 
           key = "Ranked_Top200"),
      list(path = file.path(base_path, "Filtered_and_blacklist_top200", 
                            paste0("0_", sample_name, "_", t_short, "_", cell_line, "_Deleted_hotspots_", mode_folder, "_", suffix, "_200blacklist.tsv")), 
           key = "Deleted_Top200"),
      list(path = file.path(base_path, 
                            paste0("0_", sample_name, "_", t_short, "_", cell_line, "_Ranked200hotspots_", mode_folder, "_", suffix, "_blackfilt.tsv")), 
           key = "Ranked_200Filtered")
    )
    
    for (task in tasks) {
      # Defining dynamic subfolder for each task
      task_dir <- file.path(output_parent_dir, paste0(sample_name, "_", t_short, "_", m_short, "_", task$key))
      create_igv_resources(task$path, task_dir, paste0(sample_name, "_", t_short, "_", m_short, "_", task$key))
    }
  }
}

log_message("INFO", "All IGV tasks completed.")