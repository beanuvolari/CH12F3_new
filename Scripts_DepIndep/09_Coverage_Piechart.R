#!/usr/bin/env Rscript

# =========================================================================
# 33_8_Piechart_Unified.R
# Generates pie charts for Detect-seq vs HTGTS and HTGTS vs Detect-seq
# searching through Full and Filtered hotspot directories.
# =========================================================================

library(ggplot2)
library(readr)
library(dplyr)
library(ggrepel)

# Color palette definition
pie_levels <- c(
  "Not_known_targets",
  "Targets_not_covered",
  "OFF_targets_covered",
  "ON_targets_covered",
  "OFF_targets_not_covered",
  "ON_targets_not_covered"
)
my_palette <- c(
  "Not_known_targets"      = "#8dd3c7ff",
  "Targets_not_covered"    = "#8dd3c7ff",
  "OFF_targets_covered"     = "#528cf0dc",
  "ON_targets_covered"      = "#e98697ff",
  "OFF_targets_not_covered" = "#b3faffff",
  "ON_targets_not_covered"  = "#f7e99eff" 
)

########################
#      Functions       #
########################

generate_pie_chart <- function(summary_file_path, output_dir, enlargement, sample, type, cell_line, mode, threshold, mincount) {
  
  # [INFO] Check if file exists
  if (!file.exists(summary_file_path)) {
    cat(sprintf("[WARNING] File not found: %s. Skipping.\n", summary_file_path))
    return(NULL)
  }

  # Load data
  pie_data <- read_tsv(summary_file_path, col_types = cols())
  if (nrow(pie_data) == 0) {
    cat(sprintf("[WARNING] File %s is empty. Skipping.\n", summary_file_path))
    return(NULL)
  }

 # Determine context from filename and path
  is_200      <- grepl("200", basename(summary_file_path))
  is_filtered <- grepl("Filtered_hotspots", summary_file_path)
  is_htgts    <- grepl("HTGTS", basename(summary_file_path))
  is_detailed <- grepl("Detailed", basename(summary_file_path))

  # Build Title prefix
  title_prefix <- if (is_htgts) "HTGTS" else "Detect-seq"
  
  # [INFO] Title Logic based on folder (Full vs Filtered) and file name (200)
  if (is_filtered) {
    # Case: Filtered_hotspots folder
    title_sub <- if (is_200) "Filtered Hotspots top 200" else "Filtered Hotspots"
  } else {
    # Case: Full_hotspots folder
    title_sub <- if (is_200) "All hotspots top 200" else "All hotspots"
  }
  
  ggtitle_text <- sprintf("%s Target region coverage\n%s %s - %s\n%s", 
                          title_prefix, sample, type, cell_line, title_sub)

  # Prepare Data for plotting
  pie_data$Percentage <- round(pie_data$Count / sum(pie_data$Count) * 100, 1)
  pie_data$Target <- factor(pie_data$Target, levels = pie_levels)
  
  pie_data <- pie_data %>% 
    arrange(desc(Target)) %>% 
    mutate(
      fraction = Count / sum(Count),
      ymax = cumsum(fraction),
      ymin = c(0, head(ymax, n = -1)),
      labelPosition = (ymax + ymin) / 2,
      labelText = paste0(Percentage, "%")
    )

  # Plotting
  pie_plot <- ggplot(pie_data, aes(ymax = ymax, ymin = ymin, xmax = 1, xmin = 0, fill = Target)) +
    geom_rect() +
    coord_polar(theta = "y") +
    theme_void() +
    scale_fill_manual(values = my_palette) +
    labs(title = ggtitle_text) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))

  # Smart labeling logic
  large_labels <- pie_data %>% filter(Percentage >= 1.5)
  small_labels <- pie_data %>% filter(Percentage < 1.5)

  pie_plot <- pie_plot +
    geom_text(data = large_labels, aes(x = 1.1, y = labelPosition, label = labelText), size = 3.5)
  
  if (nrow(small_labels) > 0) {
    pie_plot <- pie_plot +
      geom_text_repel(data = small_labels, aes(x = 1.3, y = labelPosition, label = labelText),
                      nudge_x = 0.5, size = 3, segment.size = 0.2, inherit.aes = FALSE)
  }

  # Define Output filename
  # enlargement_sample_type_[Detailed]CoveragePie[200]_METHOD_mode_thresholdmincount.png
  detail_str <- if (is_detailed) "Detailed" else ""
  rank_str   <- if (is_200) "200" else ""
  method_str <- if (is_htgts) "HTGTS" else "Detect"
  
  out_filename <- sprintf("%s_%s_%s_%sCoveragePie%s_%s_%s_%d%d.png", 
                          enlargement, sample, type, detail_str, rank_str, method_str, mode, threshold, mincount)
  
  # Add blackfilt suffix if applicable
  if (is_filtered) out_filename <- gsub(".png", "_blackfilt.png", out_filename)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  ggsave(file.path(output_dir, out_filename), plot = pie_plot, width = 6, height = 6, bg = "white")
  cat(sprintf("[INFO] Saved: %s\n", out_filename))
}

########################
#        Main          #
########################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("[ERROR] Missing parameters: cell_line, threshold_num, min_count, enlargement, sample. Optional: --Dep, --Indep, --mut, --dens")
}

cell_line <- args[1]
threshold_num <- as.integer(args[2])
min_count <- as.integer(args[3])
enlargement <- as.integer(args[4])
sample <- args[5]

# Parse flags for Dep/Indep and mode (similar to your previous scripts)
extra_flags <- if (length(args) >= 6) args[6:length(args)] else character(0)
RUN_DEP <- "--Dep" %in% extra_flags || length(extra_flags) == 0
RUN_INDEP <- "--Indep" %in% extra_flags || length(extra_flags) == 0
RUN_MUT <- "--mut" %in% extra_flags || length(extra_flags) == 0
RUN_DENS <- "--dens" %in% extra_flags || length(extra_flags) == 0

types <- c(); if (RUN_DEP) types <- c(types, "Dep"); if (RUN_INDEP) types <- c(types, "Indep")
modes <- c(); if (RUN_MUT) modes <- c(modes, "nmutations"); if (RUN_DENS) modes <- c(modes, "density")


# Main loop
for (t in types) {
  type_folder <- paste0(t, "endent_hotspots")
  for (m in modes) {
    # Base directory for the combination
    path <- "/scratch/DepIndep_dataset"
    base_path <- file.path(path, paste0("Threshold_", threshold_num), 
                           paste0("Enlargement_", enlargement), type_folder, sample, 
                           paste0("Ranking_", t, "_hotspots"), m)
    
    # Iterate through Full and Filtered folders
    for (set_type in c("Full_hotspots", "Filtered_hotspots")) {
      
      # Target subdirectories where tables were saved
      methods <- c("Detect-seqvsHTGTS", "HTGTSvsDetect-seq")
      
      for (meth in methods) {
        table_dir <- file.path(base_path, set_type, "Pie_chart", meth)
        
        # In this directory, we look for ALL .tsv files
        summary_files <- list.files(table_dir, pattern = "\\.tsv$", full.names = TRUE)
        
        if (length(summary_files) == 0) {
          cat(sprintf("[WARNING] No summary tables found in %s. Skipping folder.\n", table_dir))
          next
        }
        
        for (s_file in summary_files) {
          generate_pie_chart(s_file, table_dir, enlargement, sample, t, cell_line, m, threshold_num, min_count)
        }
      }
    }
  }
}

cat("[INFO] Script finished successfully!\n")