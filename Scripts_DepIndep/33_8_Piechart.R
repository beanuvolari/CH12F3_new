#!/usr/bin/env Rscript

# nohup Rscript 33_8_Piechart.R 3 5 0 Dep UM_TKO_CIT_Duv_CLEAN nmutations HTGTSvsDetect-seq > nohup_350Dep_UMTKODuv_Piechart_HTGTS_blackfilt.out 2>&1 &
# nohup Rscript 33_8_Piechart.R 3 5 0 Dep UM_TKO_CIT_Duv_CLEAN nmutations Detect-seqvsHTGTS > nohup_350Dep_UMTKODuv_Piechart_Detect_blackfilt.out 2>&1 &
# nohup Rscript 33_8_Piechart.R 3 5 0 Dep UM_TKO_CIT_Duv_CLEAN nmutations Detect-seqvsHTGTS > nohup_350Dep_UMTKODuv_Piechart_Detect.out 2>&1 &

library(ggplot2)
library(readr)
library(dplyr)
library(ggrepel)

# Color palette definition
pie_levels <- c(
  "Not_known_targets",
  "Not_covered_targets",
  "OFF_target_covered",
  "ON_target_covered",
  "OFF_target_not_covered",
  "ON_target_not_covered"
)
my_palette <- c(
  "Not_known_targets"      = "#8dd3c7ff",
  "Not_covered_targets"    = "#8dd3c7ff",
  "OFF_target_covered"     = "#528cf0dc",
  "ON_target_covered"      = "#e98697ff",
  "OFF_target_not_covered" = "#b3faffff",
  "ON_target_not_covered"  = "#e4d2d5ff" 
)

########################
#      Functions       #
########################

generate_pie_chart <- function(summary_file_path, output_dir, prefix, sample, type, title_prefix, method = "simple") {
  # Check if the file exists
  if (!file.exists(summary_file_path)) {
    warning(paste("[WARNING] Missing file:", summary_file_path))
    return(NULL)
  } else {
    cat("[INFO] Summary file found:", summary_file_path, "\n")
  }

  # Create the output directory if it does not exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("[INFO] Output directory created:", output_dir, "\n")
  }

  # Read the summary file
  pie_data <- read_tsv(summary_file_path, col_types = cols())

  # Calculate percentages
  pie_data$Percentage <- round(pie_data$Count / sum(pie_data$Count) * 100, 1)

  # Determine if the file is for the "200" version
  title_suffix <- if (grepl("200", basename(summary_file_path))) {
    "Top 200 ranked hotspots"
  } else {
    "All hotspots"
  }

  # Create dynamic title
  ggtitle_text <- sprintf("%s\nTarget region coverage:\n%s %s\n%s", 
                          title_prefix, sample, type, title_suffix)

  # Set factor levels for consistent coloring
  pie_data$Target <- factor(pie_data$Target, levels = pie_levels)

  # Calculate label positions
  pie_data <- pie_data %>% 
    arrange(desc(Target)) %>% 
    mutate(
      fraction = Count / sum(Count),
      ymax = cumsum(fraction),
      ymin = c(0, head(ymax, n = -1)),
      labelPosition = (ymax + ymin) / 2,
      labelText = paste0(Percentage, "%")
    )

  # Base plot
  pie_plot <- ggplot(pie_data, aes(ymax = ymax, ymin = ymin, xmax = 1, xmin = 0, fill = Target)) +
    geom_rect() +
    coord_polar(theta = "y") +
    theme_void() +
    ggtitle(ggtitle_text) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = my_palette)

  # Labeling method
  if (method == "simple") {
    pie_plot <- pie_plot +
      geom_text(aes(
        x = 1.2, 
        y = labelPosition,
        label = labelText
      ), size = 4)
  } else if (method == "smart") {
    large_labels <- pie_data %>% filter(Percentage >= 1.5)
    small_labels <- pie_data %>%
      filter(Percentage < 1.5) %>%
      mutate(
        angle = 360 * labelPosition,
        side = ifelse(angle > 90 & angle < 270, "left", "right"),
        x_pos = ifelse(side == "left", -1.2, 1.2),
        nudge_x = ifelse(side == "left", -0.5, 0.5),
        hjust_val = ifelse(side == "left", 1, 0)
      )

    pie_plot <- pie_plot +
      geom_text(data = large_labels, aes(
        x = 1.1, y = labelPosition, label = labelText
      ), size = 4, inherit.aes = FALSE) +
      geom_text_repel(
        data = small_labels,
        aes(x = x_pos, y = labelPosition, label = labelText),
        nudge_x = small_labels$nudge_x,
        hjust = small_labels$hjust_val,
        size = 4,
        segment.size = 0.3,
        direction = "x",
        inherit.aes = FALSE,
        show.legend = FALSE
      )
  }

  # Save the plot
  output_file <- file.path(output_dir, prefix)
  ggsave(output_file, plot = pie_plot, width = 6, height = 6, bg = "white")
  cat("[INFO] Plot saved as:", output_file, "\n")
}

########################
#        Main          #
########################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 7) {
  stop("[ERROR]: provide at least 6 parameters (threshold_num, min_count, enlargement, type, sample, dens_mut, comparison_method [HTGTSvsDetect-seq or Detect-seqvsHTGTS]).")
}

threshold_num <- as.integer(args[1])
min_count <- as.integer(args[2])
enlargement <- as.integer(args[3])
type <- args[4]
sample <- args[5]
dens_mut <- args[6]
comparison_method <- args[7]

file_dir <- file.path("Dep_Indep_dataset", 
                      paste0("Threshold_", threshold_num),
                      paste0("Enlargement_", enlargement), sample, 
                      paste0("Ranking_", type, "_hotspots"), dens_mut)

# Decide parameters based on the comparison method
if (comparison_method == "HTGTSvsDetect-seq") {
  title_prefix <- "HTGTS"

  file_dir_full <- file.path(file_dir, "HTGTSvsDetect-seq", "FULL")
  file_dir_200 <- file.path(file_dir, "HTGTSvsDetect-seq", "first_200")

  output_dir_full <- file.path(file_dir_full, "Pie_chart")
  output_dir_200 <- file.path(file_dir_200, "Pie_chart")

  prefix_full   <- sprintf("%s_%s_%s_%s_TargetCoveragePieFULL_HTGTS_%d_%d.png", enlargement, sample, type, dens_mut, threshold_num, min_count)
  prefix_full2  <- sprintf("%s_%s_%s_%s_2TargetCoveragePieFULL_HTGTS_%d_%d.png", enlargement, sample, type, dens_mut, threshold_num, min_count)
  prefix_200    <- sprintf("%s_%s_%s_%s_TargetCoveragePie200_HTGTS_%d_%d.png", enlargement, sample, type, dens_mut, threshold_num, min_count)
  prefix_200_2  <- sprintf("%s_%s_%s_%s_2TargetCoveragePie200_HTGTS_%d_%d.png", enlargement, sample, type, dens_mut, threshold_num, min_count)
  files <- setNames(
    c(prefix_full, prefix_full2, prefix_200, prefix_200_2),
    c(
      file.path(file_dir_full, sprintf("%s_%s_%s_%s_TargetCoverageSummaryFULL_HTGTS_%d%d.tsv", enlargement, sample, type, dens_mut, threshold_num, min_count)),
      file.path(file_dir_full, sprintf("%s_%s_%s_%s_2TargetCoverageSummaryFULL_HTGTS_%d%d.tsv", enlargement, sample, type, dens_mut, threshold_num, min_count)),
      file.path(file_dir_200, sprintf("%s_%s_%s_%s_TargetCoverageSummary200_HTGTS_%d%d.tsv", enlargement, sample, type, dens_mut, threshold_num, min_count)),
      file.path(file_dir_200, sprintf("%s_%s_%s_%s_2TargetCoverageSummary200_HTGTS_%d%d.tsv", enlargement, sample, type, dens_mut, threshold_num, min_count))
    )
  )
  method <- "smart"

} else if (comparison_method == "Detect-seqvsHTGTS") {
  title_prefix <- "Detect-seq"
  file_dir_full <- file.path(file_dir, "Detect-seqvsHTGTS", "FULL")
  file_dir_200 <- file.path(file_dir, "Detect-seqvsHTGTS", "first_200")
  
  output_dir_full <- file.path(file_dir_full, "Pie_chart")
  output_dir_200 <- file.path(file_dir_200, "Pie_chart")

  prefix_full   <- sprintf("%s_%s_%s_%s_TargetCoveragePieFULL_Detect_%d_%d.png", enlargement, sample, type, dens_mut, threshold_num, min_count)
  prefix_200    <- sprintf("%s_%s_%s_%s_TargetCoveragePie200_Detect_%d_%d.png", enlargement, sample, type, dens_mut, threshold_num, min_count)
  files <- setNames(
    c(prefix_full, prefix_200),
    c(
      file.path(file_dir_full, sprintf("%s_%s_%s_%s_TargetCoverageSummaryFULL_Detect_%d%d.tsv", enlargement, sample, type, dens_mut, threshold_num, min_count)),
      file.path(file_dir_200, sprintf("%s_%s_%s_%s_TargetCoverageSummary200_Detect_%d%d.tsv", enlargement, sample, type, dens_mut, threshold_num, min_count))
    )
  )
  method <- "smart"
} else {
  stop("[ERROR] comparison_method must be 'HTGTSvsDetect-seq' or 'Detect-seqvsHTGTS'")
}

# Generate the plots
for (f in names(files)) {
  # Choose output directory based on file type
  if (grepl("FULL", f)) {
    out_dir <- output_dir_full
  } else {
    out_dir <- output_dir_200
  }
  generate_pie_chart(f, out_dir, files[[f]], sample, type, title_prefix, method)
}

cat("[INFO] The script has finished successfully!\n")
