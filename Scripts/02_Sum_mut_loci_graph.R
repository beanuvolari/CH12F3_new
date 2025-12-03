#!/usr/bin/env Rscript

# Load required library
library(ggplot2)

###########################################
#              FUNCTIONS                  #
###########################################

# Function to extract threshold from filename 
extract_threshold <- function(name) {
  # Extract number after the last underscore and before ".pmat"
  match <- regmatches(name, regexpr("_(\\d+)\\.pmat$", name, perl=TRUE))
  if(length(match) == 0) return(NA_real_)
  as.numeric(sub("_(\\d+)\\.pmat", "\\1", match, perl=TRUE))
}

# Plotting function 
plot_metric <- function(x, y, output_path, plot_title, y_label) {
  df <- data.frame(x = x, y = y)

  plot <- ggplot(df, aes(x = as.factor(x), y = y)) +
    geom_line(group = 1, color = "red", linewidth = 1) +
    geom_point(size = 2, color = "blue") +
    geom_text(aes(label = y), vjust = -1, size = 4) +
    ggtitle(plot_title) +
    xlab("Thresholds") +
    ylab(y_label) +
    coord_cartesian(ylim = c(0, max(y))) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.8),
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )

  ggsave(output_path, plot = plot, width = 8, height = 6, bg = "white")
}

###########################################
#                 MAIN                    #
###########################################

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript 02_Sum_mut_loci_graph.R <input_file> <output_dir> <start_threshold>")
}

input_file <- args[1]
output_dir <- args[2]
start_threshold <- as.numeric(args[3])
filebase <- args[4]

# Load and process data
data <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Extract thresholds from filenames
thresholds <- sapply(data[, 1], extract_threshold)

valid_indices <- !is.na(thresholds)
if (sum(valid_indices) == 0) {
  stop("[ERROR] No valid thresholds found in the input file.")
}
data <- data[valid_indices, , drop=FALSE]
thresholds <- thresholds[valid_indices]

# Filter rows by start_threshold
data_filtered <- data[thresholds >= start_threshold, , drop = FALSE]
thresholds_filtered <- thresholds[thresholds >= start_threshold]

# Stop if no data remains after filtering
if (nrow(data_filtered) == 0) {
  stop(paste("[ERROR] No data points found with threshold >=", start_threshold, ". Please provide another threshold."))
}

# Prepare data for plotting
x <- thresholds_filtered
y <- data_filtered[, 2]

# Determine plot label and title
if (grepl("mutnum", input_file)) {
  ylab <- "Number of mutations"
  plot_type <- "Mutations"
} else if (grepl("locnum", input_file)) {
  ylab <- "Number of loci"
  plot_type <- "Loci"
} else {
  ylab <- "Value"
  plot_type <- "Metric"
}

# Build output filename
base_name <- tools::file_path_sans_ext(basename(input_file))
output_file <- file.path(output_dir, paste0(base_name, "_from_", start_threshold, ".png"))


# Plot the data
plot_metric(
  x = x,
  y = y,
  output_path = output_file,
  plot_title = paste(filebase, "sum of", plot_type, "from threshold:", start_threshold),
  y_label = ylab
)
