#!/usr/bin/env Rscript

# --- Load Libraries ---
suppressPackageStartupMessages({
  if (!require("eulerr")) stop("[ERROR] Package 'eulerr' is not installed.")
})

# --- 1. CONFIGURATION ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("[ERROR] Insufficient arguments. Required: InputFile, OutputDir, Threshold, MinCount, Enlargement")
}

input_file  <- args[1]
output_dir  <- args[2]
threshold   <- args[3]
min_count   <- args[4]
enlargement <- args[5]

cat(paste0("[INFO] Parameters received: Input=", input_file, "\n"))

# --- 2. DATA PROCESSING ---

if (!file.exists(input_file)) stop(paste0("[ERROR] Input file not found: ", input_file))

# Read Data
data <- read.delim(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Clean whitespace
data[] <- lapply(data, function(x) if(is.character(x)) trimws(x) else x)

# Identify Sample Columns
metadata_cols <- c("region_id", "chrom", "start", "end")
sample_cols <- setdiff(colnames(data), metadata_cols)

if (length(sample_cols) < 1) stop("[ERROR] No sample columns found.")
cat(paste0("[INFO] Detected ", length(sample_cols), " sample columns.\n"))

# Convert "+" to Boolean Matrix
sample_data <- data[, sample_cols, drop = FALSE]

# ROBUST CONVERSION: Handle NAs explicitly
# If cell is NA, treat as FALSE. If "+", treat as TRUE.
data_logical <- sample_data == "+"
data_logical[is.na(data_logical)] <- FALSE

# Filter empty rows
non_empty_rows <- rowSums(data_logical) > 0
data_logical <- data_logical[non_empty_rows, , drop = FALSE]

# --- 3. PLOTTING ---

cat("[INFO] Calculating intersections...\n")

# Use tryCatch because euler() might fail with high dimensions (e.g. 7 sets) if complex
tryCatch({
  fit <- euler(data_logical)
  
  # Colors
  n_samples <- length(sample_cols)
  
  if (n_samples <= 3) {
    # CASE 1: if samples are 3 or less, use PASTEL colors
    my_fills <- c("#FF8B94", "#BEDDF1", "#fff2ae")[1:n_samples]
  } else {
    # CASE 2: if samples are more than 3, create a green gradient
    #green_gradient_func <- colorRampPalette(c("#006400", "#2E8B57", "#C1FFC1"))
    # Colors generation
    #my_fills <- green_gradient_func(n_samples)

    my_fills <- hcl(h = seq(15, 375, length = n_samples + 1), l = 65, c = 100)[1:n_samples]
  }

  # --- STYLE SETTINGS (UNIFIED) ---
  # Change these values to affect BOTH the numbers and the legend
  font_style <- 2       # 1 = Plain, 2 = Bold
  font_size  <- 1.0     # Text size (cex)
  text_color <- "#000000" # Text color
  
  # --- 4. SAVING (PNG) ---
  output_filename <- paste0(enlargement, "_EulerPlot_", threshold, "_", min_count, ".png")
  output_path <- file.path(output_dir, output_filename)
  
  # Open PNG device
  png_type <- ifelse(capabilities("cairo"), "cairo", "Xlib")
  png(output_path, width = 10, height = 8, units = "in", res = 300, bg = "white", type = png_type)
  
  print(plot(fit,
             fills = list(fill = my_fills, alpha = 0.7),
             edges = list(col = "transparent", lwd = 0), 
             
             # 1. QUANTITIES (Numbers inside the plot)
             quantities = list(
               font = font_style, 
               cex = font_size, 
               col = text_color, 
               fontfamily = "sans"
             ),
             
             labels = FALSE, # Turn off labels inside circles
             
             # 2. LEGEND (Sample names on the right)
             legend = list(
               side = "right", 
               font = font_style, 
               cex = font_size, 
               col = text_color, 
               alpha = 1,
               border = "transparent"
             )
  ))
  
  invisible(dev.off()) 
  cat(paste0("[INFO] Saved PNG to: '", output_path, "'\n"))
  
}, error = function(e) {
  cat(paste0("[ERROR] Euler calculation failed (likely too many sets or complex overlaps): ", e$message, "\n"))
})