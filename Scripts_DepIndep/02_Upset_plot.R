#!/usr/bin/env Rscript

# --- Load Libraries ---
suppressPackageStartupMessages({
  if (!require("UpSetR")) stop("[ERROR] Package 'UpSetR' is not installed.")
  if (!require("grid")) stop("[ERROR] Package 'grid' is not installed.")
})

# --- 1. CONFIGURATION ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("[ERROR] Insufficient arguments.")

input_file  <- args[1]
output_dir  <- args[2]
threshold   <- args[3]
min_count   <- args[4]
enlargement <- args[5]

cat(paste0("[INFO] Parameters received: Input=", input_file, "\n"))

# --- 2. DATA PROCESSING ---
if (!file.exists(input_file)) stop("[ERROR] Input file not found.")
data <- read.delim(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Clean whitespace
data[] <- lapply(data, function(x) if(is.character(x)) trimws(x) else x)

metadata_cols <- c("region_id", "chrom", "start", "end")
sample_cols <- setdiff(colnames(data), metadata_cols)
n_samples <- length(sample_cols)

if (n_samples < 1) stop("[ERROR] No sample columns.")

# ROBUST CONVERSION to Binary (0/1)
sample_data <- data[, sample_cols, drop = FALSE]

# Logic: If x is NA, result is 0. If x is "+", result is 1. Else 0.
upset_data <- as.data.frame(lapply(sample_data, function(x) {
  ifelse(!is.na(x) & x == "+", 1L, 0L)
}))

# --- DEBUG: CHECK DATA ---
total_ones <- sum(upset_data, na.rm = TRUE)
cat(paste0("[INFO] Total intersections (1s) found in data: ", total_ones, "\n"))

if (total_ones == 0) {
    cat("[WARNING] DATA IS EMPTY! No '+' converted to 1. The plot will be blank.\n")
}

# --- COLORS ---
n_samples <- length(sample_cols)

if (n_samples == 3) {
  # CASE 1: Exactly 3 samples -> Specific colors
  set_colors <- c("#fff2ae", "#BEDDF1", "#FF8B94")
} else {
  # CASE 2: More than 3 samples -> Green gradient
  green_gradient_func <- colorRampPalette(c("#C1FFC1", "#006400"))
  set_colors <- green_gradient_func(n_samples)
}

# --- 3. PLOTTING (PNG) ---
output_path <- file.path(output_dir, paste0(enlargement, "_UpsetPlot_", threshold, "_", min_count, ".png"))
cat("[INFO] Generating UpSet PNG...\n")

png_type <- ifelse(capabilities("cairo"), "cairo", "Xlib")
png(output_path, width = 12, height = 8, units = "in", res = 300, bg = "white", type = png_type)

upset(
    upset_data,
    sets = sample_cols,
    keep.order = TRUE,
    mainbar.y.label = "Hotspots Intersections",
    sets.x.label = "Set Size",

    # --- COLORS ---
    sets.bar.color = set_colors,  # Variable colors per set
    main.bar.color = "#006400",   # Dark green for main bars
    matrix.color = "#006400",     # Dark green for matrix points
    # --------------
    
    text.scale = c(1.3, 1.3, 1, 1, 1.3, 1)
)


invisible(dev.off())

if (file.exists(output_path)) {
    cat(paste0("[INFO] Saved PNG to: '", output_path, "'\n"))
} else {
    cat("[ERROR] PNG file was not created.\n")
}