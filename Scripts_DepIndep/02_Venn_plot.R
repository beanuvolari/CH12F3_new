#!/usr/bin/env Rscript

# --- Load Libraries ---
suppressPackageStartupMessages({
  if (!require("VennDiagram")) stop("[ERROR] Package 'VennDiagram' is not installed.")
  if (!require("grid")) stop("[ERROR] Package 'grid' is not installed.")
})

# Disable logs
invisible(futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))

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

data[] <- lapply(data, function(x) if(is.character(x)) trimws(x) else x)

metadata_cols <- c("region_id", "chrom", "start", "end")
sample_cols <- setdiff(colnames(data), metadata_cols)
n_samples <- length(sample_cols)

if (n_samples < 1) stop("[ERROR] No sample columns.")

# CHECK LIMITS
if (n_samples > 5) {
  cat(paste0("[WARNING] Detected ", n_samples, " sets. VennDiagram supports a maximum of 5 sets.\n"))
  cat("[WARNING] Skipping Venn plot generation for this dataset.\n")
  quit(save = "no", status = 0) # Exit cleanly without error
}

# Prepare List
venn_list <- list()
for (col in sample_cols) {
  # Robust check: value is "+" AND not NA
  venn_list[[col]] <- which(!is.na(data[[col]]) & data[[col]] == "+")
}

# Colors
if (n_samples <= 3) {
  my_fills <- c("#FF8B94", "#BEDDF1", "#fff2ae")[1:n_samples]
  cat_cols <- c("#FF8B94", "#BEDDF1", "#F6CF71")[1:n_samples]
} else {
  my_fills <- rainbow(n_samples, alpha = 0.7)
  cat_cols <- rainbow(n_samples)
}

# --- 3. PLOTTING (PNG) ---
output_path <- file.path(output_dir, paste0(enlargement, "_VennPlot_", threshold, "_", min_count, ".png"))
cat("[INFO] Generating Venn PNG...\n")

png(output_path, width = 8, height = 8, units = "in", res = 300, bg = "white")

tryCatch({
  venn_obj <- venn.diagram(
    x = venn_list,
    filename = NULL, 
    fill = my_fills,
    alpha = 0.7,
    lwd = 2,
    lty = "blank",
    cex = 1.0,
    fontfamily = "sans",
    fontface = "bold",
    cat.cex = 1.2,
    cat.col = cat_cols,
    cat.fontfamily = "sans",
    cat.fontface = "bold",
    margin = 0.1
  )
  
  grid.draw(venn_obj)
  cat(paste0("[INFO] Saved PNG to: '", output_path, "'\n"))
}, error = function(e) {
  cat(paste0("[ERROR] Venn generation failed: ", e$message, "\n"))
})

invisible(dev.off())