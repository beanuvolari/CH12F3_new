library(GenomicRanges)
library(readr)
library(dplyr)
library(stringr)

# --- Modular functions ---

# Reads a BED file into a GRanges object
read_bed_to_granges <- function(filepath) {
  if (!file.exists(filepath)) {
    stop(paste("File not found:", filepath))
  }
  df <- read_tsv(filepath, col_names = FALSE, col_types = cols(), col_select = 1:3)
  colnames(df) <- c("chrom", "start", "end")
  GRanges(seqnames = df$chrom, ranges = IRanges(start = df$start, end = df$end))
}

# Enlarges GRanges regions according to the enlargement value
enlarge_region <- function(gr, enlargement) {
  enlargement_values <- c(0, 1000, 2000, 4000)
  if (!(enlargement %in% 0:3)) {
    stop("Error: Enlargement value should be 0, 1, 2, or 3.")
  }
  enlargement_bp <- enlargement_values[enlargement + 1]
  new_start <- pmax(start(gr) - enlargement_bp, 1)
  new_end <- end(gr) + enlargement_bp
  ranges(gr) <- IRanges(start = new_start, end = new_end)
  gr
}

# Computes presence matrix of overlaps between merged regions and GRanges list
compute_presence_matrix <- function(merged_regions, gr_list, sample_names) {
  results <- data.frame(
    region_id = seq_along(merged_regions),
    chrom = as.character(seqnames(merged_regions)),
    start = start(merged_regions),
    end = end(merged_regions)
  )
  for (i in seq_along(gr_list)) {
    presence_vec <- vapply(seq_along(merged_regions), function(idx) {
      region <- merged_regions[idx]
      overlaps <- countOverlaps(region, gr_list[[i]])
      overlaps > 0
    }, logical(1))
    col_name <- sample_names[i]
    results[[col_name]] <- ifelse(presence_vec, "+", "-")
  }
  results
}


# Writes the results dataframe to a file
write_results <- function(results, output_file) {
  write.table(results, output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

# --- Main ---

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript script.R <input_file1> [input_file2 ...] <output_file> <enlargement>")
}

output_file <- args[length(args) - 1]
enlargement <- as.numeric(args[length(args)])
input_files <- args[1:(length(args) - 2)]

# Extract sample names from the input files
file_names <- basename(input_files)
sample_names <- str_remove(file_names, "_merged_CTGA_.*")

# Remove unwanted words from the sample names
sample_names <- sample_names %>%
  str_remove("^Primary_B_cells_") %>%
  str_remove("(_merged_sorted)?_CLEAN$")
print("[DEBUG] Short sample names:")
print(sample_names)

# Read and convert input files to GRanges, then enlarge regions
print("[INFO] Loading and enlarging regions...")

gr_list <- lapply(input_files, read_bed_to_granges)
print("[DEBUG] First gr_list:")
print(class(gr_list))
print(length(gr_list))

gr_list <- lapply(gr_list, enlarge_region, enlargement = enlargement)
print("[DEBUG] Second gr_list:")
print(class(gr_list))
print(length(gr_list))


# Combine and merge all regions
print("[INFO] Merging regions...")
all_gr <- do.call(c, gr_list)
merged_regions <- reduce(sort(all_gr))
print("[DEBUG] merged regions:")
print(merged_regions)

# Compute presence matrix with sample names as column headers
print("[INFO] Generating the database...")
results <- compute_presence_matrix(merged_regions, gr_list, sample_names)

# Save the results to the output file
write_results(results, output_file)
