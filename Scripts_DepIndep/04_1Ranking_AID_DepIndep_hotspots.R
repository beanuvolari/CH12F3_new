#!/usr/bin/env Rscript

# Load required libraries
library(GenomicRanges)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)

######################
#      FUNCTIONS     #
######################

# Read and filter hotspot regions for the selected sample, add region_id and length
read_hotspots <- function(file, sample_col) {
  hotspots <- read_tsv(file, col_names = TRUE, col_select = -1, show_col_types = FALSE)
  # message("[DEBUG]: problems in file:") # Uncomment for debug
  # print(problems(hotspots))
  filtered <- hotspots[hotspots[[sample_col]] == "+", ]
  filtered %>%
    mutate(region_id = seq_len(n()),
           Length = end - start)
}

# Read mutation file and format columns
read_mutations <- function(file) {
  mut <- read_tsv(file, col_names = FALSE, col_types = cols(), col_select = 1:3, show_col_types = FALSE)
  colnames(mut) <- c("chrom", "start", "end")
  # message("[DEBUG]: problems in file:") # Uncomment for debug
  # print(problems(mut))
  mut
}

# Annotate hotspots with mutation counts and density
annotate_hotspots <- function(hotspots_df, mutations) {
  gr_hotspots <- GRanges(seqnames = hotspots_df$chrom,
                         ranges = IRanges(start = hotspots_df$start, end = hotspots_df$end))
  gr_mutations <- GRanges(seqnames = mutations$chrom,
                          ranges = IRanges(start = mutations$start, end = mutations$end))
  mutation_counts <- countOverlaps(gr_hotspots, gr_mutations)
  hotspots_df %>%
    mutate(Mutations = mutation_counts,
           Density = Mutations / Length)
}

# Output ranking tables based on density and mutation count
# ---- UPDATED FUNCTION ----
output_results <- function(hotspots_df, file_density, file_mutations) {

  # ----- DENSITY -----
  if (!is.null(file_density) && file_density != "") {
    hotspots_density <- hotspots_df %>%
      arrange(desc(Density), desc(Mutations)) %>%
      mutate(Rank = row_number())

    write_tsv(hotspots_density, file_density)
    message("[INFO] Density-ranked output saved to: ", file_density)
  } else {
    message("[INFO] Skipping density output (empty filename).")
  }

  # ----- MUTATIONS -----
  if (!is.null(file_mutations) && file_mutations != "") {
    hotspots_mutations <- hotspots_df %>%
      arrange(desc(Mutations), desc(Density)) %>%
      mutate(Rank = row_number())

    write_tsv(hotspots_mutations, file_mutations)
    message("[INFO] Mutation-count output saved to: ", file_mutations)
  } else {
    message("[INFO] Skipping mutation output (empty filename).")
  }

}

###################
#      MAIN       #
###################

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("[ERROR] You must provide at least 5 parameters: sample, hotspots_file, mutations_file, output_density, output_mutation.")
}

sample <- args[1]
hotspots_file <- args[2]
mutations_file <- args[3]
output_density <- args[4]
output_mutation <- args[5]

# Remove unwanted words from the sample names
short_sample <- sample %>%
  str_remove("^Primary_B_cells_") %>%
  str_remove("(_merged_sorted)?_CLEAN$")

# Read and filter hotspots
hotspots_df <- read_hotspots(hotspots_file, short_sample)
message("[INFO] Filtered hotspots loaded: ", nrow(hotspots_df), " regions.")
if (nrow(hotspots_df) == 0) {
  warning("[WARNING] No hotspots found for sample: ", short_sample)
}
#print(head(hotspots_df))   # Uncomment for debug

# Read mutations
mutations <- read_mutations(mutations_file)
message("[INFO] Mutations loaded: ", nrow(mutations), " records.")
if (nrow(mutations) == 0) {
  warning("[WARNING] No mutations found in file: ", mutations_file)
}
#print(head(mutations))      # Uncomment for debug

# [INFO] Annotate hotspots
hotspots_annotated <- annotate_hotspots(hotspots_df, mutations)
message("[INFO] Hotspots annotated with mutation counts.")
#print(head(hotspots_annotated))   # Uncomment for debug

# Output results
output_results(hotspots_annotated, output_density, output_mutation)
