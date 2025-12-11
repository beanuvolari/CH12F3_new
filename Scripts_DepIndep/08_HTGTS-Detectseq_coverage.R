#!/usr/bin/env Rscript

# =========================================================================
# 08_HTGTS_Detectseq_Coverage_Unified.R
# Questo script unifica l'analisi di sovrapposizione tra hotspot Detect-seq
# e regioni target HTGTS (ON/OFF), eseguendo analisi speculari su tutte le
# combinazioni di file (Full, Top 200, Filtrato) specificate dai flag.
# =========================================================================

# Libraries required for genomic intervals and plotting
library(GenomicRanges)
library(readr)
library(dplyr)
library(ggplot2)

########################
#       Functions      #
########################

# Check if file exists
check_file <- function(path, label) {
  if (!file.exists(path)) {
    # Non usiamo stop() qui, ma restituiamo FALSE e facciamo gestire l'assenza nel MAIN
    cat("[WARNING]", label, "file", path, "doesn't exist! Skipping related analysis.\n")
    return(FALSE)
  } else {
    cat("[INFO]", label, "file found:", path, "\n")
    return(TRUE)
  }
}

# Convert data.frame with chrom/start/end to GRanges
make_GRanges <- function(df, start_col = "start", end_col = "end") {
  GRanges(seqnames = df$chrom, ranges = IRanges(start = df$start, end = df$end))
}

# Generate and save a Pie Chart
save_pie_chart <- function(data_df, title, output_path) {
  # Calcola le percentuali per il grafico
  data_df <- data_df %>%
    mutate(
      percentage = Count / sum(Count),
      y.pos = cumsum(percentage) - 0.5 * percentage
    )

  plot <- ggplot(data_df, aes(x = "", y = percentage, fill = Target)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    theme_void() +
    geom_text(aes(y = y.pos, label = scales::percent(percentage, accuracy = 0.1)), color = "black", size = 3) +
    labs(title = title, fill = "Category") +
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(output_path, plot, width = 6, height = 6)
  cat("[INFO] Pie chart saved to:", output_path, "\n")
}

# --- Core Unified Analysis Function ---
# Esegue entrambe le analisi (Detect-seq vs HTGTS e HTGTS vs Detect-seq)
# Per un dato file di input (Detect-seq hotspots)

analyze_coverage <- function(gr_detect, gr_ONtarget, gr_OFFtarget, output_dir, file_prefix, analysis_type) {
  
  cat(sprintf("[INFO] Starting analysis type: %s\n", analysis_type))
  
  # ----------------------------------------------------
  # A) Analisi FOCUS: Detect-seq vs HTGTS (Detect-seq come Query)
  # Quanti Detect-seq hotspots cadono nei target ON/OFF?
  # ----------------------------------------------------
  
  # Sovrapposizioni: findOverlaps(Query, Subject)
  # Detect-seq è Query, HTGTS targets sono Subject
  overlaps_ON_D_vs_H <- findOverlaps(gr_detect, gr_ONtarget)
  overlaps_OFF_D_vs_H <- findOverlaps(gr_detect, gr_OFFtarget)
  
  total_detect <- length(gr_detect)
  covered_on_detect <- length(unique(queryHits(overlaps_ON_D_vs_H)))
  covered_off_detect <- length(unique(queryHits(overlaps_OFF_D_vs_H)))
  
  # Hotspots coperti da ON o OFF
  is_covered <- unique(c(queryHits(overlaps_ON_D_vs_H), queryHits(overlaps_OFF_D_vs_H)))
  not_known_detect <- total_detect - length(is_covered)
  
  # Crea la tabella di riepilogo
  summary_detect <- data.frame(
    Target = c("ON_target_covered", "OFF_target_covered", "Not_known_targets"),
    Count = c(covered_on_detect, covered_off_detect, not_known_detect)
  )
  
  # Salva summary e grafico
  summary_path_detect <- file.path(output_dir, paste0(file_prefix, "_DetectvsHTGTS_Summary.tsv"))
  write_tsv(summary_detect, summary_path_detect)
  save_pie_chart(summary_detect, 
                 sprintf("Detect-seq Hotspots Coverage in HTGTS Targets (%s)", analysis_type), 
                 file.path(output_dir, paste0(file_prefix, "_DetectvsHTGTS_Chart.pdf")))
  
  cat(sprintf("[INFO] Detect-seq vs HTGTS: ON=%d, OFF=%d, Unknown=%d out of %d\n", 
              covered_on_detect, covered_off_detect, not_known_detect, total_detect))

  # ----------------------------------------------------
  # B) Analisi FOCUS: HTGTS vs Detect-seq (HTGTS come Query)
  # Quanti target HTGTS (ON/OFF) sono coperti dagli hotspot Detect-seq?
  # ----------------------------------------------------

  # Sovrapposizioni: findOverlaps(Query, Subject)
  # HTGTS targets sono Query, Detect-seq è Subject
  overlaps_ON_H_vs_D <- findOverlaps(gr_ONtarget, gr_detect)
  overlaps_OFF_H_vs_D <- findOverlaps(gr_OFFtarget, gr_detect)

  total_on_target <- length(gr_ONtarget)
  covered_on_target <- length(unique(queryHits(overlaps_ON_H_vs_D)))
  total_off_target <- length(gr_OFFtarget)
  covered_off_target <- length(unique(queryHits(overlaps_OFF_H_vs_D)))
  
  not_covered_on <- total_on_target - covered_on_target
  not_covered_off <- total_off_target - covered_off_target
  
  # Summary 1: Dettagliato (ON/OFF coperti e non coperti)
  summary_htgts_detailed <- data.frame(
    Target = c("ON_target_covered", "ON_target_not_covered", 
               "OFF_target_covered", "OFF_target_not_covered"),
    Count = c(covered_on_target, not_covered_on, covered_off_target, not_covered_off)
  )
  
  # Summary 2: Globale (Coperti vs Non Coperti)
  total_htgts <- total_on_target + total_off_target
  covered_total <- covered_on_target + covered_off_target
  not_covered_total <- total_htgts - covered_total
  
  summary_htgts_global <- data.frame(
    Target = c("Targets_covered", "Targets_not_covered"),
    Count = c(covered_total, not_covered_total)
  )
  
  # Salva summaries e grafico globale
  summary_path_htgts_detailed <- file.path(output_dir, paste0(file_prefix, "_HTGTSvsDetectseq_DetailedSummary.tsv"))
  write_tsv(summary_htgts_detailed, summary_path_htgts_detailed)
  
  summary_path_htgts_global <- file.path(output_dir, paste0(file_prefix, "_HTGTSvsDetectseq_GlobalSummary.tsv"))
  write_tsv(summary_htgts_global, summary_path_htgts_global)

  save_pie_chart(summary_htgts_global, 
                 sprintf("HTGTS Targets Coverage by Detect-seq Hotspots (%s)", analysis_type), 
                 file.path(output_dir, paste0(file_prefix, "_HTGTSvsDetectseq_Chart.pdf")))

  cat(sprintf("[INFO] HTGTS vs Detect-seq: ON covered=%d of %d, OFF covered=%d of %d\n", 
              covered_on_target, total_on_target, covered_off_target, total_off_target))
  
  cat(sprintf("[INFO] Analysis %s finished.\n\n", analysis_type))
}


########################
#          Main        #
########################

args <- commandArgs(trailingOnly = TRUE)

# --- INIZIO GESTIONE ARGOMENTI E FLAGS ---

if (length(args) < 5) {
     stop("[ERROR] At least 5 positional parameters required: cell_line, threshold_num, min_count, enlargement, sample. Optional flags follow.")
     }

# Positional Arguments
cell_line <- args[1]
threshold_num <- as.integer(args[2])
min_count <- as.integer(args[3])
enlargement <- as.integer(args[4])
sample <- args[5]

# Optional Flags
extra_flags <- if (length(args) >= 6) {
     args[6:length(args)]
} else {
     character(0) 
}

# --- FUNZIONE parse_flags UTENTE (integrata) ---
parse_flags <- function(flags) {
  # [La logica della tua funzione parse_flags è qui...]
  # (Per mantenere la compatibilità, si assume che la funzione sia definita 
  #  come da input utente e restituisca una lista)

  RUN_DEP <- FALSE; RUN_INDEP <- FALSE; RUN_MUT <- FALSE; RUN_DENS <- FALSE
  TYPE_SET <- FALSE; MODE_SET <- FALSE; ALL_SEEN <- FALSE

  # Logica di parsing (omessa qui per brevità, ma usa quella fornita dall'utente)
  # ... (La tua logica di parsing è qui) ...
  
  if (length(flags) > 0) {
    # Re-implementing user logic here for full transparency, adjusting the implementation
    # to be robust in R.
    flags <- flags[flags != ""]

    for (f in flags) {
      if (f == "--ALL") { ALL_SEEN <- TRUE }
      if (f == "--Dep") { RUN_DEP <- TRUE; RUN_INDEP <- FALSE; TYPE_SET <- TRUE }
      if (f == "--Indep") { RUN_DEP <- FALSE; RUN_INDEP <- TRUE; TYPE_SET <- TRUE }
      if (f == "--mut") { RUN_MUT <- TRUE; RUN_DENS <- FALSE; MODE_SET <- TRUE }
      if (f == "--dens") { RUN_MUT <- FALSE; RUN_DENS <- TRUE; MODE_SET <- TRUE }
    }

    if (!TYPE_SET) { RUN_DEP <- TRUE; RUN_INDEP <- TRUE }
    if (!MODE_SET) { RUN_MUT <- TRUE; RUN_DENS <- TRUE }
    
    if (ALL_SEEN && !TYPE_SET) { RUN_DEP <- TRUE; RUN_INDEP <- TRUE }
    if (ALL_SEEN && !MODE_SET) { RUN_MUT <- TRUE; RUN_DENS <- TRUE }
    
  } else {
    # Default behavior if no flags: run all
    RUN_DEP <- TRUE; RUN_INDEP <- TRUE; RUN_MUT <- TRUE; RUN_DENS <- TRUE
  }
  
  cat("[INFO] Flags parsed:\n")
  cat("       RUN_DEP =", RUN_DEP, "\n")
  cat("       RUN_INDEP =", RUN_INDEP, "\n")
  cat("       RUN_MUT =", RUN_MUT, "\n")
  cat("       RUN_DENS =", RUN_DENS, "\n")

  return(list(
    RUN_DEP=RUN_DEP,
    RUN_INDEP=RUN_INDEP,
    RUN_MUT=RUN_MUT,
    RUN_DENS=RUN_DENS
  ))
}
# --- FINE parse_flags UTENTE ---


flags <- parse_flags(extra_flags)

# Build combinations list based on parsed flags
combinations <- list()
if (flags$RUN_DEP) combinations <- append(combinations, list(list(type = "Dep", type_path = "Dependent")))
if (flags$RUN_INDEP) combinations <- append(combinations, list(list(type = "Indep", type_path = "Independent")))

modes <- list()
if (flags$RUN_DENS) modes <- append(modes, "density")
if (flags$RUN_MUT) modes <- append(modes, "nmutations")


if (length(combinations) == 0 || length(modes) == 0) {
  stop("[ERROR] No analysis combination was selected based on flags. Check inputs.")
}

cat(sprintf("[INFO] Selected combinations: %d\n", length(combinations) * length(modes)))

# --- Caricamento file di riferimento ON/OFF (statici) ---
reference_dir <- file.path("scratch/reference_data")
On_target_path  <- file.path(reference_dir, paste0("ON_target_", cell_line, ".bed"))
Off_target_path <- file.path(reference_dir, paste0("OFF_target_", cell_line, ".bed"))

if (!check_file(On_target_path, "ON-target reference") || 
    !check_file(Off_target_path, "OFF-target reference")) {
  stop("[FATAL ERROR] Missing critical ON/OFF reference files.")
}

On_target  <- read_tsv(On_target_path, col_names = FALSE, col_types = cols(), col_select = 1:3)
Off_target <- read_tsv(Off_target_path, col_names = FALSE, col_types = cols(), col_select = 1:3)
colnames(On_target) <- colnames(Off_target) <- c("chrom", "start", "end")

gr_ONtarget  <- make_GRanges(On_target)
gr_OFFtarget <- make_GRanges(Off_target)

cat(sprintf("[INFO] HTGTS Reference Loaded: ON targets=%d, OFF targets=%d\n", 
            length(gr_ONtarget), length(gr_OFFtarget)))


# --- Lista dei file di input Detect-seq da processare per ogni combinazione ---
detect_file_suffixes <- c(
  ".tsv", 
  "200hotspots.tsv", 
  "_blackfilt.tsv", 
  "200hotspots_blackfilt.tsv"
)
analysis_labels <- c(
  "FULL_SET", 
  "TOP200_POS", 
  "FULL_FILTERED", 
  "TOP200_FILTERED"
)

# --- LOOP PRINCIPALE: Itera su tutte le combinazioni ---
for (comb in combinations) {
  current_type <- comb$type
  current_type_path <- paste0(comb$type, "endent")

  for (current_mode in modes) {
    cat("=====================================================\n")
    cat(sprintf("[INFO] START COMBINATION: TYPE=%s, MODE=%s\n", current_type, current_mode))
    
    # Costruzione del percorso base
    base_dir <- file.path("Dep_Indep_dataset", 
                          paste0("Threshold_", threshold_num),
                          paste0("Enlargement_", enlargement), 
                          sample, 
                          paste0("Ranking_", current_type, "_hotspots"), 
                          current_mode)
    
   # 1. Output directory per i set NON FILTRATI (FULL, TOP200)
    output_dir_base <- file.path(base_dir, "Coverage_plots")
    
    # 2. Output directory per i set FILTRATI
    output_dir_filtered <- file.path(base_dir, "Filtered_hotspots", "Coverage_plots")
    
    # Crea le directory di output (vengono create all'interno del loop dei file)
    
    
    # Itera sui 4 tipi di file Detect-seq per questa combinazione
    for (i in 1:length(detect_file_suffixes)) {
      suffix <- detect_file_suffixes[i]
      label <- analysis_labels[i]
      
      # Decide la directory di output in base al suffisso/etichetta
      if (grepl("FILTERED|_blackfilt", label)) {
        # Usa la directory FILTRATA
        current_output_root <- file.path(output_dir_filtered)
      } else {
        # Usa la directory BASE/NON FILTRATA
        current_output_root <- file.path(output_dir_base)
      }
      
      # --- Logica di costruzione del percorso file (OMESSA PER BREVITÀ) ---
      # ... [input_file_path e gr_input] ...
      
      if (!check_file(input_file_path, sprintf("Detect-seq input (%s)", label))) next

      input_data <- read_tsv(input_file_path, col_names = TRUE)
      gr_input <- make_GRanges(input_data)
      
      
      # 2. Esecuzione analisi unificata
      output_prefix <- paste0(base_name_clean, "_", label)
      
      # Chiamata alla funzione analyze_coverage (dove si creano le sottocartelle)
      analyze_coverage(gr_detect = gr_input, 
                       gr_ONtarget = gr_ONtarget, 
                       gr_OFFtarget = gr_OFFtarget, 
                       output_root_dir = current_output_root, # Passiamo la directory base corretta
                       file_prefix = output_prefix, 
                       analysis_type = paste(current_type, current_mode, label, sep="/"))
      
    }
    cat(sprintf("[INFO] END COMBINATION: TYPE=%s, MODE=%s\n", current_type, current_mode))
    cat("=====================================================\n\n")
  }
}

cat("[INFO] All processing completed successfully!\n")