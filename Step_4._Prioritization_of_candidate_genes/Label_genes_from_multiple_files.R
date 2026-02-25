################################################################################################
## Topic: Label genes from multiple .txt files and export master tables                         ##
## Author: Olga Andrea Hernandez Miranda (Miranda H)                                            ##
## Date: 08/09/2025                                                                            ##                                                                          ##
## Description:                                                                               ##
## This script reads multiple TXT files (1 gene ID per line; optional header "Gene"),          ##
## assigns a label based on the filename (without extension), and exports:                     ##
##   (1) One CSV per label (Gene, Label)                                                       ##
##   (2) One combined CSV with all genes and labels                                            ##
##   (3) A simple report with gene counts per label                                            ##
################################################################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

# ==============================================================================================
# [1] CONFIGURATION
# ==============================================================================================

# Directory containing the .txt files
input_dir <- "."

# Output folder
out_dir <- file.path(input_dir, "salidas")

# Output files
combined_name <- "all_genes_with_labels.csv"
report_name   <- "reporte_etiquetas.csv"

# Explicit list of expected files (comment this block if you prefer automatic detection)
txt_files <- c(
  "Fer.txt","FT.txt","GED.txt","GNM.txt","OCDT.txt",
  "Pol.txt","PostPol.txt","PrePol.txt","RE.txt"
)

# Alternative: process all .txt files in input_dir
# txt_files <- list.files(input_dir, pattern = "\\.txt$", ignore.case = TRUE)

# ==============================================================================================
# [2] HELPERS
# ==============================================================================================

read_and_label <- function(filepath) {
  label <- tools::file_path_sans_ext(basename(filepath))
  
  lines <- tryCatch(
    readr::read_lines(filepath),
    error = function(e) {
      message("Could not read: ", filepath, " -> ", conditionMessage(e))
      character(0)
    }
  )
  
  if (length(lines) == 0) return(NULL)
  
  lines <- str_trim(lines)
  lines <- lines[nzchar(lines)]
  
  # Drop optional header (Gene / label)
  if (length(lines) > 0) {
    first <- lines[1]
    if (tolower(first) == "gene" || identical(first, label)) {
      lines <- lines[-1]
    }
  }
  
  # Remove stray "gene" lines
  lines <- lines[tolower(lines) != "gene"]
  
  tibble(
    Gene  = as.character(lines),
    Label = label
  ) %>%
    filter(!is.na(Gene), Gene != "") %>%
    distinct()
}

# ==============================================================================================
# [3] INPUT CHECKS
# ==============================================================================================

if (!dir.exists(input_dir)) stop("Input folder does not exist: ", input_dir)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

paths <- file.path(input_dir, txt_files)
paths <- paths[file.exists(paths)]

if (length(paths) == 0) {
  stop("No .txt files were found. Check input_dir and txt_files.")
}

message("Processing files:\n", paste0(" - ", basename(paths), collapse = "\n"))

# ==============================================================================================
# [4] PROCESS + EXPORT
# ==============================================================================================

label_dfs <- list()

for (f in paths) {
  df_lab <- read_and_label(f)
  
  if (is.null(df_lab) || nrow(df_lab) == 0) {
    warning("No valid genes in: ", basename(f))
    next
  }
  
  etiqueta <- unique(df_lab$Label)
  out_label_path <- file.path(out_dir, paste0(etiqueta, ".csv"))
  
  write_csv(df_lab, out_label_path)
  message(sprintf(" - %s: %d genes -> %s", etiqueta, nrow(df_lab), out_label_path))
  
  label_dfs[[etiqueta]] <- df_lab
}

if (length(label_dfs) == 0) stop("No valid data to combine.")

df_all <- bind_rows(label_dfs) %>%
  arrange(Label, Gene)

write_csv(df_all, file.path(out_dir, combined_name))

reporte <- df_all %>%
  count(Label, name = "n_genes") %>%
  arrange(desc(n_genes))

write_csv(reporte, file.path(out_dir, report_name))

message("Done.\n",
        " - Combined CSV: ", file.path(out_dir, combined_name), "\n",
        " - Report:       ", file.path(out_dir, report_name), "\n")