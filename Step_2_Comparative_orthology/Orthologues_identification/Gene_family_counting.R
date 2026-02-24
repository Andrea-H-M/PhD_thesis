################################################################################################
## Topic: Arabidopsis column extraction + gene-family counting                                ##
## Author: Olga Andrea Hernandez Miranda (Miranda H)                                          ##
## Date: 09/03/2024                                                                           ##                                         ##
## Description:                                                                               ##
## This script:                                                                               ##
##  (1) Reads a “master table” (Orthogroups/TablaMaestra-like).                               ##
##  (2) Extracts the Arabidopsis column.                                                      ##
##  (3) Cleans entries, splits multi-IDs, extracts gene names (before "_").                   ##
##  (4) Assigns each gene to a canonical family (canonical list + aliases).                   ##
##  (5) Produces a family count table: number of genes per family.                            ##
################################################################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
})

# ==============================================================================================
# 1) Paths & setup
# ==============================================================================================

# Working directory (optional)
# setwd("C:/Users/andii/OneDrive/Documents/DoctoradoEnCiencias/Proyecto/Tutorial5/Results_Sep07/Orthogroups")

# Input table (uploaded / your path)
input_table <- "C:/Users/andii/OneDrive/Documents/DoctoradoEnCiencias/Proyecto/Tutorial5/Results_Sep07/Orthogroups/TablaMaestra.csv"
# If you want to use the uploaded file directly, replace with your local path after you download it.
# input_table <- "325707bf-40b3-4a00-b734-be1d4c1a0672.csv"

# Canonical families file (must contain a column named: "Familia")
canon_families_file <- "familias_canonicas_planttfdb_like.csv"

# Outputs
out_arab_genes      <- "Arabidopsis_genes_extracted.csv"
out_arab_fam_counts <- "Arabidopsis_family_counts.csv"
out_unassigned      <- "Arabidopsis_unassigned_genes.csv"

# ==============================================================================================
# 2) Helpers
# ==============================================================================================

clean_id <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace_all("^\\.+", "") %>%
    str_replace_all("\\s+", "") %>%
    toupper()
}

# Extract gene name from IDs like "AGO10_ARATH" -> "AGO10"
extract_gene_name <- function(id) {
  id <- clean_id(id)
  # keep only the part before the first underscore
  str_split_fixed(id, "_", n = 2)[, 1]
}

read_canon_families <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    transmute(Familia = toupper(str_replace_all(Familia, "\\s+", ""))) %>%
    distinct() %>%
    pull(Familia)
}

assign_family <- function(gene, familias, aliases) {
  if (is.na(gene) || gene == "") return(NA_character_)
  
  # 1) canonical prefix match (choose the LONGEST)
  cand <- familias[str_detect(gene, paste0("^", familias))]
  if (length(cand) > 0) return(cand[which.max(nchar(cand))])
  
  # 2) alias prefix match (choose the LONGEST alias)
  keys <- names(aliases)
  a_cand <- keys[str_detect(gene, paste0("^", keys))]
  if (length(a_cand) > 0) {
    best_alias <- a_cand[which.max(nchar(a_cand))]
    return(unname(aliases[[best_alias]]))
  }
  
  NA_character_
}

# ==============================================================================================
# 3) Read table and detect Arabidopsis column
# ==============================================================================================

tabla <- read_csv(input_table, show_col_types = FALSE)

# Try to find an Arabidopsis column robustly
# (matches: "Arabidopsis", "ARATH", "A_thaliana", etc.)
arab_col <- names(tabla)[str_detect(names(tabla), regex("arabidopsis|arath|thaliana|a\\.thaliana", ignore_case = TRUE))]

if (length(arab_col) == 0) {
  stop("No Arabidopsis column detected. Rename the Arabidopsis column or adjust the regex in arab_col.")
}
arab_col <- arab_col[1]  # if multiple matches, take the first

# ==============================================================================================
# 4) Extract Arabidopsis IDs -> split -> extract gene names
# ==============================================================================================

arab_genes <- tabla %>%
  transmute(Arabidopsis_raw = .data[[arab_col]]) %>%
  mutate(Arabidopsis_raw = as.character(Arabidopsis_raw)) %>%
  filter(!is.na(Arabidopsis_raw), Arabidopsis_raw != "") %>%
  # If cells contain multiple IDs separated by ";" (common in orthogroup tables)
  separate_rows(Arabidopsis_raw, sep = ";") %>%
  mutate(
    Arabidopsis_ID   = clean_id(Arabidopsis_raw),
    Arabidopsis_Gene = extract_gene_name(Arabidopsis_ID)
  ) %>%
  filter(Arabidopsis_ID != "", Arabidopsis_Gene != "") %>%
  distinct(Arabidopsis_ID, Arabidopsis_Gene)

write_csv(arab_genes, out_arab_genes)

# ==============================================================================================
# 5) Family assignment + counts
# ==============================================================================================

# Aliases (edit to match your biology / naming conventions)
aliases <- c(
  # Common TF shorthand seen in datasets
  "WRK"   = "WRKY",
  "MB"    = "MYB",
  "MYBS"  = "MYB",
  "MYBA"  = "MYB",
  "MYBL"  = "MYB",
  "MYBP"  = "MYB",
  "BH"    = "BHLH",
  "BZP"   = "BZIP",
  "TGA"   = "BZIP",
  "TGAL"  = "BZIP",
  "HDZIP" = "HD-ZIP"
)

familias_canon <- read_canon_families(canon_families_file)

arab_classified <- arab_genes %>%
  mutate(
    Familia = map_chr(Arabidopsis_Gene, assign_family,
                      familias = familias_canon,
                      aliases  = aliases)
  )

# Keep unassigned to curate later
unassigned <- arab_classified %>% filter(is.na(Familia) | Familia == "")
if (nrow(unassigned) > 0) write_csv(unassigned, out_unassigned)

# Count genes per family (unique genes)
family_counts <- arab_classified %>%
  filter(!is.na(Familia), Familia != "") %>%
  distinct(Arabidopsis_Gene, Familia) %>%
  count(Familia, name = "No. of genes", sort = TRUE)

print(family_counts, n = Inf)
write_csv(family_counts, out_arab_fam_counts)