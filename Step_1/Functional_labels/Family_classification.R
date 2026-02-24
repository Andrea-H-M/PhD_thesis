################################################################################################
## TF & Epigenetic Regulator Family Classification and Counting                               ##
## Author: Olga Andrea Hernández Miranda, Miranda H                                           ##                                              ##
## Date: 23/02/2026                                                                           ##                                                                              ##
## Description:                                                                               ##
## This script classifies transcription factors (TF) and epigenetic regulators (ER) into      ##
## canonical gene families based on ID prefix matching and computes their family distribution.##
## Unassigned genes are also exported for manual curation.                                    ##
################################################################################################

allowWGCNAThreads()

setwd("C:/Users/andii/OneDrive/Documents/DoctoradoEnCiencias/Proyecto/Tutorial7/4FT")

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(purrr)
})

# ==============================================================================================
# General utilities
# ==============================================================================================

read_genes_1col <- function(path) {
  df <- tryCatch(
    read_csv(path, show_col_types = FALSE),
    error = function(e) read_csv(path, col_names = "Gene", show_col_types = FALSE)
  )
  if (!"Gene" %in% names(df)) names(df)[1] <- "Gene"
  df %>% select(Gene)
}

clean_id <- function(x) {
  x %>%
    str_trim() %>%
    str_replace_all("^\\.+", "") %>%
    str_replace_all("\\s+", "") %>%
    toupper()
}

read_canon_families <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    transmute(Familia = toupper(str_replace_all(Familia, "\\s+", ""))) %>%
    distinct() %>%
    pull(Familia)
}

assign_family <- function(id, families, aliases) {
  if (is.na(id) || id == "") return(NA_character_)
  
  cand <- families[str_detect(id, paste0("^", families))]
  if (length(cand) > 0) return(cand[which.max(nchar(cand))])
  
  keys <- names(aliases)
  a_cand <- keys[str_detect(id, paste0("^", keys))]
  if (length(a_cand) > 0) {
    best_alias <- a_cand[which.max(nchar(a_cand))]
    return(unname(aliases[[best_alias]]))
  }
  
  NA_character_
}

run_family_count <- function(genes_path, canon_path, aliases, out_counts, out_unassigned) {
  
  df_raw   <- read_genes_1col(genes_path)
  families <- read_canon_families(canon_path)
  
  df_classified <- df_raw %>%
    transmute(Gene = clean_id(Gene)) %>%
    mutate(Family = map_chr(
      Gene, assign_family,
      families = families,
      aliases  = aliases
    ))
  
  counts <- df_classified %>%
    filter(!is.na(Family), Family != "") %>%
    count(Family, name = "No. of genes", sort = TRUE)
  
  print(counts, n = Inf)
  write_csv(counts, out_counts)
  
  unassigned <- df_classified %>% filter(is.na(Family) | Family == "")
  if (nrow(unassigned) > 0) write_csv(unassigned, out_unassigned)
}

# ==============================================================================================
# Transcription Factors (TF)
# ==============================================================================================

aliases_TF <- c(
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

run_family_count(
  "genes_FT.csv",
  "familias_canonicas_planttfdb_like.csv",
  aliases_TF,
  "counts_TF_families.csv",
  "TF_unassigned_genes.csv"
)

# ==============================================================================================
# Epigenetic Regulators (ER)
# ==============================================================================================

aliases_ER <- c(
  "ASHH" = "SET",
  "ASHR" = "SET",
  "SDG"  = "SET",
  "SET"  = "SET",
  "ATXR" = "ATXR",
  "ATX"  = "ATX",
  "SUVH" = "SUVH",
  "SUVR" = "SUVR",
  "CHR"  = "CHR",
  "BRM"  = "BRM",
  "SYD"  = "SYD",
  "INO80"= "INO80",
  "DDM1" = "DDM1",
  "EBS"  = "EBS",
  "ROS1" = "ROS1",
  "MYB"  = "MYB"
)

run_family_count(
  "genes_RE.csv",
  "familias_epigeneticas_canonicas_v2.csv",
  aliases_ER,
  "counts_ER_families.csv",
  "ER_unassigned_genes.csv"
)