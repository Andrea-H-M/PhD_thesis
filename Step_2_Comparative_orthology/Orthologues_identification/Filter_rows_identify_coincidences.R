######################################################
## Topic: Master table cleanup + matching rows       ##
## (Fill blanks with NA, then remove NA to keep      ##
## only fully-matched/complete rows)                 ##
## Author: Olga Andrea Hernandez Miranda (Miranda H) ##
## Date: 08/09/2024                                  ##
## Note: Fill empty cells with NA, then filter rows  ##
##       without NA to identify coincidences         ##
######################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# =========================
# 1) Paths & setup
# =========================

# Folder where your orthogroup results are located
directorio <- "C:/Users/andii/OneDrive/Documents/DoctoradoEnCiencias/Proyecto/Tutorial5/Results_Sep07/Orthogroups"
setwd(directorio)

# Input master table (preliminary)
input_file  <- "TablaMaestraPrelimnar.csv"

# Outputs
output_master <- "TablaMaestra.csv"
output_hits   <- "TablaMaestra_sin_NA_coincidencias.csv"

# =========================
# 2) Read master table
# =========================

tabla <- read_csv(input_file, show_col_types = FALSE)

# =========================
# 3) Fill blanks with NA (first 4 columns)
#    and standardize as character to avoid type issues
# =========================

tabla <- tabla %>%
  mutate(across(1:4, as.character)) %>%
  mutate(across(1:4, ~ na_if(., "")))  # converts "" to NA

# Save cleaned master table
write_csv(tabla, output_master)

# =========================
# 4) Identify "coincidences" by removing NA rows
#    (keeps only rows where columns 1:4 are all present)
# =========================

tabla_coincidencias <- tabla %>%
  filter(if_all(1:4, ~ !is.na(.)))   # no NA in the first 4 columns

# Save matching rows
write_csv(tabla_coincidencias, output_hits)

# =========================
# 5) Quick checks
# =========================

cat("\n--- Rows (original): ", nrow(tabla), "\n", sep = "")
cat("--- Rows (coincidences, no NA in cols 1:4): ", nrow(tabla_coincidencias), "\n", sep = "")

cat("\n--- Preview (coincidences) ---\n")
print(head(tabla_coincidencias, 10))