################################################################################################
## Topic: UpSet plot for label combinations + intersection counts export                         ##
## Author: Olga Andrea Hernandez Miranda (Miranda H)                                            ##
## Date: 22/09/2026                                                                                                                      ##
## Description:                                                                               ##
## This script reads a binary label table (0/1, TRUE/FALSE, yes/no, etc.),                      ##
## builds an UpSet combination matrix (ComplexHeatmap), plots intersections ordered by size,   ##
## and exports a CSV with intersection names and counts.                                       ##
################################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(viridisLite)
})

# ==============================================================================================
# [1] BASIC CONFIG
# ==============================================================================================

setwd("C:/Users/andii/OneDrive/Documents/DoctoradoEnCiencias/Proyecto/Tutorial7/RedConEtiquetas/RedConID/Tablas_por_moduloo_1")

infile <- "Resumen_maestro_2oMasEtiquetas_binario.csv"

# Desired set order (only those present in the file will be used)
sets_pref <- c("GNM","GED","FT","OCDT","Fer","RE","PrePol","Pol","PostPol")

# Output (intersection counts)
out_counts <- sub("\\.csv$", "", infile, ignore.case = TRUE) %>%
  paste0("_intersections_counts.csv")

# ==============================================================================================
# [2] READ DATA + CLEAN HEADERS
# ==============================================================================================

df <- read_csv(infile, show_col_types = FALSE)
names(df) <- trimws(names(df))

sets <- intersect(sets_pref, names(df))
if (length(sets) < 2) stop("Less than 2 label columns detected. Check column names.")
stopifnot("RE" %in% sets)   # keep this if RE is mandatory; otherwise comment it out

# ==============================================================================================
# [3] COERCE TO LOGICAL (robust to 0/1, T/F, yes/no, x, si/sí)
# ==============================================================================================

to_logical <- function(x) {
  if (is.logical(x)) return(replace(x, is.na(x), FALSE))
  if (is.numeric(x)) return(replace(x == 1, is.na(x), FALSE))
  
  s <- tolower(trimws(as.character(x)))
  replace(s %in% c("1","true","t","yes","y","si","sí","x"), is.na(s), FALSE)
}

df_bin <- df %>% mutate(across(all_of(sets), to_logical))

# ==============================================================================================
# [4] BUILD COMBINATION MATRIX + ORDER BY INTERSECTION SIZE
# ==============================================================================================

m <- make_comb_mat(df_bin[, sets])

# Remove empty combinations if any
m <- m[, comb_size(m) > 0]

# Order intersections (largest -> smallest)
comb_ord <- order(comb_size(m), decreasing = TRUE)
m <- m[, comb_ord]

# ==============================================================================================
# [5] COLORS + ANNOTATIONS
# ==============================================================================================

# Set colors (edit freely)
pal <- c(
  GNM    = "#1b9e77",
  GED    = "#7570b3",
  FT     = "#d95f02",
  OCDT   = "#e6ab02",
  Fer    = "#e7298a",
  RE     = "#66a61e",
  PrePol = "#1f78b4",
  Pol    = "#a6761d",
  PostPol= "#666666"
)

# Keep only sets actually used
pal <- pal[rownames(m)]

cs <- comb_size(m)
yticks <- pretty(c(0, max(cs)), n = 6)
col_fun <- circlize::colorRamp2(seq(min(cs), max(cs), length.out = 6), viridis(6))

top_anno <- HeatmapAnnotation(
  "Intersection" = anno_barplot(
    cs, which = "column",
    gp = gpar(fill = col_fun(cs), col = NA),
    border = FALSE,
    height = unit(36, "mm"),
    bar_width = 0.65,
    ylim = c(0, max(cs)),
    axis_param = list(side = "left", at = yticks, labels = yticks)
  ),
  annotation_name_side = "left"
)

ss <- set_size(m)
lticks <- pretty(c(0, max(ss)), n = 5)

left_bars <- rowAnnotation(
  "Set size" = anno_barplot(
    -ss, which = "row",
    gp = gpar(fill = pal, col = NA),
    border = FALSE,
    width = unit(36, "mm"),
    bar_width = 0.90,
    ylim = c(-max(ss), 0),
    axis_param = list(side = "bottom", at = -lticks, labels = lticks)
  ),
  annotation_name_side = "top"
)

label_anno <- rowAnnotation(
  Labels = anno_text(
    rownames(m), which = "row",
    just = "left", location = 0.5,
    gp = gpar(fontsize = 10),
    width = unit(26, "mm")
  )
)

# ==============================================================================================
# [6] UPSET PLOT
# ==============================================================================================

ht <- UpSet(
  m,
  comb_order       = seq_len(ncol(m)),   # already ordered
  top_annotation   = top_anno,
  right_annotation = NULL,
  bg_pt_col        = "#D0D0D0",
  comb_col         = pal,
  pt_size          = unit(2.4, "mm"),
  lwd              = 1,
  show_row_names   = FALSE,
  height           = unit(length(rownames(m)) * 6, "mm")
)

lg_inter <- Legend(
  title = "Count",
  col_fun = col_fun, at = yticks, labels = yticks,
  direction = "vertical",
  legend_height = unit(45, "mm")
)

draw(
  left_bars + label_anno + ht,
  padding = unit(c(6, 6, 6, 6), "mm"),
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  annotation_legend_list = list(lg_inter),
  merge_legends = FALSE
)

# ==============================================================================================
# [7] EXPORT: INTERSECTION COUNTS TABLE
# ==============================================================================================

intersections_counts <- tibble(
  intersection = ComplexHeatmap::comb_name(m),     # e.g., "GNM&FT"
  count        = as.integer(ComplexHeatmap::comb_size(m))
) %>%
  arrange(desc(count))

print(intersections_counts, n = 20)
write_csv(intersections_counts, out_counts)

cat("\nDone.\n- Intersection counts: ", out_counts, "\n", sep = "")