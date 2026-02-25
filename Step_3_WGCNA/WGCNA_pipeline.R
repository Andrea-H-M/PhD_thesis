##################################################################################################
## Topic: WGCNA pipeline (DESeq2 VST + variance filtering + module detection + functional cores)##
## Author: Olga Andrea Hernandez Miranda (Miranda H)                                            ##
## Date: 02/09/2025                                                                             ##                                                                                              ##
## Description:                                                                                 ##
## This script performs:                                                                        ##
##  1) Read count matrix, QC genes (goodSamplesGenes)                                           ##
##  2) Optional PCA (raw counts, exploratory)                                                   ##
##  3) DESeq2 filtering + VST normalization                                                     ##
##  4) Variance-based gene ranking + elbow (inflection)                                         ##
##  5) Final gene set = top-N variable genes + keep-list IDs                                    ##
##  6) WGCNA: soft threshold, TOM, dynamic modules, module merging                              ##
##  7) Export: Cytoscape edges/nodes, geneInfo (MM/p.MM), kWithin                               ##
##  8) Functional cores: per module (MM threshold + kWithin quantile)                           ##
##  9) Core biplots per module (MM vs kWithin)                                                  ##
##################################################################################################

allowWGCNAThreads()
cor <- stats::cor

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(readr)
  library(tidyr)
  library(WGCNA)
  library(inflection)
  library(scales)
  library(viridis)
  library(igraph)
})

# ==============================================================================================
# [1] CONFIGURATION
# ==============================================================================================

setwd("C:/Users/andii/OneDrive/Documents/DoctoradoEnCiencias/Proyecto/Tutorial7/RedConEtiquetas")

# --- Input files ---
counts_file          <- "Red.csv"               # counts matrix (genes x samples), comma-separated, rownames = gene IDs
keep_ids_file        <- "lista_sin_duplicados.txt"
annotated_counts_file<- "Red_Ann.csv"           # counts filtered + annotated (genes x samples), rownames = annotated gene IDs

# --- Sample metadata (EDIT if needed) ---
sampleCondition <- factor(c("G1","G1","G1","G1","G6","G6","G6","G6"))
sampleTable <- data.frame(condition = sampleCondition)

# --- Filtering parameters ---
min_count      <- 15     # DESeq2 prefilter: counts >= min_count
min_samples    <- 7      # in at least min_samples samples
top_n_var      <- 4637   # number of top-variance genes to keep
softPower      <- 16     # WGCNA soft power (adjust after pickSoftThreshold)
minModuleSize  <- 50
mergeCutHeight <- 0.40   # mergeCloseModules cutHeight

# --- Outlier detection ---
cutHeight_outliers <- 235

# --- Cytoscape export ---
cyto_threshold <- 0.40

# --- Functional core thresholds ---
mm_threshold     <- 0.80
kwithin_quantile <- 0.90

# --- Outputs ---
out_var_table         <- "VarianzaGenes_extendido.csv"
out_geneInfo          <- "geneInfo.csv"
out_kWithin           <- "ConectividadIntromodular.csv"
out_cores             <- "NucleosFuncionales.csv"
out_biplot_dir        <- "biplots_nucleo"
out_cyto_edge         <- "CytoscapeEdgeFile.txt"
out_cyto_node         <- "CytoscapeNodeFile.txt"
out_modTOM_rds        <- "modTOM.rds"
out_moduleColors_rds  <- "moduleColors.rds"

# ==============================================================================================
# [2] HELPERS
# ==============================================================================================

stop_if_missing <- function(files) {
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) stop("Missing file(s):\n - ", paste(missing, collapse = "\n - "))
}

make_safe_name <- function(x) str_replace_all(x, "[^[:alnum:]_\\-]", "_")

mm_col_for_module <- function(tbl, mod) {
  # Try common patterns (MM.mod, kME.mod, etc.)
  cand <- c(sprintf("MM.%s", mod),
            sprintf("MM%s",  mod),
            sprintf("kME.%s", mod),
            sprintf("kME%s",  mod))
  hits <- cand[cand %in% names(tbl)]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

# ==============================================================================================
# [3] INPUT CHECKS
# ==============================================================================================

stop_if_missing(c(counts_file, keep_ids_file, annotated_counts_file))

# ==============================================================================================
# [4] READ COUNTS + BASIC QC
# ==============================================================================================

counts_raw <- read.csv(counts_file, sep = ",", row.names = 1, check.names = FALSE)
gsg0 <- goodSamplesGenes(t(counts_raw))
counts_raw <- counts_raw[gsg0$goodGenes == TRUE, , drop = FALSE]

# Keep-list IDs
keep_ids <- readLines(keep_ids_file)
keep_ids <- keep_ids[keep_ids != ""]

# ==============================================================================================
# [5] OPTIONAL PCA (RAW COUNTS)
# ==============================================================================================

pca <- prcomp(t(counts_raw))
pca.dat <- as.data.frame(pca$x)
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var / sum(pca.var) * 100, 2)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(aes(label = rownames(pca.dat)), vjust = -0.6, size = 3) +
  labs(
    x = paste0("PC1: ", pca.var.percent[1], " %"),
    y = paste0("PC2: ", pca.var.percent[2], " %")
  ) +
  theme_minimal()

# ==============================================================================================
# [6] DESeq2 VST NORMALIZATION
# ==============================================================================================

dds <- DESeqDataSetFromMatrix(countData = counts_raw, colData = sampleTable, design = ~ condition)

dds_filt <- dds[rowSums(counts(dds) >= min_count) >= min_samples, ]
dds_norm <- varianceStabilizingTransformation(dds_filt, blind = FALSE)

expressiondata <- assay(dds_norm) %>% t()
expression <- as.data.frame(expressiondata)

gsg <- goodSamplesGenes(expression, verbose = 3)
if (!gsg$allOK) {
  expression <- expression[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
}

# ==============================================================================================
# [7] VARIANCE RANKING + ELBOW (INFLECTION)
# ==============================================================================================

gene_vars  <- apply(expression, 2, var)
gene_means <- apply(expression, 2, mean)
gene_sd    <- apply(expression, 2, sd)
coef_var   <- gene_sd / gene_means

var_df <- data.frame(
  Gen      = names(gene_vars),
  Varianza = gene_vars,
  Media    = gene_means,
  SD       = gene_sd,
  CoefVar  = coef_var
) %>%
  arrange(desc(Varianza)) %>%
  mutate(Rank = row_number())

write.csv(var_df, out_var_table, row.names = FALSE)

var_curve <- var_df %>%
  transmute(
    Gen = Gen,
    Varianza = Varianza
  ) %>%
  mutate(
    Varianza_acumulada = cumsum(Varianza) / sum(Varianza),
    Rango = row_number()
  )

# Elbow point (inflection)
el <- inflection::uik(var_curve$Rango, var_curve$Varianza_acumulada)
elbow_idx <- if (length(el) == 1) as.integer(el) else as.integer(el[2])
elbow_idx <- max(1, min(elbow_idx, nrow(var_curve)))
elbow_row <- var_curve[elbow_idx, , drop = FALSE]

cat("Turning point (Rank): ", elbow_row$Rango, "\n", sep = "")
cat("Accumulated variance at turning point: ", elbow_row$Varianza_acumulada, "\n", sep = "")

ggplot(var_curve, aes(x = Rango, y = Varianza_acumulada)) +
  geom_line() +
  geom_point(data = elbow_row, aes(x = Rango, y = Varianza_acumulada), size = 3) +
  geom_text(
    data = elbow_row,
    aes(label = paste0("Turning point = ", Rango, " genes\n",
                       percent(Varianza_acumulada, accuracy = 0.1))),
    vjust = -0.8
  ) +
  labs(x = "Number of genes", y = "Fraction of accumulated variance") +
  theme_minimal()

# ==============================================================================================
# [8] FINAL GENE SET: TOP-N VARIANCE + KEEP LIST
# ==============================================================================================

top_n <- min(top_n_var, nrow(var_df))
genes_top_var <- var_df$Gen[1:top_n]
genes_finales <- unique(c(genes_top_var, keep_ids))

expression_filtrada <- expression[, colnames(expression) %in% genes_finales, drop = FALSE]

# ==============================================================================================
# [9] LOAD ANNOTATED (POST-PERL) MATRIX + OUTLIER REMOVAL
# ==============================================================================================

expression_filtrada_ann <- read.csv(annotated_counts_file, sep = ",", row.names = 1, check.names = FALSE)

sampleTree <- hclust(dist(expression_filtrada_ann), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers")
abline(h = cutHeight_outliers, col = "red")

clust <- cutreeStatic(sampleTree, cutHeight = cutHeight_outliers, minSize = 1)
expression_wgcna <- expression_filtrada_ann[clust == 1, , drop = FALSE]

# ==============================================================================================
# [10] WGCNA: SOFT POWER + TOM + MODULES + MERGE
# ==============================================================================================

powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(expression_wgcna, powerVector = powers, verbose = 5)

# Adjacency + TOM
adjacency_mat <- adjacency(expression_wgcna, power = softPower)
TOM <- TOMsimilarity(adjacency_mat)
dissTOM <- 1 - TOM

geneTree <- hclust(as.dist(dissTOM), method = "average")
plot(geneTree, main = "Gene clustering (TOM dissimilarity)", labels = FALSE)

dynamicMods <- cutreeDynamic(
  dendro = geneTree, distM = dissTOM,
  deepSplit = 1, pamRespectsDendro = FALSE,
  minClusterSize = minModuleSize
)
dynamicColors <- labels2colors(dynamicMods)

plotDendroAndColors(geneTree, dynamicColors, "Initial modules", dendroLabels = FALSE)

# Merge close modules
MEList0 <- moduleEigengenes(expression_wgcna, colors = dynamicColors)
MEs0 <- MEList0$eigengenes
MEDiss <- 1 - cor(MEs0)
METree <- hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Eigengene clustering", labels = FALSE)
abline(h = mergeCutHeight, col = "red", lty = 2)

merge <- mergeCloseModules(expression_wgcna, dynamicColors, cutHeight = mergeCutHeight)
moduleColors <- merge$colors
MEs <- merge$newMEs

plotDendroAndColors(
  geneTree,
  cbind(dynamicColors, moduleColors),
  c("Initial", "Merged"),
  dendroLabels = FALSE
)

cat("Final modules (cutHeight = ", mergeCutHeight, "): ",
    length(unique(moduleColors)), "\n", sep = "")

# Save objects for downstream network plots
saveRDS(moduleColors, out_moduleColors_rds)
# Save a TOM restricted to non-grey genes for Cytoscape convenience
genes_all <- colnames(expression_wgcna)
inModule <- moduleColors != "grey"
modGenes <- genes_all[inModule]
modTOM <- TOM[inModule, inModule]
dimnames(modTOM) <- list(modGenes, modGenes)
saveRDS(modTOM, out_modTOM_rds)

# ==============================================================================================
# [11] EXPORT NETWORK TO CYTOSCAPE (NON-GREY)
# ==============================================================================================

cyt <- exportNetworkToCytoscape(
  modTOM,
  edgeFile   = out_cyto_edge,
  nodeFile   = out_cyto_node,
  weighted   = TRUE,
  threshold  = cyto_threshold,
  nodeNames  = modGenes,
  nodeAttr   = moduleColors[inModule]
)

# ==============================================================================================
# [12] GENE INFO: MODULE MEMBERSHIP (MM) + p.MM
# ==============================================================================================

modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(expression_wgcna, MEs, use = "p"))
nSamples <- nrow(expression_wgcna)
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue)             <- paste("p.MM", modNames, sep = "")

geneInfo0 <- data.frame(Gene = colnames(expression_wgcna), moduleColor = moduleColors)

for (i in seq_len(ncol(geneModuleMembership))) {
  geneInfo0 <- cbind(geneInfo0,
                     geneModuleMembership[, i, drop = FALSE],
                     MMPvalue[, i, drop = FALSE])
  colnames(geneInfo0)[(ncol(geneInfo0)-1):ncol(geneInfo0)] <-
    c(paste0("MM.", modNames[i]), paste0("p.MM.", modNames[i]))
}

geneInfo <- geneInfo0[order(geneInfo0$moduleColor), ]
write.csv(geneInfo, file = out_geneInfo, row.names = FALSE)

# ==============================================================================================
# [13] INTRAMODULAR CONNECTIVITY (kWithin)
# ==============================================================================================

ADJ <- adjacency(expression_wgcna, power = softPower)
intraModConn <- intramodularConnectivity(ADJ, moduleColors)
write.csv(intraModConn, file = out_kWithin, row.names = FALSE)

# ==============================================================================================
# [14] FUNCTIONAL CORES (MM >= threshold AND kWithin >= quantile)
# ==============================================================================================

kWithin_df <- read_csv(out_kWithin, show_col_types = FALSE)
MM_df      <- read_csv(out_geneInfo, show_col_types = FALSE)

# Ensure key column naming
if (!("Gene" %in% names(kWithin_df))) names(kWithin_df)[1] <- "Gene"
if (!("Gene" %in% names(MM_df)))      names(MM_df)[1]      <- "Gene"

df_core_base <- inner_join(kWithin_df, MM_df, by = "Gene")

mods <- unique(df_core_base$moduleColor)
mods <- mods[mods != "grey"]

cores_list <- list()

for (mod in mods) {
  sub <- df_core_base %>% filter(moduleColor == mod)
  
  mm_col <- mm_col_for_module(sub, mod)
  if (is.na(mm_col)) {
    message("MM column not found for module ", mod, ". Skipping.")
    next
  }
  
  k_cut <- quantile(sub$kWithin, kwithin_quantile, na.rm = TRUE)
  
  core_mod <- sub %>%
    mutate(MM = .data[[mm_col]]) %>%
    filter(MM >= mm_threshold, kWithin >= k_cut)
  
  if (nrow(core_mod) > 0) cores_list[[mod]] <- core_mod
}

cores_df <- bind_rows(cores_list)

# Keep only minimal columns if you prefer (uncomment)
# cores_df <- cores_df %>% select(Gene, moduleColor, kWithin, starts_with("MM."))

write_csv(cores_df, out_cores)

# ==============================================================================================
# [15] BIPLOTS PER MODULE (MM vs kWithin)
# ==============================================================================================

if (!dir.exists(out_biplot_dir)) dir.create(out_biplot_dir, recursive = TRUE)

summary_list <- list()

for (mod in mods) {
  sub <- df_core_base %>% filter(moduleColor == mod)
  
  mm_col <- mm_col_for_module(sub, mod)
  if (is.na(mm_col) || nrow(sub) == 0) next
  
  k_cut <- quantile(sub$kWithin, kwithin_quantile, na.rm = TRUE)
  
  plot_df <- sub %>%
    mutate(
      MM = .data[[mm_col]],
      Core = ifelse(MM >= mm_threshold & kWithin >= k_cut, "Core", "Non-core")
    )
  
  p <- ggplot(plot_df, aes(x = MM, y = kWithin, color = Core)) +
    geom_point(alpha = 0.7) +
    labs(
      title = paste("Functional core biplot - module", mod),
      subtitle = sprintf("MM ≥ %.2f | kWithin ≥ p%.0f", mm_threshold, 100 * kwithin_quantile),
      x = paste0("Module membership (", mod, ")"),
      y = "Intramodular connectivity (kWithin)"
    ) +
    scale_color_manual(values = c("Core" = mod, "Non-core" = "gray85")) +
    theme_minimal(base_family = "Arial") +
    theme(legend.title = element_blank())
  
  safe_mod <- make_safe_name(mod)
  png_file <- file.path(out_biplot_dir, paste0("biplot_core_", safe_mod, ".png"))
  csv_file <- file.path(out_biplot_dir, paste0("core_data_", safe_mod, ".csv"))
  
  ggsave(png_file, p, width = 7, height = 5, dpi = 300)
  write_csv(plot_df, csv_file)
  
  summary_list[[mod]] <- tibble(
    module = mod,
    n_genes = nrow(plot_df),
    n_core = sum(plot_df$Core == "Core"),
    k_cutoff = as.numeric(k_cut),
    mm_threshold = mm_threshold,
    mm_column = mm_col
  )
}

if (length(summary_list) > 0) {
  write_csv(bind_rows(summary_list), file.path(out_biplot_dir, "summary_biplots_core.csv"))
}

cat("\nDone.\n",
    "- geneInfo: ", out_geneInfo, "\n",
    "- kWithin: ", out_kWithin, "\n",
    "- cores:   ", out_cores, "\n",
    "- biplots: ", out_biplot_dir, "\n", sep = "")