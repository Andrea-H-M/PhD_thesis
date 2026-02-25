##################################################################################################
## Topic: KEGG functional enrichment (KO list)                                                  ##
## Author: Olga Andrea Hernandez Miranda (Miranda H)                                            ##
## Date: 05/09/20250                                                                            ##                                                                        ##
## Description:                                                                                 ##
## This script reads a KO list (e.g., from GhostKOALA output), normalizes KO identifiers,       ##    ##
## maps KO -> KEGG pathways, runs enrichment with clusterProfiler::enricher(),                  ##
## exports a results table, and generates barplot + dotplot (top 20 pathways).                  ##
##################################################################################################

allowWGCNAThreads()

setwd("C:/Users/andii/OneDrive/Documents/DoctoradoEnCiencias/Proyecto/Tutorial7/RedConEtiquetas/GenesPorModulo_ID1/kegg/turquoise")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(KEGGREST)
  library(clusterProfiler)
  library(ggplot2)
})

# ==============================================================================================
# [1] INPUTS / PARAMETERS
# ==============================================================================================

ko_file <- "turquoise_kobas.txt"      # must contain a column named "KO"

p_cutoff  <- 0.05
q_cutoff  <- 0.05
padj_meth <- "BH"
min_gs    <- 5
max_gs    <- 5000
topN_plot <- 20

out_table <- "KEGG_turquoise_enrichment.csv"
out_bar   <- "KEGG_enrichment_turquoise_barplot_NAMES.png"
out_dot   <- "KEGG_enrichment_turquoise_dotplot_NAMES.png"

# ==============================================================================================
# [2] READ KO LIST + NORMALIZE
# ==============================================================================================

fg <- read_tsv(ko_file, show_col_types = FALSE)
stopifnot("KO" %in% names(fg))

ko_fg <- fg$KO %>%
  as.character() %>%
  str_trim() %>%
  toupper() %>%
  str_replace("^KO:", "") %>%   # KO:Kxxxxx -> Kxxxxx
  str_replace("^KO",  "K") %>%  # KOxxxxx  -> Kxxxxx
  keep(~ grepl("^K\\d{5}$", .x)) %>%
  unique()

cat("Unique valid KOs:", length(ko_fg), "\n")
if (length(ko_fg) == 0) stop("No valid KOs were detected. Check the KO column format.")

# ==============================================================================================
# [3] KEGG DICTIONARIES (KO -> pathway, pathway names)
# ==============================================================================================

message("Downloading KO -> pathway links and pathway names from KEGG...")

links <- keggLink("pathway", "ko")   # names=ko:Kxxxxx ; values=path:koXXXXX

term2gene <- tibble(
  ID = sub("^path:", "", unname(links)),  # koXXXXX
  KO = sub("^ko:",   "", names(links))    # Kxxxxx
)

plist <- keggList("pathway")

name_lut <- tibble(
  ID       = sub("^path:", "", names(plist)),                 # koXXXXX / mapXXXXX
  TermName = sub(" - Reference pathway$", "", as.character(plist))
) %>%
  filter(grepl("^(ko|map)\\d{5}$", ID)) %>%
  mutate(ID = sub("^map", "ko", ID)) %>%                      # normalize map -> ko
  distinct(ID, .keep_all = TRUE)

n_overlap <- sum(ko_fg %in% term2gene$KO)
cat("KOs mapping to at least one pathway:", n_overlap, "\n")
if (n_overlap == 0) stop("Your KOs do not intersect KEGG pathways. Check KO format or input file.")

# ==============================================================================================
# [4] ENRICHMENT (no custom background)
# ==============================================================================================

ekegg <- enricher(
  gene          = ko_fg,
  TERM2GENE     = term2gene,
  TERM2NAME     = NULL,     # avoids Description collisions; we add names later safely
  pvalueCutoff  = p_cutoff,
  qvalueCutoff  = q_cutoff,
  pAdjustMethod = padj_meth,
  minGSSize     = min_gs,
  maxGSSize     = max_gs
)

# If no significant terms, relax thresholds for diagnostics
if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
  message("No terms at q<=0.05. Relaxing thresholds for diagnostics...")
  ekegg <- enricher(
    gene          = ko_fg,
    TERM2GENE     = term2gene,
    TERM2NAME     = NULL,
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    pAdjustMethod = padj_meth,
    minGSSize     = 3
  )
}

res <- as.data.frame(ekegg)
if (nrow(res) == 0) stop("No enriched terms found even after relaxing thresholds.")

# ==============================================================================================
# [5] ADD SAFE PATHWAY NAMES + EXPORT TABLE
# ==============================================================================================

res$ID <- sub("^map", "ko", res$ID)

res2 <- res %>%
  left_join(name_lut, by = "ID") %>%
  mutate(
    Label = ifelse(!is.na(TermName), TermName, ID),
    PathwayURL = paste0("https://www.kegg.jp/pathway/", ID)
  ) %>%
  arrange(p.adjust)

write.csv(res2, out_table, row.names = FALSE)
cat("Rows in result table:", nrow(res2), "\n")

# ==============================================================================================
# [6] PLOTS (TOP N)
# ==============================================================================================

res2 <- res2 %>%
  mutate(
    GeneRatio_num =
      as.numeric(sub("/.*", "", GeneRatio)) /
      as.numeric(sub(".*/", "", GeneRatio))
  )

topN <- min(topN_plot, nrow(res2))
plot_df <- res2 %>% slice_head(n = topN)

# Barplot (Count)
p_bar <- ggplot(plot_df, aes(x = Count, y = reorder(Label, Count), fill = p.adjust)) +
  geom_col() +
  labs(x = "Number of genes", y = "KEGG pathways") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9))

ggsave(out_bar, p_bar, width = 9, height = 6, dpi = 300)

# Dotplot (GeneRatio)
p_dot <- ggplot(plot_df, aes(x = GeneRatio_num, y = reorder(Label, GeneRatio_num),
                             size = Count, color = p.adjust)) +
  geom_point() +
  labs(x = "GeneRatio", y = "KEGG pathways") +
  theme_bw()

ggsave(out_dot, p_dot, width = 9, height = 6, dpi = 300)

cat("Done.\n- Table: ", out_table, "\n- Barplot: ", out_bar, "\n- Dotplot: ", out_dot, "\n", sep = "")