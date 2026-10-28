#!/usr/bin/env Rscript
# 06_visualise.R
# Figures and tables for the cholera genomics project.
#
# Run from the project root:   Rscript R/06_visualise.R
#
# Inputs  (written by scripts 01-05):
#   config/metadata.tsv
#   results/qc/read_summary.tsv
#   results/assembly_qc/assembly_qc_summary.tsv
#   results/amr/amr_long.tsv, results/amr/amr_matrix.tsv
#   results/phylogeny/iqtree/cholera.treefile
#   results/phylogeny/snp_dists_raw.tsv, snp_dists_masked.tsv
# Outputs: results/figures/*.png|pdf, results/tables/*.tsv, logs/R_sessionInfo.txt
#
# Wording rule for all figures: these show resistance-ASSOCIATED determinants
# (genotype). They do not show phenotypic resistance.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ape)
  library(ggtree)
  library(ggnewscale)
  library(pheatmap)
})

# ---- tunable layout settings --------------------------------------------------
LABEL_SPACE <- 0.45   # space (fraction of tree width) reserved for tip labels
STRIP_W     <- 0.04   # width of each metadata strip (fraction of tree width)
AMR_W_PER   <- 0.035  # width per AMR column (fraction of tree width)
GAP         <- 0.01

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("logs",            recursive = TRUE, showWarnings = FALSE)

save_fig <- function(plot, name, w = 8, h = 6) {
  ggsave(file.path("results/figures", paste0(name, ".png")), plot, width = w, height = h, dpi = 300)
  ggsave(file.path("results/figures", paste0(name, ".pdf")), plot, width = w, height = h)
}
have <- function(path) { ok <- file.exists(path); if (!ok) message("SKIP (missing): ", path); ok }

# ---- metadata -----------------------------------------------------------------
meta <- read_tsv("config/metadata.tsv", show_col_types = FALSE,
                 na = c("", "NA", "TODO", "unknown")) |>
  distinct(isolate_id, .keep_all = TRUE)

ref_id      <- meta |> filter(role == "reference") |> pull(isolate_id) |> first()
outgroup_id <- meta |> filter(role == "outgroup")  |> pull(isolate_id) |> first()

# =============================================================================
# 1. Read QC
# =============================================================================
if (have("results/qc/read_summary.tsv")) {
  rs <- read_tsv("results/qc/read_summary.tsv", show_col_types = FALSE)
  p_reads <- rs |>
    pivot_longer(c(pairs_before, pairs_after), names_to = "stage", values_to = "pairs") |>
    mutate(stage = recode(stage, pairs_before = "Before trimming", pairs_after = "After trimming"),
           stage = factor(stage, levels = c("Before trimming", "After trimming"))) |>
    ggplot(aes(reorder(isolate, pairs), pairs / 1e6, fill = stage)) +
    geom_col(position = "dodge") + coord_flip() +
    scale_fill_manual(values = c("grey70", "steelblue")) +
    labs(x = NULL, y = "Read pairs (millions)", fill = NULL,
         title = "Sequencing yield per isolate") +
    theme_minimal(base_size = 11)
  save_fig(p_reads, "fig_qc_reads", 7, max(3, 0.35 * nrow(rs) + 2))

  p_q30 <- ggplot(rs, aes(reorder(isolate, q30_after_pct), q30_after_pct)) +
    geom_col(fill = "steelblue") + coord_flip() +
    labs(x = NULL, y = "Q30 bases after trimming (%)", title = "Read quality per isolate") +
    theme_minimal(base_size = 11)
  save_fig(p_q30, "fig_qc_q30", 7, max(3, 0.35 * nrow(rs) + 2))
}

# =============================================================================
# 2. Assembly QC
# =============================================================================
if (have("results/assembly_qc/assembly_qc_summary.tsv")) {
  aq <- read_tsv("results/assembly_qc/assembly_qc_summary.tsv", show_col_types = FALSE,
                 na = c("", "NA", "None"))
  p_asm <- aq |>
    select(isolate, status, contigs, N50, completeness, contamination, ANI_vs_ref) |>
    pivot_longer(-c(isolate, status), names_to = "metric", values_to = "value") |>
    filter(!is.na(value)) |>
    ggplot(aes(reorder(isolate, value), value, fill = status)) +
    geom_col() + coord_flip() +
    facet_wrap(~metric, scales = "free", nrow = 1) +
    scale_fill_manual(values = c(PASS = "steelblue", FLAG = "firebrick")) +
    labs(x = NULL, y = NULL, fill = "QC status",
         title = "Assembly quality metrics",
         subtitle = "FLAG = outside the thresholds set in 03_assembly.sh") +
    theme_minimal(base_size = 10)
  save_fig(p_asm, "fig_assembly_qc", 12, max(3.5, 0.35 * nrow(aq) + 2))
}

# =============================================================================
# 3. SNP distances (before vs after recombination masking)
# =============================================================================
read_dist <- function(path) {
  m <- as.matrix(read.delim(path, row.names = 1, check.names = FALSE))
  rn <- rownames(m); rn[rn == "Reference"] <- ref_id
  rownames(m) <- rn; colnames(m) <- rn
  m
}
if (have("results/phylogeny/snp_dists_masked.tsv")) {
  dm <- read_dist("results/phylogeny/snp_dists_masked.tsv")
  pheatmap(dm, display_numbers = TRUE, number_format = "%.0f",
           main = "Pairwise core SNP distances (recombination masked)",
           filename = "results/figures/fig_snp_distance_heatmap.png",
           width = 8, height = 7)
  pheatmap(dm, display_numbers = TRUE, number_format = "%.0f",
           main = "Pairwise core SNP distances (recombination masked)",
           filename = "results/figures/fig_snp_distance_heatmap.pdf",
           width = 8, height = 7)
  write_tsv(as_tibble(dm, rownames = "isolate"), "results/tables/snp_distances_masked.tsv")

  if (have("results/phylogeny/snp_dists_raw.tsv")) {
    dr <- read_dist("results/phylogeny/snp_dists_raw.tsv")
    common <- intersect(rownames(dr), rownames(dm))
    idx <- which(upper.tri(dm[common, common]), arr.ind = TRUE)
    cmp <- tibble(pair = paste(common[idx[, 1]], common[idx[, 2]], sep = " vs "),
                  raw = dr[common, common][idx], masked = dm[common, common][idx])
    p_cmp <- ggplot(cmp, aes(raw, masked)) +
      geom_abline(linetype = "dashed", colour = "grey50") +
      geom_point(alpha = 0.7, colour = "steelblue") +
      labs(x = "SNP distance before masking", y = "SNP distance after masking",
           title = "Effect of recombination masking on pairwise distances",
           subtitle = "Points below the dashed line lost SNPs to masked recombinant regions") +
      theme_minimal(base_size = 11)
    save_fig(p_cmp, "fig_recombination_effect", 6, 5)
  }
}

# =============================================================================
# 4. AMR summaries
# =============================================================================
if (have("results/amr/amr_long.tsv")) {
  al <- read_tsv("results/amr/amr_long.tsv", show_col_types = FALSE, na = c("", "NA"))

  af <- al |> filter(tool == "AMRFinderPlus")

  # determinants by drug class (AMRFinderPlus only, to avoid double counting)
  amr_hits <- af |> filter(is.na(type) | type == "AMR")
  if (nrow(amr_hits) > 0) {
    p_class <- amr_hits |>
      mutate(class = replace_na(class, "unspecified")) |>
      count(isolate, class) |>
      ggplot(aes(isolate, n, fill = class)) +
      geom_col() + coord_flip() +
      scale_fill_viridis_d(option = "turbo") +
      labs(x = NULL, y = "Resistance-associated determinants detected", fill = "Drug class",
           title = "Resistance-associated determinants by class (AMRFinderPlus)",
           subtitle = "Genotype only: presence of a determinant does not prove phenotypic resistance") +
      theme_minimal(base_size = 11)
    save_fig(p_class, "fig_amr_by_class", 9, max(3.5, 0.35 * n_distinct(amr_hits$isolate) + 2))
  }

  # virulence / stress elements reported by --plus, kept separate from AMR
  vir <- af |> filter(!is.na(type), type != "AMR")
  if (nrow(vir) > 0) write_tsv(vir, "results/tables/amrfinder_non_amr_elements.tsv")

  write_tsv(al, "results/tables/amr_all_hits.tsv")
}

# =============================================================================
# 5. Phylogeny with metadata strips and AMR heatmap
# =============================================================================
tree_path <- "results/phylogeny/iqtree/cholera.treefile"
if (have(tree_path)) {
  tree <- read.tree(tree_path)
  tree$tip.label[tree$tip.label == "Reference"] <- ref_id

  rooted <- FALSE
  if (!is.na(outgroup_id) && outgroup_id %in% tree$tip.label) {
    tree <- root(tree, outgroup = outgroup_id, resolve.root = TRUE)
    rooted <- TRUE
  }
  tree <- ladderize(tree)

  tips <- tree$tip.label
  idx  <- match(tips, meta$isolate_id)
  lab  <- function(x) factor(replace_na(as.character(x), "unknown"))

  df_country <- data.frame(Country = lab(meta$country[idx]),            row.names = tips)
  df_year    <- data.frame(Year    = lab(meta$collection_year[idx]),    row.names = tips)

  p <- ggtree(tree) + geom_tiplab(size = 3, offset = 0) + geom_treescale(fontsize = 3)

  # IQ-TREE writes node labels as "SH-aLRT/UFboot"; show UFboot at internal nodes
  if (!is.null(tree$node.label)) {
    p$data$ufboot <- suppressWarnings(as.numeric(sub(".*/", "", p$data$label)))
    p <- p + geom_text2(aes(subset = !isTip & !is.na(ufboot), label = ufboot),
                        size = 2.3, hjust = 1.15, vjust = -0.5, colour = "grey30")
  }

  xmax <- max(p$data$x)
  o1 <- xmax * LABEL_SPACE
  o2 <- o1 + xmax * (STRIP_W + GAP)
  o3 <- o2 + xmax * (STRIP_W + GAP)

  g <- gheatmap(p, df_country, offset = o1, width = STRIP_W, colnames_angle = 90,
                colnames_position = "bottom", font.size = 2.8, hjust = 1) +
    scale_fill_viridis_d(name = "Country", option = "D")

  g <- g + new_scale_fill()
  g <- gheatmap(g, df_year, offset = o2, width = STRIP_W, colnames_angle = 90,
                colnames_position = "bottom", font.size = 2.8, hjust = 1) +
    scale_fill_viridis_d(name = "Year", option = "C")

  n_amr <- 0
  if (have("results/amr/amr_matrix.tsv")) {
    aw <- read_tsv("results/amr/amr_matrix.tsv", col_types = cols(.default = "c"))
    if (ncol(aw) > 1) {
      amr_df <- aw |> column_to_rownames("isolate")
      # most frequent determinants first
      amr_df <- amr_df[, order(-colSums(amr_df != "-", na.rm = TRUE)), drop = FALSE]
      amr_tip <- amr_df[match(tips, rownames(amr_df)), , drop = FALSE]
      rownames(amr_tip) <- tips
      amr_tip[is.na(amr_tip)] <- "not assessed"
      lv <- c("AB+AF", "AF", "AB", "-", "not assessed")
      amr_tip[] <- lapply(amr_tip, function(x) factor(x, levels = lv))
      n_amr <- ncol(amr_tip)

      g <- g + new_scale_fill()
      g <- gheatmap(g, amr_tip, offset = o3, width = AMR_W_PER * n_amr,
                    colnames_angle = 90, colnames_position = "bottom",
                    font.size = 2.3, hjust = 1) +
        scale_fill_manual(name = "Detected by", drop = FALSE,
                          values = c("AB+AF" = "#b2182b", "AF" = "#ef8a62", "AB" = "#fddbc7",
                                     "-" = "#f7f7f7", "not assessed" = "#bdbdbd"),
                          labels = c("AB+AF" = "Both tools", "AF" = "AMRFinderPlus only",
                                     "AB" = "ABRicate only", "-" = "Not detected",
                                     "not assessed" = "Not assessed"))
    }
  }

  g <- g + coord_cartesian(clip = "off") +
    labs(title = "Core-genome phylogeny with resistance-associated determinants",
         subtitle = paste0(if (rooted) paste0("Rooted on ", outgroup_id, ". ") else "Unrooted (no outgroup in tree). ",
                           "Node labels: ultrafast bootstrap. Genotype only; not phenotype.")) +
    theme(plot.margin = margin(5, 5, 70, 5), legend.position = "right")

  save_fig(g, "fig_tree_amr", 13, max(5, 0.4 * length(tips) + 3))
}

writeLines(capture.output(sessionInfo()), "logs/R_sessionInfo.txt")
message("Done. Figures in results/figures/, tables in results/tables/.")
message("Check the tree figure by eye: adjust LABEL_SPACE / STRIP_W at the top if labels overlap strips.")
