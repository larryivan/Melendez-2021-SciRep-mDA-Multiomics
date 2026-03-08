#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(VennDiagram)
  library(pheatmap)
  library(grid)
})

# ============================================================
# Figure 2 reproduction script (interactive / RStudio friendly)
# Single entry point for Fig2A, Fig2B, Fig2C, and the Fig2C right panel.
# Edit the parameters in this block, then Run All.
# ============================================================

count_file <- "data/GSE153005_count_mat.txt"
fpkm_file <- "data/GSE153005_FPKM_mat.txt"
outdir <- "results/fig2"

# DEG criteria from paper
padj_cutoff <- 0.05
fc_cutoff <- 4
min_fpkm <- 1

# Requested DESeq2 settings (no batch term)
fit_type <- "local"
sf_type <- "poscounts"
independent_filtering <- TRUE
cooks_cutoff <- TRUE

cluster_method <- "ward.D2"
set.seed(20260303)
export_fig2c_right_panel <- TRUE
fig2c_right_terms_per_cluster <- c(`1` = 3L, `2` = 2L, `3` = 2L, `4` = 4L, `5` = 1L, `6` = 7L)
fig2c_include_reactome <- FALSE

if (!file.exists(count_file) || !file.exists(fpkm_file)) {
  stop(
    "Input matrices not found. Please place files at: ",
    count_file, " and ", fpkm_file,
    " (or edit count_file / fpkm_file at top of script)."
  )
}

# ============================================================
# Helper functions
# ============================================================

sanitize_names <- function(x) {
  x <- gsub("\r", "", x, fixed = TRUE)
  trimws(x)
}

read_expr_matrix <- function(path, matrix_name = "count") {
  message("[I] Reading ", matrix_name, " matrix: ", path)
  df <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  names(df) <- sanitize_names(names(df))

  required_cols <- c("geneId", "geneName")
  miss <- setdiff(required_cols, names(df))
  if (length(miss) > 0) {
    stop(matrix_name, " matrix missing required columns: ", paste(miss, collapse = ", "))
  }

  sample_cols <- grep("^d(0|14|28)_r\\d+$", names(df), value = TRUE, ignore.case = TRUE)
  if (length(sample_cols) == 0) {
    stop(matrix_name, " matrix has no sample columns matching d0_r1 / d14_r1 / d28_r1.")
  }
  sample_cols <- sample_cols[order(tolower(sample_cols))]

  for (sc in sample_cols) {
    df[[sc]] <- suppressWarnings(as.numeric(df[[sc]]))
  }

  df <- df %>%
    mutate(
      geneId = as.character(geneId),
      gene_id = sub("\\..*$", "", geneId),
      geneName = as.character(geneName)
    )

  list(df = df, sample_cols = sample_cols)
}

build_sample_info <- function(sample_cols) {
  sample_info <- tibble(sample = sample_cols) %>%
    mutate(
      sample_lower = tolower(sample),
      day_raw = str_extract(sample_lower, "^d\\d+"),
      condition = case_when(
        day_raw == "d0" ~ "D0",
        day_raw == "d14" ~ "D14",
        day_raw == "d28" ~ "D28",
        TRUE ~ NA_character_
      ),
      replicate = as.integer(str_extract(sample_lower, "(?<=_r)\\d+"))
    ) %>%
    dplyr::select(-sample_lower) %>%
    filter(!is.na(condition)) %>%
    mutate(condition = factor(condition, levels = c("D0", "D14", "D28"))) %>%
    arrange(condition, replicate)

  if (nrow(sample_info) != length(sample_cols)) {
    stop("Some sample names could not be parsed into D0/D14/D28 conditions.")
  }
  sample_info
}

classify_biotype <- function(raw_biotype) {
  case_when(
    raw_biotype == "protein_coding" ~ "Protein-coding",
    str_detect(
      raw_biotype,
      regex(
        "lncrna|linc|antisense|processed_transcript|sense_intronic|sense_overlapping|3prime_overlapping_ncrna|bidirectional_promoter_lncrna|macro_lncRNA|non_coding",
        ignore_case = TRUE
      )
    ) ~ "lncRNA",
    TRUE ~ "Other"
  )
}

build_biotype_map <- function(count_df) {
  if (!"type" %in% names(count_df)) {
    stop("Count matrix is missing 'type' column, cannot classify protein-coding vs lncRNA.")
  }

  count_df %>%
    transmute(
      geneId,
      gene_id,
      geneName,
      raw_biotype = as.character(type),
      biotype = classify_biotype(raw_biotype)
    ) %>%
    distinct(geneId, .keep_all = TRUE)
}

run_pair_deseq <- function(label, group_a, group_b, counts_df, fpkm_df, sample_info, biotype_map) {
  message("[I] Running DESeq2: ", label)

  pair_samples <- sample_info %>%
    filter(condition %in% c(group_a, group_b)) %>%
    mutate(condition = factor(as.character(condition), levels = c(group_a, group_b))) %>%
    arrange(condition, replicate)

  count_mat <- counts_df %>%
    dplyr::select(geneId, all_of(pair_samples$sample)) %>%
    as.data.frame()
  rownames(count_mat) <- count_mat$geneId
  count_mat <- as.matrix(count_mat[, pair_samples$sample, drop = FALSE])
  storage.mode(count_mat) <- "integer"

  coldata <- data.frame(
    condition = pair_samples$condition,
    row.names = pair_samples$sample,
    stringsAsFactors = TRUE
  )

  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = coldata,
    design = ~ condition
  )

  dds <- DESeq(
    dds,
    fitType = fit_type,
    sfType = sf_type,
    quiet = TRUE
  )

  res <- results(
    dds,
    contrast = c("condition", group_b, group_a),
    alpha = padj_cutoff,
    independentFiltering = independent_filtering,
    cooksCutoff = cooks_cutoff
  )

  res_df <- as.data.frame(res) %>%
    rownames_to_column("geneId") %>%
    mutate(
      gene_id = sub("\\..*$", "", geneId),
      comparison = label
    )

  norm_counts <- counts(dds, normalized = TRUE)
  group_a_samples <- pair_samples %>% filter(condition == group_a) %>% pull(sample)
  group_b_samples <- pair_samples %>% filter(condition == group_b) %>% pull(sample)
  mean_a <- rowMeans(norm_counts[, group_a_samples, drop = FALSE])
  mean_b <- rowMeans(norm_counts[, group_b_samples, drop = FALSE])

  fc_tbl <- tibble(
    geneId = rownames(norm_counts),
    mean_norm_a = as.numeric(mean_a),
    mean_norm_b = as.numeric(mean_b)
  ) %>%
    mutate(
      fc_norm = if_else(mean_norm_a == 0 & mean_norm_b == 0, 1, mean_norm_b / mean_norm_a),
      abs_fc_norm = pmax(fc_norm, 1 / fc_norm),
      direction = case_when(
        mean_norm_b > mean_norm_a ~ "Up",
        mean_norm_b < mean_norm_a ~ "Down",
        TRUE ~ "NS"
      )
    )

  pair_fpkm <- fpkm_df %>%
    dplyr::select(geneId, gene_id, all_of(pair_samples$sample))

  pair_fpkm$max_fpkm_pair <- apply(
    as.matrix(pair_fpkm[, pair_samples$sample, drop = FALSE]),
    1,
    function(x) max(as.numeric(x), na.rm = TRUE)
  )
  pair_fpkm <- pair_fpkm %>% dplyr::select(geneId, gene_id, max_fpkm_pair)

  full_tbl <- res_df %>%
    left_join(fc_tbl, by = "geneId") %>%
    left_join(pair_fpkm, by = c("geneId", "gene_id")) %>%
    left_join(biotype_map, by = c("geneId", "gene_id")) %>%
    mutate(
      fc_for_filter = abs_fc_norm,
      lfc_equiv_for_filter = log2(fc_for_filter),
      is_deg = !is.na(padj) &
        padj < padj_cutoff &
        !is.na(fc_for_filter) &
        fc_for_filter > fc_cutoff &
        !is.na(max_fpkm_pair) &
        max_fpkm_pair >= min_fpkm
    )

  deg_tbl <- full_tbl %>%
    filter(is_deg) %>%
    filter(biotype %in% c("Protein-coding", "lncRNA")) %>%
    arrange(padj, desc(abs(log2FoldChange)))

  list(full = full_tbl, deg = deg_tbl)
}

save_fig2a_venn <- function(gene_sets, out_png, out_pdf) {
  message("[I] Plotting Fig2A")
  venn_grob <- venn.diagram(
    x = gene_sets,
    filename = NULL,
    imagetype = "png",
    fill = c("#43C8F2", "#C2D0F2", "#6E8FD0"),
    alpha = c(0.70, 0.55, 0.55),
    lwd = 0,
    cex = 1.7,
    cat.cex = 1.35,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-20, 20, 180),
    cat.dist = c(0.06, 0.06, 0.08),
    cat.col = c("#26AEDB", "#6E8FD0", "#355DA9"),
    margin = 0.05
  )

  png(out_png, width = 2200, height = 1800, res = 300)
  grid.newpage()
  grid.draw(venn_grob)
  dev.off()

  pdf(out_pdf, width = 8.5, height = 7.0)
  grid.newpage()
  grid.draw(venn_grob)
  dev.off()
}

save_fig2b_barplot <- function(deg_list, out_png, out_pdf, out_tsv) {
  message("[I] Plotting Fig2B")

  comparisons <- c("D0 vs D14", "D14 vs D28", "D0 vs D28")
  biotypes <- c("Protein-coding", "lncRNA")
  directions <- c("Down", "Up")

  counts <- bind_rows(lapply(names(deg_list), function(lbl) {
    deg_list[[lbl]] %>%
      dplyr::count(comparison = lbl, biotype, direction, name = "n")
  })) %>%
    complete(
      comparison = comparisons,
      biotype = biotypes,
      direction = directions,
      fill = list(n = 0)
    ) %>%
    mutate(
      comparison = factor(comparison, levels = comparisons),
      biotype = factor(biotype, levels = biotypes),
      direction = factor(direction, levels = directions),
      value = if_else(direction == "Down", -n, n),
      label_x = value + if_else(value > 0, 28, -28),
      label_hjust = if_else(value > 0, 0, 1)
    )

  write_tsv(counts, out_tsv)

  color_map <- c("Protein-coding" = "#183DFF", "lncRNA" = "#FF1A1A")
  x_limit <- max(abs(counts$value), 1)
  x_limit <- ceiling((x_limit + 80) / 100) * 100
  label_offset <- max(42, round(x_limit * 0.04))
  title_y <- length(comparisons) + 0.78
  counts <- counts %>%
    mutate(
      label_x = value + if_else(value > 0, label_offset, -label_offset),
      label_hjust = if_else(value > 0, 0, 1)
    )
  dodge <- position_dodge(width = 0.72)

  p <- ggplot(counts, aes(x = value, y = comparison, fill = biotype)) +
    geom_vline(xintercept = 0, linewidth = 0.9, color = "black") +
    geom_hline(yintercept = c(1.5, 2.5), linewidth = 0.75, color = "black") +
    geom_col(position = dodge, width = 0.58, color = NA) +
    geom_text(
      aes(x = label_x, label = n, color = biotype, hjust = label_hjust),
      position = dodge,
      size = 4.7,
      fontface = "bold",
      show.legend = FALSE
    ) +
    annotate(
      "text",
      x = -x_limit * 0.55,
      y = title_y,
      label = "Down-regulated",
      fontface = "bold",
      size = 5.0,
      vjust = 0
    ) +
    annotate(
      "text",
      x = x_limit * 0.55,
      y = title_y,
      label = "Up-regulated",
      fontface = "bold",
      size = 5.0,
      vjust = 0
    ) +
    scale_fill_manual(values = color_map) +
    scale_color_manual(values = color_map) +
    scale_x_continuous(
      limits = c(-x_limit, x_limit),
      breaks = pretty(c(-x_limit, x_limit), n = 7),
      labels = function(x) abs(x)
    ) +
    scale_y_discrete(expand = expansion(add = c(0.35, 0.95))) +
    labs(x = "Number of genes", y = NULL, fill = NULL) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.0),
      axis.text = element_text(color = "black", face = "bold"),
      axis.title.x = element_text(face = "bold"),
      legend.position = "right",
      plot.margin = margin(40, 15, 12, 15)
    )

  ggsave(out_png, p, width = 10.5, height = 6.0, dpi = 320)
  ggsave(out_pdf, p, width = 10.5, height = 6.0)

  counts
}

save_fig2c_heatmap <- function(fpkm_df, sample_info, deg_union, out_png, out_pdf, out_clusters) {
  message("[I] Plotting Fig2C")

  sample_order <- sample_info$sample

  heat_df <- fpkm_df %>%
    filter(geneId %in% deg_union) %>%
    dplyr::select(geneId, geneName, gene_id, type, specific_type, all_of(sample_order))

  heat_df <- heat_df[match(deg_union, heat_df$geneId), ]
  heat_df <- heat_df[!is.na(heat_df$geneId), ]

  expr_mat <- as.matrix(heat_df[, sample_order, drop = FALSE])
  rownames(expr_mat) <- heat_df$geneId

  log2_mat <- log2(expr_mat + 1)
  z_mat <- t(scale(t(log2_mat)))
  z_mat[is.na(z_mat)] <- 0

  row_hc <- hclust(dist(z_mat), method = cluster_method)
  row_clusters <- cutree(row_hc, k = 6)
  cluster_order <- row_clusters[row_hc$order]
  cluster_display_order <- as.integer(unique(cluster_order))
  gaps_row <- which(diff(cluster_order) != 0)

  row_names_ordered <- names(row_clusters)[row_hc$order]
  cluster_df <- tibble(
    geneId = row_names_ordered,
    cluster = as.integer(row_clusters[row_hc$order]),
    row_order = seq_along(row_names_ordered),
    cluster_display_rank = match(as.integer(row_clusters[row_hc$order]), cluster_display_order)
  ) %>%
    left_join(
      heat_df %>% dplyr::select(geneId, geneName, gene_id, type, specific_type),
      by = "geneId"
    ) %>%
    arrange(row_order)
  write_tsv(cluster_df, out_clusters)
  attr(cluster_df, "cluster_display_order") <- cluster_display_order

  col_annotation <- data.frame(
    Day = factor(as.character(sample_info$condition), levels = c("D0", "D14", "D28")),
    Replicate = factor(sample_info$replicate),
    row.names = sample_info$sample,
    stringsAsFactors = TRUE
  )

  row_annotation <- data.frame(
    Cluster = factor(row_clusters, levels = 1:6),
    row.names = names(row_clusters),
    stringsAsFactors = TRUE
  )

  ann_colors <- list(
    Day = c(D0 = "#8FD3FF", D14 = "#FFD076", D28 = "#FFA16C"),
    Replicate = c("1" = "#F0F0F0", "2" = "#D9D9D9", "3" = "#BDBDBD", "4" = "#969696"),
    Cluster = c(
      "1" = "#FF1A1A",
      "2" = "#23C9D6",
      "3" = "#E643FF",
      "4" = "#18C018",
      "5" = "#203DFF",
      "6" = "#000000"
    )
  )

  heat_colors <- colorRampPalette(c("#0B2E8A", "#F7F7F7", "#B40426"))(256)
  heat_breaks <- seq(-2.5, 2.5, length.out = 257)

  pheatmap(
    z_mat,
    cluster_rows = row_hc,
    cluster_cols = FALSE,
    cutree_rows = 6,
    gaps_row = gaps_row,
    show_rownames = FALSE,
    show_colnames = FALSE,
    annotation_col = col_annotation,
    annotation_row = row_annotation,
    annotation_colors = ann_colors,
    color = heat_colors,
    breaks = heat_breaks,
    border_color = NA,
    treeheight_row = 70,
    filename = out_png,
    width = 11,
    height = 16
  )

  pheatmap(
    z_mat,
    cluster_rows = row_hc,
    cluster_cols = FALSE,
    cutree_rows = 6,
    gaps_row = gaps_row,
    show_rownames = FALSE,
    show_colnames = FALSE,
    annotation_col = col_annotation,
    annotation_row = row_annotation,
    annotation_colors = ann_colors,
    color = heat_colors,
    breaks = heat_breaks,
    border_color = NA,
    treeheight_row = 70,
    filename = out_pdf,
    width = 11,
    height = 16
  )

  cluster_df
}

assert_enrichment_pkgs <- function(include_reactome = FALSE) {
  pkgs <- c("clusterProfiler", "AnnotationDbi", "org.Hs.eg.db")
  if (isTRUE(include_reactome)) {
    pkgs <- c(pkgs, "ReactomePA")
  }
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Fig2C right-panel enrichment requires packages: ",
      paste(pkgs, collapse = ", "),
      ". Missing: ",
      paste(missing, collapse = ", ")
    )
  }
}

run_fig2c_cluster_enrichment <- function(cluster_df, terms_per_cluster, out_raw_tsv = NULL, include_reactome = FALSE) {
  message("[I] Running cluster enrichment for Fig2C right panel")
  assert_enrichment_pkgs(include_reactome = include_reactome)

  id_map <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(cluster_df$gene_id),
    keytype = "ENSEMBL",
    columns = c("ENSEMBL", "ENTREZID", "SYMBOL")
  ) %>%
    as_tibble() %>%
    filter(!is.na(ENSEMBL), !is.na(ENTREZID)) %>%
    distinct(ENSEMBL, ENTREZID, SYMBOL)

  cluster_ids <- cluster_df %>%
    filter(!is.na(gene_id)) %>%
    distinct(cluster, gene_id) %>%
    inner_join(id_map, by = c("gene_id" = "ENSEMBL")) %>%
    distinct(cluster, gene_id, ENTREZID, SYMBOL)

  universe_entrez <- unique(cluster_ids$ENTREZID)

  parse_enrich_result <- function(x, cluster_i, source_name) {
    if (is.null(x)) return(tibble())
    res <- as_tibble(x@result)
    if (nrow(res) == 0) return(tibble())
    res %>%
      filter(!is.na(p.adjust), p.adjust < 0.05, Count > 0) %>%
      transmute(
        cluster = as.integer(cluster_i),
        source = source_name,
        term = Description,
        n_genes = as.integer(Count),
        p.adjust = as.numeric(p.adjust),
        neglog10_fdr = -log10(p.adjust)
      )
  }

  enrich_rows <- list()
  for (cluster_i in sort(unique(cluster_ids$cluster))) {
    genes_i <- unique(cluster_ids$ENTREZID[cluster_ids$cluster == cluster_i])
    if (length(genes_i) < 10) next

    go_bp <- tryCatch(
      clusterProfiler::enrichGO(
        gene = genes_i,
        universe = universe_entrez,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      ),
      error = function(e) NULL
    )

    reactome <- NULL
    if (isTRUE(include_reactome)) {
      reactome <- tryCatch(
        ReactomePA::enrichPathway(
          gene = genes_i,
          universe = universe_entrez,
          organism = "human",
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2,
          readable = TRUE
        ),
        error = function(e) NULL
      )
    }

    one_cluster <- bind_rows(
      parse_enrich_result(go_bp, cluster_i, "GO:BP"),
      parse_enrich_result(reactome, cluster_i, "Reactome")
    ) %>%
      distinct(term, .keep_all = TRUE) %>%
      arrange(p.adjust, desc(n_genes))

    n_keep <- as.integer(terms_per_cluster[[as.character(cluster_i)]])
    if (is.na(n_keep) || n_keep <= 0) n_keep <- 3L
    one_cluster <- slice_head(one_cluster, n = n_keep)

    if (nrow(one_cluster) > 0) {
      enrich_rows[[as.character(cluster_i)]] <- one_cluster
    }
  }

  used_terms <- bind_rows(enrich_rows) %>%
    arrange(cluster, p.adjust)

  if (nrow(used_terms) == 0) {
    stop("No significant terms found for any cluster (FDR < 0.05).")
  }

  if (!is.null(out_raw_tsv)) {
    write_tsv(used_terms, out_raw_tsv)
  }

  used_terms
}

save_fig2c_right_panel <- function(terms_df, out_png, out_pdf, out_tsv, cluster_display_order = 1:6) {
  message("[I] Plotting Fig2C right panel (from own enrichment)")

  cluster_colors <- c(
    `1` = "#FF1A1A",
    `2` = "#56CFC7",
    `3` = "#E643FF",
    `4` = "#18C018",
    `5` = "#203DFF",
    `6` = "#000000"
  )

  df <- terms_df %>%
    mutate(
      cluster = factor(cluster, levels = as.character(cluster_display_order)),
      term_label = sprintf("%s (%s)", term, n_genes)
    ) %>%
    group_by(cluster) %>%
    arrange(p.adjust, desc(n_genes), .by_group = TRUE) %>%
    mutate(
      term_order = row_number(),
      term_id = paste0("C", cluster, "_", term_order),
      xmax = max(neglog10_fdr, na.rm = TRUE),
      xmax = if_else(!is.finite(xmax) | xmax < 1, 1, xmax),
      text_x = neglog10_fdr + pmax(0.22, xmax * 0.03)
    ) %>%
    ungroup()

  term_levels <- df %>%
    arrange(cluster, desc(term_order)) %>%
    pull(term_id) %>%
    unique()
  df$term_id <- factor(df$term_id, levels = term_levels)

  write_tsv(
    df %>% dplyr::select(cluster, source, term, n_genes, p.adjust, neglog10_fdr),
    out_tsv
  )

  p <- ggplot(df, aes(x = neglog10_fdr, y = term_id, fill = cluster)) +
    geom_col(width = 0.68, orientation = "y") +
    geom_text(
      aes(x = text_x, label = term_label),
      hjust = 0,
      size = 4.9,
      lineheight = 0.98,
      color = "black",
      show.legend = FALSE
    ) +
    facet_wrap(~ cluster, ncol = 1, scales = "free", strip.position = "left") +
    scale_fill_manual(values = cluster_colors) +
    scale_x_continuous(
      name = "-Log10 FDR        (number of genes)",
      position = "top",
      expand = expansion(mult = c(0, 0.34))
    ) +
    scale_y_discrete(labels = rep("", nlevels(df$term_id))) +
    coord_cartesian(clip = "off") +
    labs(y = NULL) +
    theme_classic(base_size = 13.5) +
    theme(
      plot.background = element_rect(fill = "#ECECEC", color = NA),
      panel.background = element_rect(fill = "#ECECEC", color = NA),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      axis.title.x = element_text(face = "bold", size = 17, margin = margin(b = 8)),
      axis.text.x = element_text(face = "bold", color = "black", size = 12.5),
      strip.background = element_blank(),
      strip.text = element_blank(),
      strip.placement = "outside",
      panel.spacing.y = unit(0.12, "lines"),
      legend.position = "none",
      plot.margin = margin(6, 110, 6, 8)
    )

  ggsave(out_png, p, width = 13.2, height = 15.8, dpi = 320)
  ggsave(out_pdf, p, width = 13.2, height = 15.8)
}

# ============================================================
# Main workflow
# ============================================================

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

count_obj <- read_expr_matrix(count_file, "count")
fpkm_obj <- read_expr_matrix(fpkm_file, "FPKM")
counts_df <- count_obj$df
fpkm_df <- fpkm_obj$df

sample_cols <- intersect(count_obj$sample_cols, fpkm_obj$sample_cols)
sample_info <- build_sample_info(sample_cols)

sample_cols <- sample_info$sample
counts_df <- counts_df %>%
  dplyr::select(geneId, geneName, gene_id, any_of(c("type", "specific_type")), all_of(sample_cols))
fpkm_df <- fpkm_df %>%
  dplyr::select(geneId, geneName, gene_id, any_of(c("type", "specific_type")), all_of(sample_cols))

biotype_map <- build_biotype_map(counts_df)

comparison_defs <- list(
  "D0 vs D14" = c("D0", "D14"),
  "D14 vs D28" = c("D14", "D28"),
  "D0 vs D28" = c("D0", "D28")
)

pair_results <- list()
for (lbl in names(comparison_defs)) {
  pair <- comparison_defs[[lbl]]
  pair_results[[lbl]] <- run_pair_deseq(
    label = lbl,
    group_a = pair[[1]],
    group_b = pair[[2]],
    counts_df = counts_df,
    fpkm_df = fpkm_df,
    sample_info = sample_info,
    biotype_map = biotype_map
  )

  base_lbl <- gsub(" ", "_", gsub(" vs ", "_vs_", lbl, fixed = TRUE), fixed = TRUE)
  write_tsv(pair_results[[lbl]]$full, file.path(outdir, paste0("deseq2_", base_lbl, ".tsv")))
  write_tsv(pair_results[[lbl]]$deg, file.path(outdir, paste0("degs_", base_lbl, ".tsv")))
}

# Fig2A
venn_sets <- lapply(pair_results, function(x) unique(x$deg$geneId))
save_fig2a_venn(
  gene_sets = venn_sets,
  out_png = file.path(outdir, "Fig2A_venn.png"),
  out_pdf = file.path(outdir, "Fig2A_venn.pdf")
)

# Fig2B
fig2b_counts <- save_fig2b_barplot(
  deg_list = lapply(pair_results, function(x) x$deg),
  out_png = file.path(outdir, "Fig2B_up_down_barplot.png"),
  out_pdf = file.path(outdir, "Fig2B_up_down_barplot.pdf"),
  out_tsv = file.path(outdir, "Fig2B_counts.tsv")
)

paper_counts <- tribble(
  ~comparison, ~direction, ~biotype, ~n_paper,
  "D0 vs D14", "Down", "Protein-coding", 993L,
  "D0 vs D14", "Down", "lncRNA", 201L,
  "D0 vs D14", "Up", "Protein-coding", 835L,
  "D0 vs D14", "Up", "lncRNA", 168L,
  "D14 vs D28", "Down", "Protein-coding", 441L,
  "D14 vs D28", "Down", "lncRNA", 58L,
  "D14 vs D28", "Up", "Protein-coding", 525L,
  "D14 vs D28", "Up", "lncRNA", 50L,
  "D0 vs D28", "Down", "Protein-coding", 1213L,
  "D0 vs D28", "Down", "lncRNA", 255L,
  "D0 vs D28", "Up", "Protein-coding", 1357L,
  "D0 vs D28", "Up", "lncRNA", 224L
)

fig2b_vs_paper <- fig2b_counts %>%
  dplyr::select(comparison, direction, biotype, n) %>%
  left_join(paper_counts, by = c("comparison", "direction", "biotype")) %>%
  mutate(
    delta = n - n_paper,
    pct_of_paper = round(100 * n / n_paper, 2)
  ) %>%
  arrange(comparison, direction, biotype)
write_tsv(fig2b_vs_paper, file.path(outdir, "Fig2B_vs_paper_counts.tsv"))

# Fig2C
deg_union <- unique(unlist(lapply(pair_results, function(x) x$deg$geneId)))
cluster_df <- save_fig2c_heatmap(
  fpkm_df = fpkm_df,
  sample_info = sample_info,
  deg_union = deg_union,
  out_png = file.path(outdir, "Fig2C_heatmap.png"),
  out_pdf = file.path(outdir, "Fig2C_heatmap.pdf"),
  out_clusters = file.path(outdir, "Fig2C_clusters.tsv")
)

if (isTRUE(export_fig2c_right_panel)) {
  cluster_display_order <- attr(cluster_df, "cluster_display_order")
  if (is.null(cluster_display_order) || length(cluster_display_order) == 0) {
    cluster_display_order <- 1:6
  }

  fig2c_terms <- run_fig2c_cluster_enrichment(
    cluster_df = cluster_df,
    terms_per_cluster = fig2c_right_terms_per_cluster,
    out_raw_tsv = file.path(outdir, "Fig2C_right_panel_enrichment_raw.tsv"),
    include_reactome = fig2c_include_reactome
  )
  save_fig2c_right_panel(
    terms_df = fig2c_terms,
    out_png = file.path(outdir, "Fig2C_right_panel_own_enrichment.png"),
    out_pdf = file.path(outdir, "Fig2C_right_panel_own_enrichment.pdf"),
    out_tsv = file.path(outdir, "Fig2C_right_panel_own_enrichment.tsv"),
    cluster_display_order = cluster_display_order
  )
}

summary_tbl <- tibble(
  comparison = names(pair_results),
  n_deg = as.integer(sapply(pair_results, function(x) nrow(x$deg))),
  n_protein_coding = as.integer(sapply(pair_results, function(x) sum(x$deg$biotype == "Protein-coding", na.rm = TRUE))),
  n_lncRNA = as.integer(sapply(pair_results, function(x) sum(x$deg$biotype == "lncRNA", na.rm = TRUE)))
)
write_tsv(summary_tbl, file.path(outdir, "Fig2_summary_counts.tsv"))

run_parameters <- tibble(
  parameter = c(
    "count_file", "fpkm_file", "outdir",
    "padj_cutoff", "fc_cutoff", "min_fpkm",
    "fitType", "sfType", "independentFiltering", "cooksCutoff",
    "design", "cluster_method", "export_fig2c_right_panel", "fig2c_right_terms_per_cluster", "fig2c_include_reactome"
  ),
  value = c(
    count_file, fpkm_file, outdir,
    as.character(padj_cutoff), as.character(fc_cutoff), as.character(min_fpkm),
    fit_type, sf_type, as.character(independent_filtering), as.character(cooks_cutoff),
    "~ condition", cluster_method, as.character(export_fig2c_right_panel),
    paste(names(fig2c_right_terms_per_cluster), fig2c_right_terms_per_cluster, sep = ":", collapse = ","),
    as.character(fig2c_include_reactome)
  )
)
write_tsv(run_parameters, file.path(outdir, "run_parameters.tsv"))

message("[OK] Finished.")
message("[OK] Outputs written to: ", normalizePath(outdir))
print(summary_tbl)
