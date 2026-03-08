#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(forcats)
})

# ============================================================
# Figure 3 reproduction (data-driven from your own matrices)
# ============================================================

count_file <- "data/GSE153005_count_mat.txt"
fpkm_file <- "data/GSE153005_FPKM_mat.txt"
outdir <- "results/fig3"

padj_cutoff <- 0.05
fc_cutoff <- 4
min_fpkm <- 1

fit_type <- "local"
sf_type <- "poscounts"
independent_filtering <- TRUE
cooks_cutoff <- TRUE

top_ns <- c(10, 20, 30, 40, 50)
lnc_classes <- c("lincRNA", "antisense", "sense_intronic", "TEC", "processed_transcript", "sense_overlapping")

# Fig3A strict-paper settings:
# Paper states "We used the absolute FC to analyze the top ranked DEGs".
# Here we operationalize ranking by absolute DESeq2 log2 fold-change.
fig3a_rank_mode <- "abs_log2FoldChange"
fig3a_hypergeom_universe <- "comparison_DEGs"

set.seed(20260304)

if (!file.exists(count_file) || !file.exists(fpkm_file)) {
  stop(
    "Input matrices not found. Please place files at: ",
    count_file, " and ", fpkm_file
  )
}

# ============================================================
# Helpers
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
    stop("Count matrix is missing 'type' column.")
  }

  count_df %>%
    transmute(
      geneId,
      gene_id,
      geneName,
      raw_biotype = as.character(type),
      specific_type = as.character(specific_type),
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
    arrange(desc(fc_for_filter), padj)

  list(full = full_tbl, deg = deg_tbl)
}

hypergeom_right_tail <- function(x, K, N, n) {
  if (is.na(x) || is.na(K) || is.na(N) || is.na(n) || N <= 0 || K < 0 || n < 0 || x < 0) return(NA_real_)
  phyper(q = x - 1, m = K, n = N - K, k = n, lower.tail = FALSE)
}

infer_target_gene <- function(lnc_name) {
  x <- lnc_name
  x <- sub("-AS[0-9]+$", "", x, perl = TRUE)
  x <- sub("AS[0-9]+$", "", x, perl = TRUE)
  x <- sub("-OT$", "", x, perl = TRUE)
  x <- sub("OT$", "", x, perl = TRUE)
  x <- sub("OS$", "", x, perl = TRUE)
  x
}

# ============================================================
# Data loading + DE
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
}

# ============================================================
# Fig 3A
# ============================================================

if (!fig3a_rank_mode %in% c("abs_log2FoldChange", "normalized_count_fc")) {
  stop("fig3a_rank_mode must be one of: abs_log2FoldChange, normalized_count_fc")
}
if (!fig3a_hypergeom_universe %in% c("comparison_DEGs", "all_annotated_genes")) {
  stop("fig3a_hypergeom_universe must be one of: comparison_DEGs, all_annotated_genes")
}

global_bg_tbl <- biotype_map %>%
  filter(biotype %in% c("Protein-coding", "lncRNA")) %>%
  distinct(geneId, biotype)

fig3a_ranked_genes <- bind_rows(lapply(names(pair_results), function(cmp) {
  deg_raw <- pair_results[[cmp]]$deg
  deg_ranked <- if (fig3a_rank_mode == "abs_log2FoldChange") {
    deg_raw %>%
      mutate(rank_score = abs(log2FoldChange)) %>%
      arrange(desc(rank_score), pvalue, padj, desc(fc_for_filter))
  } else {
    deg_raw %>%
      mutate(rank_score = fc_for_filter) %>%
      arrange(desc(rank_score), pvalue, padj)
  }

  deg_ranked %>%
    mutate(comparison = cmp, rank_idx = row_number())
}))

write_tsv(
  fig3a_ranked_genes %>%
    filter(rank_idx <= max(top_ns)) %>%
    dplyr::select(comparison, rank_idx, geneId, geneName, biotype, direction, rank_score, log2FoldChange, fc_for_filter, padj, pvalue),
  file.path(outdir, "Fig3A_top_ranked_genes.tsv")
)

fig3a_counts <- bind_rows(lapply(names(pair_results), function(cmp) {
  deg <- fig3a_ranked_genes %>%
    filter(comparison == cmp) %>%
    arrange(rank_idx)

  bg_cmp <- if (fig3a_hypergeom_universe == "comparison_DEGs") {
    deg %>% distinct(geneId, biotype)
  } else {
    global_bg_tbl
  }
  N_bg <- nrow(bg_cmp)
  K_lnc <- sum(bg_cmp$biotype == "lncRNA", na.rm = TRUE)
  K_pc <- sum(bg_cmp$biotype == "Protein-coding", na.rm = TRUE)

  bind_rows(lapply(top_ns, function(top_n) {
    top_tbl <- slice_head(deg, n = top_n)
    n_draw <- nrow(top_tbl)
    n_lnc <- sum(top_tbl$biotype == "lncRNA", na.rm = TRUE)
    n_pc <- sum(top_tbl$biotype == "Protein-coding", na.rm = TRUE)

    p_lnc <- hypergeom_right_tail(n_lnc, K_lnc, N_bg, n_draw)
    p_pc <- hypergeom_right_tail(n_pc, K_pc, N_bg, n_draw)

    cat_tbl <- top_tbl %>%
      mutate(
        category = case_when(
          biotype == "lncRNA" & direction == "Down" ~ "lncRNA down",
          biotype == "lncRNA" & direction == "Up" ~ "lncRNA up",
          biotype == "Protein-coding" & direction == "Up" ~ "Protein-coding up",
          biotype == "Protein-coding" & direction == "Down" ~ "Protein-coding down",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(category)) %>%
      count(category, name = "n") %>%
      mutate(
        comparison = cmp,
        top_n = top_n,
        p_lnc = p_lnc,
        p_pc = p_pc,
        star = if_else(!is.na(p_lnc) & !is.na(p_pc) & p_lnc < p_pc, "*", "")
      )

    cat_tbl
  }))
}))

fig3a_counts <- fig3a_counts %>%
  complete(
    comparison = names(comparison_defs),
    top_n = top_ns,
    category = c("lncRNA down", "lncRNA up", "Protein-coding up", "Protein-coding down"),
    fill = list(n = 0, p_lnc = NA_real_, p_pc = NA_real_, star = "")
  ) %>%
  mutate(
    comparison = factor(comparison, levels = c("D0 vs D14", "D0 vs D28", "D14 vs D28")),
    top_n = factor(top_n, levels = top_ns),
    category = factor(category, levels = c("lncRNA down", "lncRNA up", "Protein-coding up", "Protein-coding down"))
  )

fig3a_plot_df <- fig3a_counts %>%
  mutate(
    biotype_group = if_else(str_detect(as.character(category), "^lncRNA"), "lncRNA", "Protein-coding"),
    direction_group = if_else(str_detect(as.character(category), "up$"), "up", "down"),
    top_n_num = as.numeric(as.character(top_n)),
    top_n_idx = match(top_n_num, top_ns)
  ) %>%
  group_by(comparison, top_n, top_n_num, top_n_idx, biotype_group) %>%
  arrange(
    case_when(
      biotype_group == "lncRNA" & direction_group == "down" ~ 1L,
      biotype_group == "lncRNA" & direction_group == "up" ~ 2L,
      biotype_group == "Protein-coding" & direction_group == "up" ~ 1L,
      biotype_group == "Protein-coding" & direction_group == "down" ~ 2L,
      TRUE ~ 3L
    ),
    .by_group = TRUE
  ) %>%
  mutate(
    y_min = lag(cumsum(n), default = 0),
    y_max = cumsum(n)
  ) %>%
  ungroup() %>%
  mutate(
    x_center = top_n_idx + if_else(biotype_group == "lncRNA", -0.21, 0.21),
    bar_width = 0.35,
    xmin = x_center - bar_width / 2,
    xmax = x_center + bar_width / 2
  )

fig3a_totals <- fig3a_plot_df %>%
  group_by(comparison, top_n, top_n_num, top_n_idx, biotype_group) %>%
  summarise(total = max(y_max), .groups = "drop")

fig3a_stars <- fig3a_counts %>%
  group_by(comparison, top_n) %>%
  summarise(
    star = if_else(any(star != ""), first(star[star != ""]), ""),
    .groups = "drop"
  ) %>%
  left_join(
    fig3a_totals %>%
      group_by(comparison, top_n, top_n_num, top_n_idx) %>%
      summarise(y_star = max(total) + 1.2, .groups = "drop"),
    by = c("comparison", "top_n")
  )

write_tsv(fig3a_counts, file.path(outdir, "Fig3A_top_counts.tsv"))
write_tsv(fig3a_stars, file.path(outdir, "Fig3A_hypergeom_stars.tsv"))

fig3a_colors <- c(
  "lncRNA down" = "#FF4B19",
  "lncRNA up" = "#B40000",
  "Protein-coding up" = "#2A62D5",
  "Protein-coding down" = "#8ED4F9"
)

p3a_ymax <- ceiling((max(fig3a_plot_df$y_max, na.rm = TRUE) + 3) / 5) * 5

p3a <- ggplot(fig3a_plot_df) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = y_min, ymax = y_max, fill = category),
    color = NA
  ) +
  geom_text(
    data = fig3a_stars %>% filter(star != ""),
    aes(x = top_n_idx, y = y_star, label = star),
    inherit.aes = FALSE,
    size = 7,
    fontface = "bold"
  ) +
  facet_wrap(~ comparison, nrow = 1) +
  scale_fill_manual(values = fig3a_colors) +
  scale_x_continuous(
    breaks = seq_along(top_ns),
    labels = top_ns,
    limits = c(0.5, length(top_ns) + 0.5)
  ) +
  scale_y_continuous(
    limits = c(0, p3a_ymax),
    breaks = seq(0, p3a_ymax, by = 5)
  ) +
  labs(x = "Top DEGs", y = "Number of DEGs", fill = "Type") +
  theme_bw(base_size = 13) +
  theme(
    panel.background = element_rect(fill = "#F0F0F0", color = "black", linewidth = 0.8),
    panel.grid.major = element_line(color = "#DDDDDD", linewidth = 0.45),
    panel.grid.minor = element_line(color = "#E9E9E9", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(fill = "#F0F0F0", color = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

ggsave(file.path(outdir, "Fig3A_top_degs.png"), p3a, width = 12.8, height = 4.8, dpi = 320)
ggsave(file.path(outdir, "Fig3A_top_degs.pdf"), p3a, width = 12.8, height = 4.8)

# ============================================================
# Fig 3B
# ============================================================

lnc_bg <- biotype_map %>%
  filter(biotype == "lncRNA") %>%
  distinct(geneId, specific_type)

N_lnc_bg <- nrow(lnc_bg)

fig3b_counts <- bind_rows(lapply(names(pair_results), function(cmp) {
  lnc_deg <- pair_results[[cmp]]$deg %>%
    filter(biotype == "lncRNA") %>%
    distinct(geneId, specific_type)

  n_draw <- nrow(lnc_deg)

  bind_rows(lapply(lnc_classes, function(cls) {
    x <- sum(lnc_deg$specific_type == cls, na.rm = TRUE)
    K <- sum(lnc_bg$specific_type == cls, na.rm = TRUE)
    p <- hypergeom_right_tail(x, K, N_lnc_bg, n_draw)
    star <- case_when(
      !is.na(p) & p < 0.001 ~ "***",
      !is.na(p) & p < 0.05 ~ "*",
      TRUE ~ ""
    )
    tibble(
      comparison = cmp,
      specific_type = cls,
      n = x,
      p_value = p,
      star = star
    )
  }))
}))

fig3b_counts <- fig3b_counts %>%
  mutate(
    comparison = factor(comparison, levels = c("D0 vs D14", "D0 vs D28", "D14 vs D28")),
    specific_type = factor(specific_type, levels = lnc_classes)
  )

write_tsv(fig3b_counts, file.path(outdir, "Fig3B_lnc_class_counts.tsv"))

fig3b_colors <- c(
  "lincRNA" = "#143CFF",
  "antisense" = "#1EA7B8",
  "sense_intronic" = "#E4DE00",
  "TEC" = "#31C828",
  "processed_transcript" = "#2BA5A0",
  "sense_overlapping" = "#F08C00"
)

p3b <- ggplot(fig3b_counts, aes(x = specific_type, y = n, fill = specific_type)) +
  geom_col(width = 0.72) +
  geom_text(
    aes(y = n + max(n) * 0.04 + 0.8, label = star),
    size = 5,
    fontface = "bold"
  ) +
  facet_wrap(~ comparison, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = fig3b_colors) +
  scale_x_discrete(labels = c(
    "lincRNA" = "lincRNA",
    "antisense" = "Antisense",
    "sense_intronic" = "Sense intronic",
    "TEC" = "TEC",
    "processed_transcript" = "Processed transcript",
    "sense_overlapping" = "Sense overlapping"
  )) +
  labs(x = NULL, y = "Number of DE lncRNAs") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

ggsave(file.path(outdir, "Fig3B_lnc_classes.png"), p3b, width = 12.8, height = 4.8, dpi = 320)
ggsave(file.path(outdir, "Fig3B_lnc_classes.pdf"), p3b, width = 12.8, height = 4.8)

# ============================================================
# Fig 3C
# ============================================================

message("[I] Building lncRNA-target pairs for Fig3C")

gene_name_map <- biotype_map %>%
  filter(!is.na(geneName), geneName != "") %>%
  group_by(geneName) %>%
  arrange(desc(biotype == "Protein-coding"), geneId) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(geneName, geneId, biotype, gene_id)

lnc_genes <- biotype_map %>%
  filter(biotype == "lncRNA") %>%
  dplyr::select(lnc_geneId = geneId, lnc_gene_id = gene_id, lnc_geneName = geneName)

pairs_tbl <- lnc_genes %>%
  mutate(target_geneName = infer_target_gene(lnc_geneName)) %>%
  inner_join(
    gene_name_map %>%
      filter(biotype == "Protein-coding") %>%
      transmute(target_geneName = geneName, target_geneId = geneId, target_gene_id = gene_id),
    by = "target_geneName"
  ) %>%
  distinct(lnc_geneId, target_geneId, .keep_all = TRUE)

if (nrow(pairs_tbl) == 0) {
  warning("[W] No lncRNA-target pairs inferred by gene-name heuristics. Fig3C will be skipped.")
} else {
  fpkm_expr <- fpkm_df %>% dplyr::select(geneId, all_of(sample_cols))

  calc_pair_cor <- function(lnc_id, pc_id) {
    x <- as.numeric(fpkm_expr[fpkm_expr$geneId == lnc_id, sample_cols, drop = TRUE])
    y <- as.numeric(fpkm_expr[fpkm_expr$geneId == pc_id, sample_cols, drop = TRUE])
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 4) {
      return(tibble(cor = NA_real_, p_cor = NA_real_))
    }
    ct <- suppressWarnings(cor.test(log2(x[ok] + 1), log2(y[ok] + 1), method = "pearson"))
    tibble(cor = unname(ct$estimate), p_cor = ct$p.value)
  }

  cor_tbl <- bind_rows(lapply(seq_len(nrow(pairs_tbl)), function(i) {
    one <- pairs_tbl[i, ]
    cdat <- calc_pair_cor(one$lnc_geneId, one$target_geneId)
    bind_cols(one, cdat)
  })) %>%
    filter(!is.na(cor), !is.na(p_cor), p_cor < 0.05)

  write_tsv(cor_tbl, file.path(outdir, "Fig3C_pairs_all.tsv"))

  fig3c_colors <- c(
    "Positive Up" = "#F22A2A",
    "Positive Down" = "#2D58F4",
    "Negative" = "#A8A8A8"
  )
  fig3c_axis_ticks <- seq(-10, 10, by = 2)
  fig3c_tick_len <- 0.32

  build_fig3c_annotations <- function(cmp_label, points_tbl) {
    ann_spec <- if (cmp_label == "D0 vs D14") {
      tibble(
        lnc_geneName = c("SOX1-OT", "EPHA1-AS1"),
        label = c("SOX1 and SOX1-OT", "EPHA1 and\nEPHA1-AS1"),
        label_x = c(1.3, -7.3),
        label_y = c(9.7, -8.2),
        line_x = c(2.3, -5.8),
        line_y = c(9.1, -7.0),
        hjust = c(0, 0),
        vjust = c(0.5, 0.5)
      )
    } else if (cmp_label == "D0 vs D28") {
      tibble(
        lnc_geneName = c("DIO3OS", "CACNA1C-AS2", "FOXD3-AS1"),
        label = c("DIO3 and\nDIO3OS", "CACNA1C and\nCACNA1C-AS2", "FOXD3\nand\nFOXD3-AS1"),
        label_x = c(2.2, -9.0, -6.7),
        label_y = c(10.1, 6.0, -8.1),
        line_x = c(5.6, -6.4, -6.0),
        line_y = c(10.0, 5.1, -6.0),
        hjust = c(0, 0, 0),
        vjust = c(0.5, 0.5, 0.5)
      )
    } else {
      tibble()
    }

    if (nrow(ann_spec) == 0) {
      return(ann_spec)
    }

    ann_spec %>%
      left_join(
        points_tbl %>%
          dplyr::select(lnc_geneName, target_geneName, lfc_lnc, lfc_pc, corr_type) %>%
          distinct(lnc_geneName, .keep_all = TRUE),
        by = "lnc_geneName"
      ) %>%
      filter(!is.na(lfc_lnc), !is.na(lfc_pc))
  }

  make_fig3c_scatter <- function(cmp_label, out_png, out_pdf) {
    lfc_tbl <- pair_results[[cmp_label]]$full %>%
      dplyr::select(geneId, log2FoldChange, is_deg) %>%
      rename(lfc = log2FoldChange, is_deg_cmp = is_deg)

    points <- cor_tbl %>%
      left_join(
        lfc_tbl %>%
          rename(lnc_geneId = geneId, lfc_lnc = lfc, lnc_deg_cmp = is_deg_cmp),
        by = "lnc_geneId"
      ) %>%
      left_join(
        lfc_tbl %>%
          rename(target_geneId = geneId, lfc_pc = lfc, pc_deg_cmp = is_deg_cmp),
        by = "target_geneId"
      ) %>%
      filter(!is.na(lfc_lnc), !is.na(lfc_pc)) %>%
      mutate(
        both_deg = lnc_deg_cmp & pc_deg_cmp,
        corr_type = case_when(
          cor < 0 ~ "Negative",
          cor >= 0 & lfc_lnc >= 0 & lfc_pc >= 0 ~ "Positive Up",
          cor >= 0 & lfc_lnc < 0 & lfc_pc < 0 ~ "Positive Down",
          TRUE ~ "Negative"
        ),
        corr_type = factor(corr_type, levels = c("Negative", "Positive Down", "Positive Up"))
      )

    write_tsv(
      points,
      file.path(outdir, paste0("Fig3C_points_", gsub(" ", "_", gsub(" vs ", "_vs_", cmp_label, fixed = TRUE), fixed = TRUE), ".tsv"))
    )

    if (nrow(points) == 0) {
      warning("[W] No points for ", cmp_label, " in Fig3C.")
      return(NULL)
    }

    x_tick_df <- tibble(x = fig3c_axis_ticks)
    y_tick_df <- tibble(y = fig3c_axis_ticks)
    x_lab_df <- tibble(
      x = fig3c_axis_ticks[fig3c_axis_ticks != 0],
      y = -1.02,
      lab = as.character(fig3c_axis_ticks[fig3c_axis_ticks != 0])
    )
    y_lab_df <- tibble(
      y = fig3c_axis_ticks[fig3c_axis_ticks != 0],
      x = -1.02,
      lab = as.character(fig3c_axis_ticks[fig3c_axis_ticks != 0])
    )
    ann_df <- build_fig3c_annotations(cmp_label, points)

    p <- ggplot(points, aes(x = lfc_lnc, y = lfc_pc)) +
      geom_segment(
        data = tibble(x = -10, y = -10, xend = 10, yend = 10),
        aes(x = x, y = y, xend = xend, yend = yend),
        inherit.aes = FALSE,
        linetype = "22",
        linewidth = 0.85,
        color = "black"
      ) +
      geom_segment(
        data = tibble(x = -10, y = 0, xend = 10, yend = 0),
        aes(x = x, y = y, xend = xend, yend = yend),
        inherit.aes = FALSE,
        linewidth = 1.05,
        color = "black"
      ) +
      geom_segment(
        data = tibble(x = 0, y = -10, xend = 0, yend = 10),
        aes(x = x, y = y, xend = xend, yend = yend),
        inherit.aes = FALSE,
        linewidth = 1.05,
        color = "black"
      ) +
      geom_segment(
        data = x_tick_df,
        aes(x = x, y = -fig3c_tick_len, xend = x, yend = fig3c_tick_len),
        inherit.aes = FALSE,
        linewidth = 0.95,
        color = "black"
      ) +
      geom_segment(
        data = y_tick_df,
        aes(x = -fig3c_tick_len, y = y, xend = fig3c_tick_len, yend = y),
        inherit.aes = FALSE,
        linewidth = 0.95,
        color = "black"
      ) +
      geom_text(
        data = x_lab_df,
        aes(x = x, y = y, label = lab),
        inherit.aes = FALSE,
        size = 4.2
      ) +
      geom_text(
        data = y_lab_df,
        aes(x = x, y = y, label = lab),
        inherit.aes = FALSE,
        size = 4.2
      ) +
      geom_point(
        data = points %>% filter(!both_deg),
        aes(color = corr_type),
        shape = 21,
        fill = "white",
        size = 3.7,
        stroke = 1.15,
        alpha = 1
      ) +
      geom_point(
        data = points %>% filter(both_deg),
        aes(color = corr_type),
        shape = 16,
        size = 3.9,
        stroke = 0,
        alpha = 1
      ) +
      scale_color_manual(values = fig3c_colors, drop = FALSE) +
      coord_fixed(xlim = c(-10, 10), ylim = c(-10.8, 10.8), clip = "off") +
      labs(
        x = paste0("Log2 FC ", cmp_label, " (lncRNAs)"),
        y = paste0("Log2 FC ", cmp_label, " (protein coding)"),
        title = cmp_label
      ) +
      theme_void(base_size = 14) +
      theme(
        axis.title.x = element_text(size = 12.5, color = "black", margin = margin(t = 10)),
        axis.title.y = element_text(size = 12.5, color = "black", angle = 90, margin = margin(r = 12)),
        plot.title = element_text(size = 19, face = "bold", hjust = 0.5, color = "black", margin = margin(b = 4)),
        legend.position = "none",
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(t = 24, r = 28, b = 12, l = 28)
      )

    if (nrow(ann_df) > 0) {
      p <- p +
        geom_segment(
          data = ann_df,
          aes(x = lfc_lnc, y = lfc_pc, xend = line_x, yend = line_y, color = corr_type),
          inherit.aes = FALSE,
          linewidth = 0.95,
          alpha = 0.9
        ) +
        geom_text(
          data = ann_df,
          aes(x = label_x, y = label_y, label = label, color = corr_type, hjust = hjust, vjust = vjust),
          inherit.aes = FALSE,
          size = 4.6,
          lineheight = 0.96
        )
    }

    if (cmp_label == "D0 vs D14") {
      p <- p + annotate(
        "text",
        x = -14.4,
        y = 11.45,
        label = "C)",
        hjust = 0,
        vjust = 1,
        size = 8.2,
        fontface = "bold",
        color = "black"
      )
    }

    ggsave(out_png, p, width = 6.2, height = 6.0, dpi = 320, bg = "white")
    ggsave(out_pdf, p, width = 6.2, height = 6.0, bg = "white")
    p
  }

  p3c_left <- make_fig3c_scatter(
    cmp_label = "D0 vs D14",
    out_png = file.path(outdir, "Fig3C_scatter_D0_vs_D14.png"),
    out_pdf = file.path(outdir, "Fig3C_scatter_D0_vs_D14.pdf")
  )

  p3c_right <- make_fig3c_scatter(
    cmp_label = "D0 vs D28",
    out_png = file.path(outdir, "Fig3C_scatter_D0_vs_D28.png"),
    out_pdf = file.path(outdir, "Fig3C_scatter_D0_vs_D28.pdf")
  )

  if (requireNamespace("patchwork", quietly = TRUE) && !is.null(p3c_left) && !is.null(p3c_right)) {
    p3c_combined <- p3c_left + p3c_right + patchwork::plot_layout(ncol = 2)
    ggsave(file.path(outdir, "Fig3C_scatter_combined.png"), p3c_combined, width = 12.6, height = 6.0, dpi = 320, bg = "white")
    ggsave(file.path(outdir, "Fig3C_scatter_combined.pdf"), p3c_combined, width = 12.6, height = 6.0, bg = "white")
  }
}

run_parameters <- tibble(
  parameter = c(
    "count_file", "fpkm_file", "outdir",
    "padj_cutoff", "fc_cutoff", "min_fpkm",
    "fitType", "sfType", "independentFiltering", "cooksCutoff",
    "top_ns", "lnc_classes",
    "fig3a_rank_mode", "fig3a_hypergeom_universe"
  ),
  value = c(
    count_file, fpkm_file, outdir,
    as.character(padj_cutoff), as.character(fc_cutoff), as.character(min_fpkm),
    fit_type, sf_type, as.character(independent_filtering), as.character(cooks_cutoff),
    paste(top_ns, collapse = ","), paste(lnc_classes, collapse = ","),
    fig3a_rank_mode, fig3a_hypergeom_universe
  )
)
write_tsv(run_parameters, file.path(outdir, "run_parameters_fig3.tsv"))

message("[OK] Finished Fig3 reproduction.")
message("[OK] Outputs written to: ", normalizePath(outdir))
