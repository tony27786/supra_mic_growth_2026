# =============================================================================
# hetero_cnv.R -- mosdepth summaries -> normalized per-chromosome depth tables.
#
# Run from this directory, after 02_mosdepth and merge_batches.sh:
#   Rscript hetero_cnv.R
#
# Reads  : $MOS_DIR/*.mosdepth.summary.txt   (set in wgs_analysis_config.R)
# Writes : mosdepth_chr_normalized_depth_long.tsv   <- deposited as SourceData
#          mosdepth_chr_normalized_depth_wide.tsv
#          plus CNV hint/summary tables and genome-wide plots
# =============================================================================
source("./wgs_analysis_config.R")

library(tidyverse)
library(ggplot2)

############################################
### STEP 4.1a chromosome-level normalized depth
############################################

# =========================
# Directory
# =========================
mos_dir <- MOS_DIR          # from wgs_analysis_config.R
out_dir <- DEPTH_DIR

out_long <- file.path(out_dir, "mosdepth_chr_normalized_depth_long.tsv")
out_wide <- file.path(out_dir, "mosdepth_chr_normalized_depth_wide.tsv")

# =========================
# Read one mosdepth summary
# =========================
read_chr_depth_one <- function(file) {
  sample <- sub("\\.mosdepth\\.summary\\.txt$", "", basename(file))
  
  x <- read_tsv(
    file,
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    # Remove *_region line
    filter(!str_detect(chrom, "_region$")) %>%
    # Remove total line
    filter(chrom != "total") %>%
    mutate(
      length = as.numeric(length),
      bases = as.numeric(bases),
      mean = as.numeric(mean),
      min = as.numeric(min),
      max = as.numeric(max)
    ) %>%
    filter(!is.na(length), !is.na(bases), !is.na(mean)) %>%
    # Remove small contigs
    mutate(
      is_nuclear = length >= 1e5
    ) %>%
    filter(is_nuclear) %>%
    mutate(
      sample = sample
    ) %>%
    select(sample, chrom, length, bases, mean, min, max)
  
  if (nrow(x) == 0) {
    return(tibble(
      sample = sample,
      chrom = NA_character_,
      chr_length = NA_real_,
      chr_bases = NA_real_,
      chr_mean_depth = NA_real_,
      chr_min_depth = NA_real_,
      chr_max_depth = NA_real_,
      sample_mean_depth_nuclear = NA_real_,
      chr_depth_norm = NA_real_,
      chr_log2_norm = NA_real_
    ))
  }
  
  sample_mean_depth_nuclear <- sum(x$bases, na.rm = TRUE) / sum(x$length, na.rm = TRUE)
  
  x %>%
    transmute(
      sample = sample,
      chrom = chrom,
      chr_length = length,
      chr_bases = bases,
      chr_mean_depth = mean,
      chr_min_depth = min,
      chr_max_depth = max,
      sample_mean_depth_nuclear = sample_mean_depth_nuclear,
      chr_depth_norm = chr_mean_depth / sample_mean_depth_nuclear,
      chr_log2_norm = log2(chr_depth_norm)
    )
}

# =========================
# Collect files
# =========================
summary_files <- list.files(
  mos_dir,
  pattern = "\\.mosdepth\\.summary\\.txt$",
  full.names = TRUE
)

# =========================
# Parse all samples
# =========================
chr_depth_long <- purrr::map_dfr(summary_files, read_chr_depth_one) %>%
  mutate(
    sample_num = readr::parse_number(sample),
    sample_suffix = stringr::str_extract(sample, "[A-Za-z]+$"),
    chrom_num = readr::parse_number(chrom)
  ) %>%
  arrange(sample_num, sample_suffix, chrom_num, chrom) %>%
  select(-sample_num, -sample_suffix, -chrom_num)

# =========================
# Wide table
# =========================
chr_depth_wide <- chr_depth_long %>%
  select(sample, chrom, chr_depth_norm) %>%
  pivot_wider(
    names_from = chrom,
    values_from = chr_depth_norm
  ) %>%
  mutate(
    sample_num = readr::parse_number(sample),
    sample_suffix = stringr::str_extract(sample, "[A-Za-z]+$")
  ) %>%
  arrange(sample_num, sample_suffix, sample) %>%
  select(-sample_num, -sample_suffix)

# =========================
# Simple CNV flag (preliminary screening)
# =========================
chr_depth_flag <- chr_depth_long %>%
  mutate(
    chr_cnv_flag = case_when(
      is.na(chr_depth_norm) ~ NA_character_,
      chr_depth_norm >= 1.30 ~ "gain_candidate",
      chr_depth_norm <= 0.70 ~ "loss_candidate",
      TRUE ~ "neutral_like"
    )
  )

# chr abnormal in each sample
sample_chr_flag_summary <- chr_depth_flag %>%
  group_by(sample) %>%
  summarise(
    n_chr_gain_candidate = sum(chr_cnv_flag == "gain_candidate", na.rm = TRUE),
    n_chr_loss_candidate = sum(chr_cnv_flag == "loss_candidate", na.rm = TRUE),
    max_chr_depth_norm = max(chr_depth_norm, na.rm = TRUE),
    min_chr_depth_norm = min(chr_depth_norm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    chr_level_cnv_hint = case_when(
      n_chr_gain_candidate == 0 & n_chr_loss_candidate == 0 ~ "none",
      n_chr_gain_candidate + n_chr_loss_candidate == 1 ~ "single_chr_candidate",
      n_chr_gain_candidate + n_chr_loss_candidate >= 2 ~ "multi_chr_candidate",
      TRUE ~ "none"
    )
  ) %>%
  mutate(
    sample_num = readr::parse_number(sample),
    sample_suffix = stringr::str_extract(sample, "[A-Za-z]+$")
  ) %>%
  arrange(sample_num, sample_suffix, sample) %>%
  select(-sample_num, -sample_suffix)

# =========================
# Save
# =========================
write_tsv(chr_depth_long, out_long)
write_tsv(chr_depth_wide, out_wide)
write_tsv(
  sample_chr_flag_summary,
  file.path(out_dir, "mosdepth_chr_level_cnv_hint.tsv")
)

message("Wrote: ", out_long)
message("Wrote: ", out_wide)
message("Wrote: ", file.path(out_dir, "mosdepth_chr_level_cnv_hint.tsv"))

# =========================
# Quick inspection 4.1a
# =========================
chr_depth_long %>%
  group_by(chrom) %>%
  summarise(
    n = n(),
    median_norm = median(chr_depth_norm, na.rm = TRUE),
    q1 = quantile(chr_depth_norm, 0.25, na.rm = TRUE),
    q3 = quantile(chr_depth_norm, 0.75, na.rm = TRUE),
    min_norm = min(chr_depth_norm, na.rm = TRUE),
    max_norm = max(chr_depth_norm, na.rm = TRUE)
  ) %>%
  arrange(desc(median_norm))

sample_chr_flag_summary %>%
  arrange(desc(max_chr_depth_norm)) %>%
  print(n = 30)

sample_chr_flag_summary %>%
  filter(chr_level_cnv_hint != "none") %>%
  arrange(desc(n_chr_gain_candidate), desc(max_chr_depth_norm))

#############################################
# STEP 4.1b chromosome-centered robust table
#############################################

# ---- helper: robust z based on MAD ----
robust_z_from_median_mad <- function(x) {
  med <- median(x, na.rm = TRUE)
  mad_val <- mad(x, center = med, constant = 1, na.rm = TRUE)
  
  if (is.na(mad_val) || mad_val == 0) {
    return(rep(NA_real_, length(x)))
  }
  
  (x - med) / (1.4826 * mad_val)
}

# ---- chromosome-level baseline across cohort ----
chrom_baseline_summary <- chr_depth_long %>%
  group_by(chrom) %>%
  summarise(
    n = n(),
    chrom_median_norm = median(chr_depth_norm, na.rm = TRUE),
    chrom_mean_norm = mean(chr_depth_norm, na.rm = TRUE),
    chrom_sd_norm = sd(chr_depth_norm, na.rm = TRUE),
    chrom_mad_norm = mad(chr_depth_norm, center = median(chr_depth_norm, na.rm = TRUE), constant = 1, na.rm = TRUE),
    chrom_q1_norm = quantile(chr_depth_norm, 0.25, na.rm = TRUE),
    chrom_q3_norm = quantile(chr_depth_norm, 0.75, na.rm = TRUE),
    chrom_iqr_norm = IQR(chr_depth_norm, na.rm = TRUE),
    chrom_min_norm = min(chr_depth_norm, na.rm = TRUE),
    chrom_max_norm = max(chr_depth_norm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    chrom_num = readr::parse_number(chrom)
  ) %>%
  arrange(chrom_num, chrom) %>%
  select(-chrom_num)

# ---- join baseline back to each sample x chromosome ----
chr_depth_centered_long <- chr_depth_long %>%
  left_join(chrom_baseline_summary, by = "chrom") %>%
  group_by(chrom) %>%
  mutate(
    # relative to chromosome-wide median
    chrom_centered_ratio = chr_depth_norm / chrom_median_norm,
    chrom_centered_log2 = log2(chrom_centered_ratio),
    
    # robust z within each chromosome across all samples
    chrom_robust_z = robust_z_from_median_mad(chr_depth_norm),
    
    # IQR-based outlier rule
    chrom_iqr_upper = chrom_q3_norm + 1.5 * chrom_iqr_norm,
    chrom_iqr_lower = chrom_q1_norm - 1.5 * chrom_iqr_norm,
    chrom_iqr_outlier = case_when(
      is.na(chr_depth_norm) ~ NA,
      chr_depth_norm > chrom_iqr_upper ~ TRUE,
      chr_depth_norm < chrom_iqr_lower ~ TRUE,
      TRUE ~ FALSE
    ),
    
    # simple robust flags
    chrom_robust_flag = case_when(
      is.na(chrom_robust_z) ~ NA_character_,
      chrom_robust_z >= 3 ~ "gain_outlier",
      chrom_robust_z <= -3 ~ "loss_outlier",
      TRUE ~ "inlier"
    )
  ) %>%
  ungroup() %>%
  mutate(
    sample_num = readr::parse_number(sample),
    sample_suffix = stringr::str_extract(sample, "[A-Za-z]+$"),
    chrom_num = readr::parse_number(chrom)
  ) %>%
  arrange(sample_num, sample_suffix, chrom_num, chrom) %>%
  select(-sample_num, -sample_suffix, -chrom_num)

# ---- sample-level summary after chromosome-centering ----
sample_chr_centered_summary <- chr_depth_centered_long %>%
  group_by(sample) %>%
  summarise(
    n_gain_outlier = sum(chrom_robust_flag == "gain_outlier", na.rm = TRUE),
    n_loss_outlier = sum(chrom_robust_flag == "loss_outlier", na.rm = TRUE),
    max_chr_depth_norm = max(chr_depth_norm, na.rm = TRUE),
    min_chr_depth_norm = min(chr_depth_norm, na.rm = TRUE),
    max_centered_ratio = max(chrom_centered_ratio, na.rm = TRUE),
    min_centered_ratio = min(chrom_centered_ratio, na.rm = TRUE),
    max_robust_z = max(chrom_robust_z, na.rm = TRUE),
    min_robust_z = min(chrom_robust_z, na.rm = TRUE),
    n_iqr_outlier = sum(chrom_iqr_outlier %in% TRUE, na.rm = TRUE),
    centered_cnv_hint = case_when(
      n_gain_outlier == 0 & n_loss_outlier == 0 ~ "none",
      n_gain_outlier + n_loss_outlier == 1 ~ "single_chr_outlier",
      n_gain_outlier + n_loss_outlier >= 2 ~ "multi_chr_outlier",
      TRUE ~ "none"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    total_outlier = n_gain_outlier + n_loss_outlier,
    sample_num = readr::parse_number(sample),
    sample_suffix = stringr::str_extract(sample, "[A-Za-z]+$")
  ) %>%
  # rank the candidates
  arrange(
    desc(max_centered_ratio),
    desc(total_outlier),
    desc(max_robust_z),
    sample_num, sample_suffix, sample
  ) %>%
  mutate(
    candidate_priority_rank = row_number()
  ) %>%
  # reorder to natural seq
  arrange(sample_num, sample_suffix, sample) %>%
  select(-sample_num, -sample_suffix)

# ---- wide table of chromosome-centered ratio ----
chr_depth_centered_wide <- chr_depth_centered_long %>%
  select(sample, chrom, chrom_centered_ratio) %>%
  pivot_wider(
    names_from = chrom,
    values_from = chrom_centered_ratio
  ) %>%
  mutate(
    sample_num = readr::parse_number(sample),
    sample_suffix = stringr::str_extract(sample, "[A-Za-z]+$")
  ) %>%
  arrange(sample_num, sample_suffix, sample) %>%
  select(-sample_num, -sample_suffix)

# ---- save ----
write_tsv(
  chrom_baseline_summary,
  file.path(out_dir, "mosdepth_chrom_baseline_summary.tsv")
)

write_tsv(
  chr_depth_centered_long,
  file.path(out_dir, "mosdepth_chr_centered_depth_long.tsv")
)

write_tsv(
  chr_depth_centered_wide,
  file.path(out_dir, "mosdepth_chr_centered_depth_wide.tsv")
)

write_tsv(
  sample_chr_centered_summary,
  file.path(out_dir, "mosdepth_chr_centered_cnv_hint.tsv")
)

message("Wrote: ", file.path(out_dir, "mosdepth_chrom_baseline_summary.tsv"))
message("Wrote: ", file.path(out_dir, "mosdepth_chr_centered_depth_long.tsv"))
message("Wrote: ", file.path(out_dir, "mosdepth_chr_centered_depth_wide.tsv"))
message("Wrote: ", file.path(out_dir, "mosdepth_chr_centered_cnv_hint.tsv"))

# =========================
# Quick inspection 4.1b
# =========================

chrom_baseline_summary %>%
  print(n = Inf)

sample_chr_centered_summary %>%
  count(centered_cnv_hint) %>%
  print()

sample_chr_centered_summary %>%
  arrange(desc(candidate_priority_rank)) %>%
  print(n = 30)

chr_depth_centered_long %>%
  filter(chrom_robust_flag != "inlier") %>%
  select(
    sample, chrom, chr_depth_norm,
    chrom_median_norm, chrom_centered_ratio,
    chrom_robust_z, chrom_robust_flag
  ) %>%
  arrange(desc(abs(chrom_robust_z))) %>%
  print(n = 50)

############################################
### STEP 4.2 1 kb bin genome-wide coverage plot
############################################

# =========================
# Directory
# =========================
mos_dir <- MOS_DIR          # from wgs_analysis_config.R
out_dir <- DEPTH_DIR

plot_dir <- file.path(out_dir, "cnv_genomewide_plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Files from 4.1a
chr_depth_long_file <- file.path(out_dir, "mosdepth_chr_normalized_depth_long.tsv")
chr_depth_long <- read_tsv(chr_depth_long_file, show_col_types = FALSE)

# =============================
# Candidate samples list switch
# =============================
# strategy:
candidate_mode <- "all"    # "all" | "manual" | "top_n"
top_n_candidate <- 15
candidate_samples <- switch(candidate_mode,
                            "all" = sample_chr_centered_summary %>%
                              filter(centered_cnv_hint != "none") %>%
                              arrange(candidate_priority_rank) %>%
                              pull(sample),
                            
                            "manual" = c("Cn147", "Cn202", "Cn81"),
                            
                            "top_n" = sample_chr_centered_summary %>%
                              arrange(candidate_priority_rank) %>%
                              slice_head(n = top_n_candidate) %>%
                              pull(sample))
message("Candidate mode: ", candidate_mode,
        " → ", length(candidate_samples), " samples selected")

# =========================
# FLC-related loci
# =========================
flc_loci <- tibble(
  gene  = c("ERG11", "AFR1"),
  chrom = c("NC_026745.1", "NC_026745.1"),
  start = c(123979, 1900719),
  end   = c(126449, 1907783)
)

# =========================
# Read one regions.bed.gz
# =========================
read_regions_one <- function(file, chr_depth_long_tbl) {
  sample <- sub("\\.regions\\.bed\\.gz$", "", basename(file))
  
  x <- read_tsv(
    file,
    col_names = c("chrom", "start", "end", "bin_mean_depth"),
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    mutate(
      start = as.numeric(start),
      end = as.numeric(end),
      bin_mean_depth = as.numeric(bin_mean_depth),
      bin_width = end - start
    ) %>%
    filter(!is.na(start), !is.na(end), !is.na(bin_mean_depth), bin_width > 0)
  
  # Normalization using sample_mean_depth_nuclear in 4.1a
  sample_mean_depth <- chr_depth_long_tbl %>%
    filter(sample == !!sample) %>%
    distinct(sample, sample_mean_depth_nuclear) %>%
    pull(sample_mean_depth_nuclear)
  
  if (length(sample_mean_depth) == 0 || is.na(sample_mean_depth[1])) {
    sample_mean_depth <- NA_real_
  } else {
    sample_mean_depth <- sample_mean_depth[1]
  }
  
  x %>%
    mutate(
      sample = sample,
      sample_mean_depth_nuclear = sample_mean_depth,
      bin_depth_norm = bin_mean_depth / sample_mean_depth_nuclear,
      bin_log2_norm = log2(bin_depth_norm),
      bin_mid = (start + end) / 2
    ) %>%
    select(
      sample, chrom, start, end, bin_width, bin_mid,
      bin_mean_depth, sample_mean_depth_nuclear,
      bin_depth_norm, bin_log2_norm
    )
}

# =========================
# Collect region files
# =========================
region_files <- list.files(
  mos_dir,
  pattern = "\\.regions\\.bed\\.gz$",
  full.names = TRUE
)

# Only for candidates
region_files_use <- region_files[
  sub("\\.regions\\.bed\\.gz$", "", basename(region_files)) %in% candidate_samples
]

# =========================
# Parse candidate samples
# =========================
regions_long <- purrr::map_dfr(region_files_use, read_regions_one, chr_depth_long_tbl = chr_depth_long)
write_tsv(regions_long, file.path(out_dir, "candidate_regions_normalized_long.tsv")) # Long pivot

# =========================
# Build genome coordinate system
# 目的：把每条染色体拼起来，画成一条 genome-wide x 轴
# 用 chr_depth_long 里的 chr_length 作为 reference
# =========================
chrom_ref <- chr_depth_long %>%
  distinct(chrom, chr_length) %>%
  mutate(
    chrom_num = readr::parse_number(chrom)
  ) %>%
  arrange(chrom_num, chrom) %>%
  mutate(
    chrom_offset = lag(cumsum(chr_length), default = 0),
    chrom_center = chrom_offset + chr_length / 2
  ) %>%
  select(-chrom_num)

regions_plot_df <- regions_long %>%
  left_join(chrom_ref, by = "chrom") %>%
  mutate(
    genome_x = chrom_offset + bin_mid
  )

# =========================
# Quick summaries
# =========================
sample_region_summary <- regions_plot_df %>%
  group_by(sample) %>%
  summarise(
    n_bins = n(),
    median_bin_norm = median(bin_depth_norm, na.rm = TRUE),
    min_bin_norm = min(bin_depth_norm, na.rm = TRUE),
    max_bin_norm = max(bin_depth_norm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(max_bin_norm))

write_tsv(sample_region_summary, file.path(out_dir, "candidate_regions_summary.tsv"))

# =========================
# Plot function
# =========================
plot_genomewide_bins <- function(
    df_sample,
    chrom_ref_tbl,
    flc_tbl = NULL,
    mode = c("clipped", "raw", "log2"),
    clip_upper = 3
) {
  mode <- match.arg(mode)
  sample_name <- unique(df_sample$sample)
  
  df_plot <- df_sample
  
  if (mode == "raw") {
    df_plot <- df_plot %>%
      mutate(y_plot = bin_depth_norm)
    y_lab <- "1 kb bin normalized depth"
    hline_main <- 1
    hline_aux <- c(0.7, 1.3)
  }
  
  if (mode == "clipped") {
    df_plot <- df_plot %>%
      mutate(y_plot = pmin(bin_depth_norm, clip_upper))
    y_lab <- paste0("1 kb bin normalized depth (clipped at ", clip_upper, ")")
    hline_main <- 1
    hline_aux <- c(0.7, 1.3)
  }
  
  if (mode == "log2") {
    df_plot <- df_plot %>%
      mutate(y_plot = log2(bin_depth_norm))
    y_lab <- "log2(1 kb bin normalized depth)"
    hline_main <- 0
    hline_aux <- c(log2(0.7), log2(1.3))
  }
  
  p <- ggplot(df_plot, aes(x = genome_x, y = y_plot)) +
    geom_point(size = 0.25, alpha = 0.6) +
    geom_hline(yintercept = hline_main, linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = hline_aux, linetype = "dotted", linewidth = 0.3) +
    geom_vline(
      data = chrom_ref_tbl %>% filter(chrom_offset > 0),
      aes(xintercept = chrom_offset),
      inherit.aes = FALSE,
      linewidth = 0.25,
      alpha = 0.6
    ) +
    scale_x_continuous(
      breaks = chrom_ref_tbl$chrom_center,
      labels = chrom_ref_tbl$chrom
    ) +
    labs(
      title = paste0(sample_name, " (", mode, ")"),
      x = "Genome position (chromosomes concatenated)",
      y = y_lab
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)
    )
  
  if (!is.null(flc_tbl) && nrow(flc_tbl) > 0) {
    flc_ok <- flc_tbl %>%
      filter(!is.na(chrom), !is.na(start), !is.na(end)) %>%
      left_join(chrom_ref_tbl, by = "chrom") %>%
      mutate(
        flc_xmin = chrom_offset + start,
        flc_xmax = chrom_offset + end,
        flc_xmid = (flc_xmin + flc_xmax) / 2
      )
    
    if (nrow(flc_ok) > 0) {
      ymax_now <- max(df_plot$y_plot, na.rm = TRUE)
      p <- p +
        geom_rect(
          data = flc_ok,
          aes(xmin = flc_xmin, xmax = flc_xmax, ymin = -Inf, ymax = Inf),
          inherit.aes = FALSE,
          alpha = 0.08
        ) +
        geom_text(
          data = flc_ok,
          aes(x = flc_xmid, y = ymax_now, label = gene),
          inherit.aes = FALSE,
          vjust = -0.4,
          size = 3
        )
    }
  }
  
  p
}

# =========================
# Plot candidate samples
# =========================
plots <- regions_plot_df %>%
  group_split(sample)

for (df_s in plots) {
  s <- unique(df_s$sample)
  
  p1 <- plot_genomewide_bins(df_s, chrom_ref, flc_loci, mode = "clipped", clip_upper = 3)
  p2 <- plot_genomewide_bins(df_s, chrom_ref, flc_loci, mode = "log2")
  
  ggsave(
    filename = file.path(plot_dir, paste0(s, "_genomewide_clipped.png")),
    plot = p1, width = 12, height = 4.8, dpi = 300
  )
  
  ggsave(
    filename = file.path(plot_dir, paste0(s, "_genomewide_log2.png")),
    plot = p2, width = 12, height = 4.8, dpi = 300
  )
}

message("Wrote: ", file.path(out_dir, "candidate_regions_normalized_long.tsv"))
message("Wrote: ", file.path(out_dir, "candidate_regions_summary.tsv"))
message("Plots saved in: ", plot_dir)

# =========================
# Plot using karyoploteR
# =========================
library(karyoploteR)
library(GenomicRanges)
library(regioneR)

# karyoploteR custom genome
custom_genome <- toGRanges(data.frame(
  chr = chrom_ref$chrom,
  start = 1,
  end = chrom_ref$chr_length
))

# regions to GRanges
make_bin_granges <- function(df_sample) {
  gr <- GRanges(
    seqnames = df_sample$chrom,
    ranges = IRanges(start = df_sample$start + 1, end = df_sample$end)
  )
  
  mcols(gr)$bin_depth_norm <- df_sample$bin_depth_norm
  mcols(gr)$bin_log2_norm <- df_sample$bin_log2_norm
  mcols(gr)$bin_depth_clipped <- pmin(df_sample$bin_depth_norm, 3)
  
  return(gr)
}
make_flc_granges <- function(flc_tbl) {
  flc_ok <- flc_tbl %>%
    filter(!is.na(chrom), !is.na(start), !is.na(end))
  
  if (nrow(flc_ok) == 0) return(NULL)
  
  gr <- GRanges(
    seqnames = flc_ok$chrom,
    ranges = IRanges(start = flc_ok$start, end = flc_ok$end)
  )
  mcols(gr)$gene <- flc_ok$gene
  gr
}
flc_gr <- make_flc_granges(flc_loci)

# Plot using karyoplateR
plot_genomewide_karyo <- function(
    df_sample,
    genome_gr,
    flc_gr = NULL,
    mode = c("clipped", "raw", "log2"),
    clip_upper = 3,
    main_cex = 1
) {
  mode <- match.arg(mode)
  sample_name <- unique(df_sample$sample)
  
  stopifnot(length(sample_name) == 1)
  
  gr <- make_bin_granges(df_sample)
  
  if (mode == "raw") {
    y <- mcols(gr)$bin_depth_norm
    ymin <- 0
    ymax <- max(y, na.rm = TRUE)
    main_line <- 1
    aux_lines <- c(0.7, 1.3)
    ylab <- "Normalized depth"
  }
  
  if (mode == "clipped") {
    y <- mcols(gr)$bin_depth_clipped
    ymin <- 0
    ymax <- clip_upper
    main_line <- 1
    aux_lines <- c(0.7, 1.3)
    ylab <- paste0("Normalized depth (clipped at ", clip_upper, ")")
  }
  
  if (mode == "log2") {
    y <- mcols(gr)$bin_log2_norm
    ymin <- min(y[is.finite(y)], na.rm = TRUE)
    ymax <- max(y[is.finite(y)], na.rm = TRUE)
    main_line <- 0
    aux_lines <- c(log2(0.7), log2(1.3))
    ylab <- "log2(normalized depth)"
  }
  
  kp <- plotKaryotype(
    genome = genome_gr,
    plot.type = 1,
    main = paste0(sample_name, " (", mode, ")"),
    cex = main_cex, chromosomes = "all"
  )
  
  # Background and Axis
  kpDataBackground(kp, data.panel = 1)
  kpAxis(kp, ymin = ymin, ymax = ymax, data.panel = 1, cex = 0.7)
  
  # kpAbline
  kpAbline(kp, h = main_line, ymin = ymin, ymax = ymax, data.panel = 1)
  for (h in aux_lines) {
    kpAbline(kp, h = h, ymin = ymin, ymax = ymax, data.panel = 1, lty = 3)
  }
  
  # kpPoints
  kpPoints(
    kp,
    data = gr,
    y = y,
    ymin = ymin,
    ymax = ymax,
    cex = 0.2,
    pch = 16,
    data.panel = 1
  )
  
  # FLC
  if (!is.null(flc_gr) && length(flc_gr) > 0) {
    kpPlotRegions(
      kp,
      data = flc_gr,
      data.panel = 1,
      r0 = 0,
      r1 = 1,
      col = "#d95f0215",
      border = "#d95f02"
    )
    
    kpText(
      kp,
      data = flc_gr,
      labels = mcols(flc_gr)$gene,
      y = rep(ymax, length(flc_gr)),
      ymin = ymin,
      ymax = ymax,
      pos = 3,
      cex = 0.7,
      data.panel = 1
    )
  }
  
  invisible(kp)
}

candidate_df_list <- regions_plot_df %>%
  group_split(sample)
plot_dir_karyo <- file.path(out_dir, "cnv_genomewide_plots_karyo")
dir.create(plot_dir_karyo, showWarnings = FALSE, recursive = TRUE)
for (df_s in candidate_df_list) {
  s <- unique(df_s$sample)
  
  png(file.path(plot_dir_karyo, paste0(s, "_karyo_clipped.png")),
      width = 4000, height = 2800, res = 250)
  plot_genomewide_karyo(
    df_sample = df_s,
    genome_gr = custom_genome,
    flc_gr = flc_gr,
    mode = "clipped",
    clip_upper = 3
  )
  dev.off()
  
  png(file.path(plot_dir_karyo, paste0(s, "_karyo_log2.png")),
      width = 4000, height = 2800, res = 250)
  plot_genomewide_karyo(
    df_sample = df_s,
    genome_gr = custom_genome,
    flc_gr = flc_gr,
    mode = "log2"
  )
  dev.off()
}

plot_genomewide_karyo(df_sample = candidate_df_list[[7]], genome_gr = custom_genome, 
                      flc_gr = flc_gr, mode = "clipped", clip_upper = 3)

# =========================
# Quick inspection 4.2
# =========================
sample_region_summary %>%
  print(n = 30)

regions_plot_df %>%
  select(sample, chrom, start, end, bin_mean_depth, sample_mean_depth_nuclear, bin_depth_norm) %>%
  print(n = 30)

############################################
### STEP 4.3 CNV candidate review
############################################

# ---- 候选样本优先级来源 ----
# 这里优先用 4.1b 的 centered summary
# 原因：它已经校正了 chromosome-specific baseline
candidate_priority_tbl <- sample_chr_centered_summary %>%
  mutate(
    total_outlier = n_gain_outlier + n_loss_outlier
  ) %>%
  arrange(
    desc(max_centered_ratio),
    desc(total_outlier),
    desc(max_robust_z)
  )

# =========================
# 1) sample x chrom 层面的自动证据汇总
# =========================
# =========================
# 1.1) chromosome-level 证据 (来自 4.1b, 每 sample×chrom 一行)
# =========================
sample_chrom_event_chr_tbl <- chr_depth_centered_long %>%
  filter(sample %in% candidate_samples) %>%
  transmute(
    sample, chrom,
    chr_depth_norm,
    chrom_centered_ratio,
    chrom_robust_z,       # 单值, 不需要 max/min
    chrom_robust_flag
  )

# =========================
# 1.2) 1 kb bin-level 证据 (来自 4.2)
# 新增: sd 和 outlier bin count, 用于区分 uniform vs focal
# =========================
sample_chrom_event_bin_tbl <- regions_plot_df %>%
  filter(sample %in% candidate_samples) %>%
  group_by(sample, chrom) %>%
  summarise(
    n_bins                  = n(),
    median_bin_norm         = median(bin_depth_norm, na.rm = TRUE),
    max_bin_norm            = max(bin_depth_norm, na.rm = TRUE),
    min_bin_norm            = min(bin_depth_norm, na.rm = TRUE),
    sd_bin_norm             = sd(bin_depth_norm, na.rm = TRUE),
    # 高于 1.5 的 bin 占比 → focal vs whole-chr
    frac_bins_above_1.5     = sum(bin_depth_norm > 1.5, na.rm = TRUE) / n(),
    # 低于 0.5 的 bin 占比
    frac_bins_below_0.5     = sum(bin_depth_norm < 0.5, na.rm = TRUE) / n(),
    .groups = "drop"
  )

# =========================
# 1.3) 合并两层证据
# =========================
sample_chrom_event_tbl <- sample_chrom_event_chr_tbl %>%
  left_join(sample_chrom_event_bin_tbl, by = c("sample", "chrom")) %>%
  mutate(
    event_score = abs(chrom_robust_z),
    auto_pattern_hint = case_when(
      !is.na(sd_bin_norm) & sd_bin_norm >= 1.5 ~ "high_variance_artifact",  # mappability 伪迹 (chr2 类)
      chr_depth_norm >= 1.5 & sd_bin_norm <  0.5 ~ "uniform_gain",          # 0.3→0.5: 真 disomy sd~0.4
      chr_depth_norm >= 1.5 & sd_bin_norm >= 0.5 ~ "focal_or_noisy_gain",
      chr_depth_norm <= 0.7 & sd_bin_norm <  0.5 ~ "uniform_loss",
      chr_depth_norm <= 0.7 & sd_bin_norm >= 0.5 ~ "focal_or_noisy_loss",
      frac_bins_above_1.5 > 0.05 & chr_depth_norm < 1.5 ~ "focal_gain_subregion",
      TRUE ~ "no_clear_pattern"
    )
  )

# FLC 相关染色体
flc_chrom <- unique(flc_loci$chrom)  # NC_026745.1

# =========================
# per-sample shadow 地板 (稳健版, 触发与地板值分离)
# 触发: 是否有真 disomy → centered ratio >=1.5 (排除 chr2 伪迹 ~1.1 和 partial gain)
# 地板值: copy-1 染色体的稳健中心 → median, 自动吸收 disomy + partial gain 的分母抬升
# =========================
sample_floor_summary <- chr_depth_centered_long %>%
  group_by(sample) %>%
  summarise(
    n_disomy = sum(chrom_centered_ratio >= 1.50 & chr_depth_norm >= 1.40, na.rm = TRUE),
    expected_floor_raw = median(
      chr_depth_norm[chr_depth_norm >= 0.50 & chr_depth_norm < 1.15],
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    # 仅当样本确有真 disomy 且地板明显 <1 才用作 shadow 预测
    expected_floor = if_else(n_disomy >= 1 & expected_floor_raw < 0.97,
                             expected_floor_raw, NA_real_)
  ) %>%
  select(sample, n_disomy, expected_floor)

# ---- shadow 参数 (可调) ----
shadow_tol      <- 0.05
shadow_lo_guard <- 0.78
shadow_hi_guard <- 0.97
shadow_sd_guard <- 0.50   # 0.3→0.5

cnv_review_template <- sample_chrom_event_tbl %>%
  left_join(
    sample_chr_centered_summary %>% select(sample, candidate_priority_rank),
    by = "sample"
  ) %>%
  left_join(sample_floor_summary, by = "sample") %>%
  filter(abs(chrom_centered_ratio - 1) > 0.15 | abs(event_score) > 2) %>%
  group_by(sample) %>%
  arrange(desc(event_score), desc(abs(chrom_centered_ratio - 1)), desc(sd_bin_norm), .by_group = TRUE) %>%
  slice_head(n = 3) %>%
  mutate(candidate_rank = row_number()) %>%
  ungroup() %>%
  mutate(
    passenger_shadow = case_when(
      is.na(chr_depth_norm) ~ NA,
      is.na(expected_floor) ~ FALSE,
      chr_depth_norm < shadow_hi_guard &
        chr_depth_norm > shadow_lo_guard &
        (is.na(sd_bin_norm) | sd_bin_norm < shadow_sd_guard) &
        abs(chr_depth_norm - expected_floor) <= shadow_tol ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  transmute(
    sample               = sample,
    sample_priority_rank = candidate_priority_rank,
    candidate_rank       = candidate_rank,
    candidate_chrom      = chrom,
    chr_depth_norm       = round(chr_depth_norm, 3),
    chrom_centered_ratio = round(chrom_centered_ratio, 3),
    chrom_robust_z       = round(chrom_robust_z, 2),   # 带符号
    sd_bin_norm          = round(sd_bin_norm, 3),
    auto_pattern_hint    = auto_pattern_hint,
    expected_floor       = round(expected_floor, 3),
    passenger_shadow     = passenger_shadow,
    pattern_type = case_when(
      chr_depth_norm >= 1.50 ~ "whole_chr_gain",
      chr_depth_norm <= 0.70 ~ "whole_chr_loss",
      chr_depth_norm >= 1.15 ~ "mild_gain_review",
      chr_depth_norm <= 0.85 ~ "mild_loss_review",
      TRUE ~ "neutral"
    ),
    pattern_confidence = case_when(
      passenger_shadow %in% TRUE ~ "auto_shadow",
      (chrom %in% flc_chrom) & chr_depth_norm > 0.70 & chr_depth_norm < 1.50 ~ "needs_manual_review",  # FLC 永远复核
      chr_depth_norm >= 1.50 | chr_depth_norm <= 0.70 ~ "high_auto",
      chr_depth_norm > 0.90 & chr_depth_norm < 1.10 &
        auto_pattern_hint == "no_clear_pattern" ~ "auto_neutral",
      TRUE ~ "needs_manual_review"
    ),
    flc_related = if_else(chrom %in% flc_chrom, "yes", "no"),
    flc_gene    = if_else(chrom %in% flc_chrom, "ERG11/AFR1", NA_character_),
    bin_pattern         = NA_character_,
    event_span          = NA_character_,
    supporting_evidence = NA_character_,
    comment = case_when(
      passenger_shadow %in% TRUE ~ "denominator shadow of in-sample disomy",
      auto_pattern_hint == "high_variance_artifact" ~ "high bin variance, likely mappability artifact; verify on plot",
      chr_depth_norm > 0.90 & chr_depth_norm < 1.10 &
        auto_pattern_hint == "no_clear_pattern" & !(chrom %in% flc_chrom) ~ "within-noise, no focal signal",
      TRUE ~ NA_character_
    ),
    reviewer    = NA_character_,
    review_date = NA_character_
  ) %>%
  mutate(sample_num = readr::parse_number(sample),
         sample_suffix = stringr::str_extract(sample, "[A-Za-z]+$")) %>%
  arrange(sample_priority_rank, sample_num, sample_suffix, candidate_rank) %>%
  select(-sample_num, -sample_suffix)

# ---- statistic ----
message("Template rows: ", nrow(cnv_review_template),
        " (from ", n_distinct(cnv_review_template$sample), " samples)")
message("  Auto-shadow (skip):  ", sum(cnv_review_template$passenger_shadow %in% TRUE))
message("  Auto-neutral (skip): ", sum(cnv_review_template$pattern_confidence == "auto_neutral"))
message("  High-confidence g/l: ", sum(cnv_review_template$pattern_confidence == "high_auto"))
message("  Needs manual review: ", sum(cnv_review_template$pattern_confidence == "needs_manual_review"))
message("  FLC-related:         ", sum(cnv_review_template$flc_related == "yes"))

# =========================
# 3) Save
# =========================
write_tsv(cnv_review_template, file.path(out_dir, "cnv_event_review_template.tsv"))
# =========================
# 4) Quick inspection 4.3
# =========================
candidate_priority_tbl %>%
  select(sample, n_gain_outlier, n_loss_outlier, max_centered_ratio, max_robust_z, centered_cnv_hint) %>%
  print(n = 30)

cnv_review_template %>%
  print(n = 50)

############################################
### STEP 5 Het SNP analysis — mixed infection & disomy detection
### 单倍体背景下, het SNP 的全基因组分布模式区分:
###   多克隆混合感染 vs 单克隆 disomy
############################################

# =========================
# Directory
# =========================
vcf_dir  <- VCF_FILTERED_DIR    # from wgs_analysis_config.R
out_dir  <- dirname(vcf_dir)
out_dir5 <- file.path(out_dir, "step5_het_analysis")
plot_dir5 <- file.path(out_dir5, "plots")
dir.create(out_dir5, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir5, showWarnings = FALSE, recursive = TRUE)

# Step 4 的输出目录 (cross-reference 用)
depth_dir <- DEPTH_DIR          # from wgs_analysis_config.R

# =========================
# Read pre-extracted het SNPs (After extract Het SNP from batch)
# =========================
het_extracted_file <- HET_FILE  # from wgs_analysis_config.R

het_long <- read_tsv(
  het_extracted_file,
  col_types = cols(
    sample    = col_character(),
    chrom     = col_character(),
    pos       = col_integer(),
    ref       = col_character(),
    alt       = col_character(),
    qual      = col_double(),
    filter_st = col_character(),
    gt        = col_character(),
    dp        = col_integer(),
    ad        = col_character()
  )
) %>%
  mutate(
    ad_ref = as.numeric(str_extract(ad, "^[0-9]+")),
    ad_alt = as.numeric(str_extract(ad, "(?<=,)[0-9]+")),
    af     = ad_alt / (ad_ref + ad_alt)
  ) %>%
  filter(!is.na(af), is.finite(af)) %>%
  select(sample, chrom, pos, ref, alt, qual, filter_st, gt, dp, ad_ref, ad_alt, af) %>%
  mutate(
    sample_num    = readr::parse_number(sample),
    sample_suffix = stringr::str_extract(sample, "[A-Za-z]+$"),
    chrom_num     = readr::parse_number(chrom)
  ) %>%
  arrange(sample_num, sample_suffix, chrom_num, chrom, pos) %>%
  select(-sample_num, -sample_suffix, -chrom_num)

message("Total het SNPs loaded: ", nrow(het_long))

# Easy being killed
# library(progressr)
# library(furrr)
# bcf_path <- "~/.conda/envs/bwasam_snake/bin/bcftools"
# options(future.globals.maxSize = 200 * 1024^3)
# plan(multisession, workers = 8)
# 
# extract_het_bcftools <- function(vcf_file, bcftools) {
#   sample <- sub("\\.filtered\\.vcf\\.gz$", "", basename(vcf_file))
#   
#   cmd <- paste0(
#     bcftools, " query ",
#     "-i 'GT=\"0/1\" && TYPE=\"snp\" && N_ALT=1' ",
#     "-f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t[%GT]\\t[%DP]\\t[%AD]\\n' ",
#     shQuote(vcf_file), " 2>/dev/null"
#   )
#   
#   raw <- system(cmd, intern = TRUE)
#   if (length(raw) == 0) {
#     return(tibble(
#       sample = character(), chrom = character(), pos = numeric(),
#       ref = character(), alt = character(), qual = numeric(),
#       filter_st = character(), gt = character(), dp = numeric(),
#       ad_ref = numeric(), ad_alt = numeric(), af = numeric()
#     ))
#   }
#   
#   read_tsv(
#     I(raw),
#     col_names = c("chrom", "pos", "ref", "alt", "qual", "filter_st", "gt", "dp", "ad"),
#     show_col_types = FALSE,
#     progress = FALSE
#   ) %>%
#     mutate(
#       sample = sample,
#       pos    = as.numeric(pos),
#       qual   = suppressWarnings(as.numeric(qual)),
#       dp     = as.numeric(dp),
#       ad_ref = as.numeric(str_extract(ad, "^[0-9]+")),
#       ad_alt = as.numeric(str_extract(ad, "(?<=,)[0-9]+")),
#       af     = ad_alt / (ad_ref + ad_alt)
#     ) %>%
#     filter(!is.na(af), is.finite(af)) %>%
#     select(sample, chrom, pos, ref, alt, qual, filter_st, gt, dp, ad_ref, ad_alt, af)
# }
# 
# vcf_dir is VCF_FILTERED_DIR from wgs_analysis_config.R
# vcf_files <- list.files(vcf_dir, pattern = "\\.filtered\\.vcf\\.gz$", full.names = TRUE)
# message("Found ", length(vcf_files), " filtered VCF files")
# 
# handlers(handler_txtprogressbar(
#   style = 3,
#   width = 60
# ))
# 
# het_long <- with_progress({
#   p <- progressor(along = vcf_files)
#   
#   future_map_dfr(vcf_files, function(f) {
#     result <- extract_het_bcftools(f, bcftools = bcf_path)
#     p()
#     result
#   }, .options = furrr_options(seed = TRUE))
# })
# 
# het_long <- het_long %>%
#   mutate(
#     sample_num    = readr::parse_number(sample),
#     sample_suffix = stringr::str_extract(sample, "[A-Za-z]+$"),
#     chrom_num     = readr::parse_number(chrom)
#   ) %>%
#   arrange(sample_num, sample_suffix, chrom_num, chrom, pos) %>%
#   select(-sample_num, -sample_suffix, -chrom_num)
# 
# message("Total het SNPs loaded: ", nrow(het_long))
# 
# plan(sequential)

# Too slow, deprecated
# het_long <- purrr::map_dfr(vcf_files, read_het_snp_one) %>%
#   mutate(
#     sample_num    = readr::parse_number(sample),
#     sample_suffix = stringr::str_extract(sample, "[A-Za-z]+$"),
#     chrom_num     = readr::parse_number(chrom)
#   ) %>%
#   arrange(sample_num, sample_suffix, chrom_num, chrom, pos) %>%
#   select(-sample_num, -sample_suffix, -chrom_num)
# 
# message("Total het SNPs parsed: ", nrow(het_long))

# =========================
# DP + AF filter (extra layer on top of bcftools filter)
# =========================
dp_min <- 10
dp_max_quantile <- 0.995
af_min <- 0.15
af_max <- 0.85

dp_upper <- het_long %>%
  pull(dp) %>%
  quantile(dp_max_quantile, na.rm = TRUE)

het_filt <- het_long %>%
  filter(dp >= dp_min, dp <= dp_upper) %>%
  filter(af >= af_min, af <= af_max)

message("Het SNPs after DP filter: ",
        sum(het_long$dp >= dp_min & het_long$dp <= dp_upper),
        " → after AF filter [", af_min, ",", af_max, "]: ", nrow(het_filt))

############################################
### 5.1 Genome-wide het SNP density per sample
### 修正: 加入 cohort baseline, 用相对偏移判断
############################################

# Chromosome lengths from Step 4
chr_depth_long_file <- file.path(depth_dir, "mosdepth_chr_normalized_depth_long.tsv")
chr_depth_long_4 <- read_tsv(chr_depth_long_file, show_col_types = FALSE)

chrom_lengths <- chr_depth_long_4 %>%
  distinct(chrom, chr_length)

genome_size <- sum(chrom_lengths$chr_length, na.rm = TRUE)

# 全部样本列表 (从 Step 4 获取, 保证不丢样本)
all_samples <- chr_depth_long_4 %>%
  distinct(sample) %>%
  pull(sample)

# 全部 chromosome 列表
all_chroms <- chrom_lengths %>% pull(chrom)
n_chrom_total <- length(all_chroms)

# per-sample genome-wide summary (先只算有 het 的样本)
het_genome_raw <- het_filt %>%
  group_by(sample) %>%
  summarise(
    total_het        = n(),
    median_af        = median(af, na.rm = TRUE),
    mean_af          = mean(af, na.rm = TRUE),
    sd_af            = sd(af, na.rm = TRUE),
    median_dp        = median(dp, na.rm = TRUE),
    n_chrom_with_het = n_distinct(chrom),
    .groups = "drop"
  )

# 补全缺失样本 (0 het 的样本)
het_genome_summary <- tibble(sample = all_samples) %>%
  left_join(het_genome_raw, by = "sample") %>%
  mutate(
    across(c(total_het, n_chrom_with_het), ~ replace_na(.x, 0L)),
    het_per_kb = total_het / (genome_size / 1000)
  )

# cohort baseline: 中位数和 MAD
cohort_median_het <- median(het_genome_summary$total_het)
cohort_mad_het    <- mad(het_genome_summary$total_het, constant = 1.4826)

message("Cohort baseline: median het = ", cohort_median_het,
        ", MAD = ", round(cohort_mad_het, 1))

het_genome_summary <- het_genome_summary %>%
  mutate(
    het_robust_z = (total_het - cohort_median_het) / cohort_mad_het,
    genome_het_flag = case_when(
      total_het == 0                   ~ "no_het",
      het_robust_z <= 2                ~ "baseline",      # cohort 背景水平
      het_robust_z > 2 & het_robust_z <= 5 ~ "elevated",  # 明显偏高
      het_robust_z > 5                 ~ "high",           # 强烈偏高
      TRUE                             ~ "baseline"
    )
  ) %>%
  mutate(
    sample_num    = readr::parse_number(sample),
    sample_suffix = stringr::str_extract(sample, "[A-Za-z]+$")
  ) %>%
  arrange(sample_num, sample_suffix) %>%
  select(-sample_num, -sample_suffix)

############################################
### 5.2 Per-chromosome het SNP density & uniformity
### 修正: 分母用全基因组 chr 数, 不是有 het 的 chr 数
############################################

het_chrom_summary <- het_filt %>%
  group_by(sample, chrom) %>%
  summarise(
    n_het     = n(),
    median_af = median(af, na.rm = TRUE),
    mean_af   = mean(af, na.rm = TRUE),
    sd_af     = sd(af, na.rm = TRUE),
    q25_af    = quantile(af, 0.25, na.rm = TRUE),
    q75_af    = quantile(af, 0.75, na.rm = TRUE),
    median_dp = median(dp, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  left_join(chrom_lengths, by = "chrom") %>%
  mutate(
    het_per_kb_chrom = n_het / (chr_length / 1000)
  )

# 每条 chr 的 cohort baseline
chrom_het_baseline <- het_chrom_summary %>%
  group_by(chrom) %>%
  summarise(
    chrom_median_het     = median(n_het, na.rm = TRUE),
    chrom_mad_het        = mad(n_het, constant = 1.4826, na.rm = TRUE),
    chrom_median_het_kb  = median(het_per_kb_chrom, na.rm = TRUE),
    .groups = "drop"
  )

# 每个 sample x chrom 相对于 baseline 的 z
het_chrom_summary <- het_chrom_summary %>%
  left_join(chrom_het_baseline, by = "chrom") %>%
  mutate(
    chrom_het_z = if_else(
      chrom_mad_het > 0,
      (n_het - chrom_median_het) / chrom_mad_het,
      NA_real_
    ),
    chrom_het_flag = case_when(
      is.na(chrom_het_z)    ~ NA_character_,
      chrom_het_z > 3       ~ "het_enriched",
      chrom_het_z < -3      ~ "het_depleted",
      TRUE                  ~ "baseline"
    )
  )

# 均匀度 (分母修正为全基因组 chr 数)
het_uniformity <- het_chrom_summary %>%
  group_by(sample) %>%
  summarise(
    n_chrom_with_het     = n(),
    total_het            = sum(n_het),
    max_het_on_one_chrom = max(n_het),
    chrom_with_max_het   = chrom[which.max(n_het)],
    max_chrom_het_frac   = max(n_het) / sum(n_het),
    top2_chrom_het_frac  = sum(sort(n_het, decreasing = TRUE)[1:min(2, n())]) / sum(n_het),
    cv_het_per_kb        = sd(het_per_kb_chrom, na.rm = TRUE) /
      mean(het_per_kb_chrom, na.rm = TRUE),
    # 修正: 分母用全基因组 chr 数
    frac_chrom_with_het  = n_distinct(chrom[n_het >= 5]) / n_chrom_total,
    # 新增: 有多少条 chr 的 het 显著高于 baseline
    n_chrom_het_enriched = sum(chrom_het_flag == "het_enriched", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # 补全 0 het 样本
  right_join(tibble(sample = all_samples), by = "sample") %>%
  mutate(
    across(c(n_chrom_with_het, total_het, max_het_on_one_chrom,
             n_chrom_het_enriched), ~ replace_na(.x, 0L)),
    across(c(max_chrom_het_frac, top2_chrom_het_frac, cv_het_per_kb,
             frac_chrom_with_het), ~ replace_na(.x, 0)),
    chrom_with_max_het = replace_na(chrom_with_max_het, NA_character_)
  ) %>%
  mutate(
    het_distribution = case_when(
      total_het == 0    ~ "no_het",
      total_het <= 20   ~ "too_few",
      # het 高度集中在 1-2 条 chr
      max_chrom_het_frac >= 0.60 ~ "concentrated",
      top2_chrom_het_frac >= 0.75 ~ "concentrated",
      # het 分散但有个别 chr 显著富集 → 混合模式
      n_chrom_het_enriched >= 1 & frac_chrom_with_het >= 0.5 ~ "dispersed_with_hotspot",
      # het 均匀分散
      frac_chrom_with_het >= 0.5 & cv_het_per_kb < 1.5 ~ "dispersed_uniform",
      TRUE ~ "intermediate"
    )
  ) %>%
  mutate(
    sample_num    = readr::parse_number(sample),
    sample_suffix = stringr::str_extract(sample, "[A-Za-z]+$")
  ) %>%
  arrange(sample_num, sample_suffix) %>%
  select(-sample_num, -sample_suffix)

############################################
### 5.3 AF distribution shape (不变, 但补全样本)
############################################

af_shape_raw <- het_filt %>%
  group_by(sample) %>%
  summarise(
    n             = n(),
    af_median     = median(af, na.rm = TRUE),
    af_mean       = mean(af, na.rm = TRUE),
    af_sd         = sd(af, na.rm = TRUE),
    af_skewness   = (mean(af, na.rm = TRUE) - median(af, na.rm = TRUE)) / sd(af, na.rm = TRUE),
    af_q10        = quantile(af, 0.10, na.rm = TRUE),
    af_q25        = quantile(af, 0.25, na.rm = TRUE),
    af_q75        = quantile(af, 0.75, na.rm = TRUE),
    af_q90        = quantile(af, 0.90, na.rm = TRUE),
    frac_af_near_50 = sum(af >= 0.40 & af <= 0.60) / n(),
    frac_af_near_33 = sum(af >= 0.25 & af <= 0.40) / n(),
    frac_af_near_67 = sum(af >= 0.60 & af <= 0.75) / n(),
    .groups = "drop"
  ) %>%
  mutate(
    af_pattern = case_when(
      n < 20 ~ "too_few",
      frac_af_near_50 >= 0.50 ~ "peak_at_50",
      frac_af_near_33 >= 0.25 & frac_af_near_67 >= 0.25 ~ "bimodal_33_67",
      af_median < 0.35 | af_median > 0.65 ~ "skewed_peak",
      TRUE ~ "broad"
    )
  )

af_shape <- tibble(sample = all_samples) %>%
  left_join(af_shape_raw, by = "sample") %>%
  mutate(
    n = replace_na(n, 0L),
    af_pattern = replace_na(af_pattern, "no_het")
  ) %>%
  mutate(
    sample_num    = readr::parse_number(sample),
    sample_suffix = stringr::str_extract(sample, "[A-Za-z]+$")
  ) %>%
  arrange(sample_num, sample_suffix) %>%
  select(-sample_num, -sample_suffix)

############################################
### 5.4 Cross-reference with Step 4 depth
### 修正: chrom_het_z 纳入判断
############################################

chr_centered_file <- file.path(depth_dir, "mosdepth_chr_centered_depth_long.tsv")
chr_centered_4 <- read_tsv(chr_centered_file, show_col_types = FALSE)

depth_chrom_ref <- chr_centered_4 %>%
  select(sample, chrom, chr_depth_norm, chrom_centered_ratio,
         chrom_robust_z, chrom_robust_flag) %>%
  distinct()

# ---- 用 depth 作为主表, het 作为辅表 ----
# 这样 depth gain 但 0 het 的 chrom 不会被丢掉
het_depth_cross <- depth_chrom_ref %>%
  left_join(
    het_chrom_summary %>%
      select(sample, chrom, n_het, het_per_kb_chrom,
             chrom_het_z, chrom_het_flag),
    by = c("sample", "chrom")
  ) %>%
  # 没有 het SNP 的 chrom, 补 0
  mutate(
    n_het          = replace_na(n_het, 0L),
    het_per_kb_chrom = replace_na(het_per_kb_chrom, 0),
    chrom_het_flag = replace_na(chrom_het_flag, "no_het")
  ) %>%
  mutate(
    event_class = case_when(
      is.na(chr_depth_norm) ~ NA_character_,
      
      # depth 升高 + het 显著富集 → heterologous disomy
      chrom_robust_flag == "gain_outlier" & chrom_het_flag == "het_enriched" ~
        "disomy_with_het",
      
      # depth 升高 + het 无或 baseline → isogenic (同源) disomy
      chrom_robust_flag == "gain_outlier" & chrom_het_flag %in% c("no_het", "baseline") ~
        "disomy_no_het",
      
      # depth 正常但 het 显著富集 → 混合感染信号
      chrom_robust_flag == "inlier" & chrom_het_flag == "het_enriched" ~
        "mixed_infection_signal",
      
      # depth 降低 + het 少 → loss
      chrom_robust_flag == "loss_outlier" ~
        "loss_candidate",
      
      TRUE ~ "neutral"
    )
  ) %>%
  mutate(
    sample_num    = readr::parse_number(sample),
    sample_suffix = stringr::str_extract(sample, "[A-Za-z]+$"),
    chrom_num     = readr::parse_number(chrom)
  ) %>%
  arrange(sample_num, sample_suffix, chrom_num) %>%
  select(-sample_num, -sample_suffix, -chrom_num)

############################################
### 5.5 Final integrated classification
### 修正: 基于 cohort baseline 的相对判断
############################################

clonality_call <- het_genome_summary %>%
  select(sample, total_het, het_per_kb, het_robust_z, genome_het_flag) %>%
  left_join(
    het_uniformity %>%
      select(sample, het_distribution, max_chrom_het_frac,
             chrom_with_max_het, cv_het_per_kb, frac_chrom_with_het,
             n_chrom_het_enriched),
    by = "sample"
  ) %>%
  left_join(
    af_shape %>%
      select(sample, af_median, af_sd, af_pattern),
    by = "sample"
  ) %>%
  left_join(
    het_depth_cross %>%
      group_by(sample) %>%
      summarise(
        n_disomy_with_het     = sum(event_class == "disomy_with_het", na.rm = TRUE),
        n_disomy_no_het       = sum(event_class == "disomy_no_het", na.rm = TRUE),
        n_mixed_signal_chroms = sum(event_class == "mixed_infection_signal", na.rm = TRUE),
        .groups = "drop"
      ),
    by = "sample"
  ) %>%
  mutate(
    across(
      c(n_disomy_with_het, n_disomy_no_het, n_mixed_signal_chroms,
        n_chrom_het_enriched),
      ~ replace_na(.x, 0L)
    )
  ) %>%
  mutate(
    clonality = case_when(
      # 没有 het 或极少, 也没有 depth 异常 → 纯单克隆
      genome_het_flag %in% c("no_het", "baseline") &
        n_disomy_with_het == 0 & n_disomy_no_het == 0 &
        n_chrom_het_enriched == 0 ~
        "monoclonal_clean",
      
      # baseline het + depth gain (无论 het 是否富集) → 单克隆 disomy
      genome_het_flag %in% c("no_het", "baseline") &
        (n_disomy_with_het >= 1 | n_disomy_no_het >= 1) ~
        "monoclonal_with_disomy",
      
      # het 整体偏高 + 均匀分散 + 多条 chr 富集 → 混合感染
      genome_het_flag %in% c("elevated", "high") &
        het_distribution %in% c("dispersed_uniform", "dispersed_with_hotspot") &
        n_mixed_signal_chroms >= 3 ~
        "mixed_infection_likely",
      
      # het 整体偏高 + 均匀分散但信号不够强
      genome_het_flag == "elevated" &
        het_distribution %in% c("dispersed_uniform", "dispersed_with_hotspot") ~
        "mixed_infection_possible",
      
      # het baseline 但个别 chr 有 het 富集 (无 depth 支持)
      genome_het_flag == "baseline" & n_chrom_het_enriched >= 1 &
        n_disomy_with_het == 0 & n_disomy_no_het == 0 ~
        "localized_het_anomaly",
      
      # het 集中
      het_distribution == "concentrated" ~
        "concentrated_het",
      
      TRUE ~ "ambiguous_needs_review"
    ),
    
    confidence = case_when(
      clonality == "monoclonal_clean" ~ "high",
      clonality == "mixed_infection_likely" & het_robust_z > 5 ~ "high",
      clonality == "monoclonal_with_disomy" ~ "high",
      str_detect(clonality, "ambiguous") ~ "low",
      TRUE ~ "moderate"
    )
  ) %>%
  mutate(
    sample_num    = readr::parse_number(sample),
    sample_suffix = stringr::str_extract(sample, "[A-Za-z]+$")
  ) %>%
  arrange(sample_num, sample_suffix) %>%
  select(-sample_num, -sample_suffix)

# =========================
# Save all tables
# =========================
write_tsv(het_filt,           file.path(out_dir5, "het_snp_filtered_long.tsv"))
write_tsv(het_genome_summary, file.path(out_dir5, "het_genome_summary.tsv"))
write_tsv(het_chrom_summary,  file.path(out_dir5, "het_chrom_summary.tsv"))
write_tsv(het_uniformity,     file.path(out_dir5, "het_uniformity.tsv"))
write_tsv(af_shape,           file.path(out_dir5, "af_shape_summary.tsv"))
write_tsv(het_depth_cross,    file.path(out_dir5, "het_depth_cross_reference.tsv"))
write_tsv(clonality_call,     file.path(out_dir5, "clonality_classification.tsv"))

message("Wrote 7 tables to: ", out_dir5)

############################################
### 5.6 Plots
############################################

# ----- 5.6a Genome-wide AF histogram per sample -----
plot_af_histogram <- function(df_sample, bins = 50) {
  sample_name <- unique(df_sample$sample)
  n_snps <- nrow(df_sample)
  
  ggplot(df_sample, aes(x = af)) +
    geom_histogram(bins = bins, fill = "steelblue", color = "white", linewidth = 0.1) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "red", linewidth = 0.5) +
    geom_vline(xintercept = c(1/3, 2/3), linetype = "dotted", color = "darkorange", linewidth = 0.4) +
    labs(
      title = paste0(sample_name, "  (n = ", n_snps, " het SNPs)"),
      subtitle = "Red=0.5 (1:1 mix / disomy)  Orange=1/3,2/3 (trisomy)",
      x = "Alt allele frequency", y = "Count"
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank())
}

het_split <- het_filt %>% group_split(sample)

for (df_s in het_split) {
  s <- unique(df_s$sample)
  if (nrow(df_s) < 5) next
  
  p <- plot_af_histogram(df_s)
  ggsave(
    filename = file.path(plot_dir5, paste0(s, "_af_histogram.png")),
    plot = p, width = 8, height = 4, dpi = 300
  )
}

# ----- 5.6b Per-chromosome AF histogram (faceted) -----
plot_af_by_chrom <- function(df_sample) {
  sample_name <- unique(df_sample$sample)
  
  chrom_order <- df_sample %>%
    distinct(chrom) %>%
    mutate(chrom_num = readr::parse_number(chrom)) %>%
    arrange(chrom_num, chrom) %>%
    pull(chrom)
  
  df_sample <- df_sample %>%
    mutate(chrom = factor(chrom, levels = chrom_order))
  
  ggplot(df_sample, aes(x = af)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white", linewidth = 0.1) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "red", linewidth = 0.4) +
    geom_vline(xintercept = c(1/3, 2/3), linetype = "dotted", color = "darkorange", linewidth = 0.3) +
    facet_wrap(~ chrom, scales = "free_y") +
    labs(
      title = paste0(sample_name, " — AF by chromosome"),
      subtitle = "Concentrated on 1-2 chr → disomy | Dispersed → mixed infection",
      x = "Alt allele frequency", y = "Count"
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 7)
    )
}

for (df_s in het_split) {
  s <- unique(df_s$sample)
  if (nrow(df_s) < 10) next
  n_chroms <- length(unique(df_s$chrom))
  
  p <- plot_af_by_chrom(df_s)
  ggsave(
    filename = file.path(plot_dir5, paste0(s, "_af_by_chrom.png")),
    plot = p,
    width = max(10, n_chroms * 0.8),
    height = max(6, n_chroms * 0.5),
    dpi = 300
  )
}

# ----- 5.6c Het SNP density barplot per chromosome (faceted by sample) -----
plot_het_density_bar <- function(df_chrom_summary, sample_name) {
  df_plot <- df_chrom_summary %>%
    filter(sample == sample_name) %>%
    mutate(
      chrom_num = readr::parse_number(chrom),
      chrom = fct_reorder(chrom, chrom_num)
    )
  
  ggplot(df_plot, aes(x = chrom, y = het_per_kb_chrom)) +
    geom_col(fill = "steelblue") +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.4) +
    labs(
      title = paste0(sample_name, " — het SNP density per chromosome"),
      subtitle = "Red line = 0.05 het/kb threshold",
      x = "Chromosome", y = "Het SNPs per kb"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 60, hjust = 1)
    )
}

for (s in unique(het_chrom_summary$sample)) {
  if (sum(het_chrom_summary$sample == s) == 0) next
  
  p <- plot_het_density_bar(het_chrom_summary, s)
  ggsave(
    filename = file.path(plot_dir5, paste0(s, "_het_density_bar.png")),
    plot = p, width = 8, height = 4.5, dpi = 300
  )
}

# ----- 5.6d AF genome-wide scatter (reuse chrom_ref from Step 4.2) -----
# chrom_ref 需要在环境中, 如果不在则重建
if (!exists("chrom_ref")) {
  chrom_ref <- chrom_lengths %>%
    mutate(chrom_num = readr::parse_number(chrom)) %>%
    arrange(chrom_num, chrom) %>%
    mutate(
      chrom_offset = lag(cumsum(chr_length), default = 0),
      chrom_center = chrom_offset + chr_length / 2
    ) %>%
    select(-chrom_num)
}

plot_af_genome_scatter <- function(df_sample, chrom_ref_tbl) {
  sample_name <- unique(df_sample$sample)
  
  df_plot <- df_sample %>%
    left_join(chrom_ref_tbl, by = "chrom") %>%
    mutate(genome_x = chrom_offset + pos)
  
  ggplot(df_plot, aes(x = genome_x, y = af)) +
    geom_point(size = 0.3, alpha = 0.4) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", linewidth = 0.4) +
    geom_hline(yintercept = c(1/3, 2/3), linetype = "dotted", color = "darkorange", linewidth = 0.3) +
    geom_vline(
      data = chrom_ref_tbl %>% filter(chrom_offset > 0),
      aes(xintercept = chrom_offset),
      inherit.aes = FALSE, linewidth = 0.25, alpha = 0.5
    ) +
    scale_x_continuous(
      breaks = chrom_ref_tbl$chrom_center,
      labels = chrom_ref_tbl$chrom
    ) +
    labs(
      title = paste0(sample_name, " — AF across genome"),
      subtitle = "Uniform scatter → mixed infection | Clustered on 1 chr → disomy",
      x = "Genome position", y = "Alt allele frequency"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)
    )
}

for (df_s in het_split) {
  s <- unique(df_s$sample)
  if (nrow(df_s) < 10) next
  
  p <- plot_af_genome_scatter(df_s, chrom_ref)
  ggsave(
    filename = file.path(plot_dir5, paste0(s, "_af_genome_scatter.png")),
    plot = p, width = 12, height = 4.8, dpi = 300
  )
}

message("All plots saved to: ", plot_dir5)

# =========================
# Quick inspection 5
# =========================
het_genome_summary %>%
  arrange(desc(het_per_kb)) %>%
  print(n = 30)

het_uniformity %>%
  filter(het_distribution != "too_few") %>%
  arrange(het_distribution, desc(cv_het_per_kb)) %>%
  print(n = 30)

clonality_call %>%
  count(clonality, confidence) %>%
  print(n = 20)

clonality_call %>%
  filter(clonality != "monoclonal_clean") %>%
  arrange(clonality, desc(het_per_kb)) %>%
  print(n = 50)

het_depth_cross %>%
  filter(event_class != "neutral") %>%
  arrange(event_class, sample, chrom) %>%
  print(n = 50)

# cohort baseline 多少
message("Cohort median het: ", cohort_median_het)

# 分类分布
clonality_call %>% count(clonality, confidence)

# 真正偏高的样本
het_genome_summary %>%
  filter(genome_het_flag != "baseline") %>%
  arrange(desc(het_robust_z))

############################################
### Working-set view: 19 株 (20 株工作集 - H99)
### Subset for 19 isolates
############################################

stopifnot(exists("chr_depth_long"), exists("cnv_review_template"), exists("clonality_call"))

# ---- 19 株 + 分组 (H99 standard 本批不含, 排除) ----
working_groups <- tibble::tribble(
  ~id,  ~group,
  101, "HR", 106, "HR", 126, "HR", 128, "HR", 160, "HR",
  162, "HR", 208, "HR", 213, "HR", 215, "HR", 217, "HR",
  147, "control", 150, "control", 169, "control",
  174, "control", 225, "control", 226, "control",
  196, "R", 197, "R", 198, "R"
)
working_samples <- paste0("Cn", working_groups$id)   # 命名约定 CnNNN
group_lookup <- working_groups %>%
  mutate(sample = paste0("Cn", id)) %>% select(sample, group)

# ---- 自检: 19 株是否都在数据里 (命名对不上会 0 匹配, 这里直接报警) ----
all_seen <- unique(chr_depth_long$sample)
missing  <- setdiff(working_samples, all_seen)
message("Working set matched: ", length(intersect(working_samples, all_seen)),
        "/", length(working_samples))
if (length(missing) > 0)
  message("  NOT found (check prefix/naming): ", paste(missing, collapse = ", "))

# ---- 输出目录 ----
depth_dir <- DEPTH_DIR          # from wgs_analysis_config.R
ws_dir    <- file.path(depth_dir, "review_ws19")
dir.create(ws_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1) CNV review template 子集 (带 group) ----
cnv_review_template_ws <- cnv_review_template %>%
  filter(sample %in% working_samples) %>%
  left_join(group_lookup, by = "sample") %>%
  relocate(group, .after = sample)
write_tsv(cnv_review_template_ws, file.path(ws_dir, "cnv_event_review_template_ws19.tsv"))

# ---- 2) Step5 clonality 子集 (带 group) ----
clonality_call_ws <- clonality_call %>%
  filter(sample %in% working_samples) %>%
  left_join(group_lookup, by = "sample") %>%
  relocate(group, .after = sample)
write_tsv(clonality_call_ws, file.path(ws_dir, "clonality_classification_ws19.tsv"))

# ---- 2b) 可选: per-chrom event 证据子集 (给 Fig6 2x2 用) ----
if (exists("het_depth_cross")) {
  het_depth_cross %>%
    filter(sample %in% working_samples) %>%
    left_join(group_lookup, by = "sample") %>%
    relocate(group, .after = sample) %>%
    write_tsv(file.path(ws_dir, "het_depth_cross_ws19.tsv"))
}

# ---- 3) 一屏 triage: 19 株要看什么 ----
message("\n--- CNV review load (19 strains) ---")
print(count(cnv_review_template_ws, pattern_confidence))

message("\n--- per-strain rows still needing manual review ---")
cnv_review_template_ws %>%
  filter(pattern_confidence == "needs_manual_review") %>%
  count(sample, group, sort = TRUE) %>% print(n = Inf)

message("\n--- Step5 clonality (sorted by het_robust_z) ---")
clonality_call_ws %>%
  select(sample, group, total_het, het_robust_z, genome_het_flag, clonality, confidence) %>%
  arrange(desc(het_robust_z)) %>% print(n = Inf)

# ---- 4) 把 19 株已生成的图收集到一个文件夹, 只浏览这些 ----
plot_src_dirs <- path.expand(c(
  file.path(depth_dir, "cnv_genomewide_plots"),
  file.path(depth_dir, "cnv_genomewide_plots_karyo"),
  file.path(dirname(VCF_FILTERED_DIR), "step5_het_analysis", "plots")
))
fig_out <- file.path(ws_dir, "plots")
dir.create(fig_out, showWarnings = FALSE, recursive = TRUE)

copied <- 0L
for (d in plot_src_dirs) {
  if (!dir.exists(d)) next
  pngs <- list.files(d, pattern = "\\.png$", full.names = TRUE)
  for (s in working_samples) {
    hit <- pngs[startsWith(basename(pngs), paste0(s, "_"))]
    if (length(hit)) { file.copy(hit, fig_out, overwrite = TRUE); copied <- copied + length(hit) }
  }
}
message("\nCopied ", copied, " plot PNGs for 19 strains → ", fig_out)

no_plots <- working_samples[!vapply(working_samples, function(s)
  any(startsWith(list.files(fig_out), paste0(s, "_"))), logical(1))]
if (length(no_plots) > 0)
  message("  No plots for: ", paste(no_plots, collapse = ", "),
          " (likely non-candidate / few het SNPs)")

