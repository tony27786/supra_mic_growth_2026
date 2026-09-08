library(tidyverse)
source("./wgs_analysis_config.R")  # shared batch paths (edit BATCH_TAG there)

##########################
### STEP 01 mapping QC ###
##########################

# =========================
# Directory
# =========================
map_dir <- MAP_DIR
metrics_dir <- file.path(map_dir, "qc_metrics_map")
flag_dir <- file.path(metrics_dir, "flagstat")
stats_dir <- file.path(metrics_dir, "stats")
idx_dir <- file.path(metrics_dir, "idxstats")

# =========================
# Functions
# =========================
read_flagstat_one <- function(file) {
  x <- readLines(file, warn = FALSE)
  sample <- sub("\\.flagstat\\.txt$", "", basename(file))
  
  get_n1 <- function(line) {
    suppressWarnings(as.numeric(gsub(",", "", sub(" .*", "", line))))
  }
  
  get_pct <- function(line) {
    out <- sub(".*\\(([^%]+)%.*", "\\1", line)
    suppressWarnings(as.numeric(out))
  }
  
  total <- get_n1(x[1])
  
  primary_line <- x[grepl(" primary$", x)]
  primary <- if (length(primary_line)) get_n1(primary_line[1]) else NA_real_
  
  secondary_line <- x[grepl(" secondary$", x)]
  secondary <- if (length(secondary_line)) get_n1(secondary_line[1]) else NA_real_
  
  supp_line <- x[grepl(" supplementary$", x)]
  supplementary <- if (length(supp_line)) get_n1(supp_line[1]) else NA_real_
  
  dup_line <- x[grepl(" duplicates$", x)]
  duplicates <- if (length(dup_line)) get_n1(dup_line[1]) else NA_real_
  
  pm_line <- x[grepl(" primary mapped \\(", x)]
  m_line  <- x[grepl(" mapped \\(", x) & !grepl("primary mapped", x)]
  
  mapped_primary_reads <- if (length(pm_line)) get_n1(pm_line[1]) else NA_real_
  mapped_primary_pct <- if (length(pm_line)) get_pct(pm_line[1]) else NA_real_
  
  mapped_reads <- if (length(m_line)) get_n1(m_line[1]) else NA_real_
  mapped_pct <- if (length(m_line)) get_pct(m_line[1]) else NA_real_
  
  pp_line <- x[grepl(" properly paired \\(", x)]
  properly_paired_reads <- if (length(pp_line)) get_n1(pp_line[1]) else NA_real_
  properly_paired_pct <- if (length(pp_line)) get_pct(pp_line[1]) else NA_real_
  
  sing_line <- x[grepl(" singletons \\(", x)]
  singleton_reads <- if (length(sing_line)) get_n1(sing_line[1]) else NA_real_
  singleton_pct <- if (length(sing_line)) get_pct(sing_line[1]) else NA_real_
  
  diffchr_line <- x[grepl(" with mate mapped to a different chr$", x)]
  diffchr <- if (length(diffchr_line)) get_n1(diffchr_line[1]) else NA_real_
  
  diffchr_q5_line <- x[grepl(" with mate mapped to a different chr \\(mapQ>=5\\)$", x)]
  diffchr_mapq5 <- if (length(diffchr_q5_line)) get_n1(diffchr_q5_line[1]) else NA_real_
  
  tibble(
    sample = sample,
    total_reads = total,
    primary_reads = primary,
    secondary_reads = secondary,
    supplementary_reads = supplementary,
    duplicate_reads = duplicates,
    mapped_reads = mapped_reads,
    mapped_pct = mapped_pct,
    mapped_primary_reads = mapped_primary_reads,
    mapped_primary_pct = mapped_primary_pct,
    properly_paired_reads = properly_paired_reads,
    properly_paired_pct = properly_paired_pct,
    singleton_reads = singleton_reads,
    singleton_pct = singleton_pct,
    mate_diff_chr_reads = diffchr,
    mate_diff_chr_mapq5_reads = diffchr_mapq5
  )
}

read_stats_one <- function(file) {
  sample <- sub("\\.stats\\.txt$", "", basename(file))
  
  x <- read_tsv(
    file,
    comment = "#",
    col_names = FALSE,
    show_col_types = FALSE,
    progress = FALSE
  )
  
  sn <- x %>%
    filter(X1 == "SN") %>%
    transmute(
      key = stringr::str_trim(stringr::str_remove(X2, ":$")),
      value = stringr::str_trim(X3)
    )
  
  get_num <- function(key_name) {
    y <- sn %>% filter(key == key_name) %>% pull(value)
    if (length(y) == 0) return(NA_real_)
    suppressWarnings(as.numeric(y[1]))
  }
  
  tibble(
    sample = sample,
    raw_total_sequences = get_num("raw total sequences"),
    filtered_sequences = get_num("filtered sequences"),
    sequences_stats = get_num("sequences"),
    reads_mapped_stats = get_num("reads mapped"),
    reads_unmapped_stats = get_num("reads unmapped"),
    reads_properly_paired_stats = get_num("reads properly paired"),
    reads_paired_stats = get_num("reads paired"),
    reads_duplicated_stats = get_num("reads duplicated"),
    error_rate = get_num("error rate"),
    average_length = get_num("average length"),
    average_quality = get_num("average quality"),
    insert_size_average = get_num("insert size average"),
    insert_size_sd = get_num("insert size standard deviation")
  )
}

read_idxstats_one <- function(file) {
  sample <- sub("\\.idxstats\\.txt$", "", basename(file))
  
  x <- read_tsv(
    file,
    col_names = c("chrom", "chrom_length", "mapped", "unmapped"),
    show_col_types = FALSE,
    progress = FALSE
  )
  
  x_chr <- x %>% filter(chrom != "*")
  x_star <- x %>% filter(chrom == "*")
  
  tibble(
    sample = sample,
    n_reference_seqs = nrow(x_chr),
    ref_total_length = sum(x_chr$chrom_length, na.rm = TRUE),
    mapped_chr_reads = sum(x_chr$mapped, na.rm = TRUE),
    unmapped_reads_idxstats = ifelse(nrow(x_star) > 0, x_star$unmapped[1], NA_real_),
    max_chr_mapped_reads = max(x_chr$mapped, na.rm = TRUE),
    min_chr_mapped_reads = min(x_chr$mapped, na.rm = TRUE)
  )
}

# =========================
# Read all mapping QC
# =========================
flag_files <- list.files(flag_dir, pattern = "\\.flagstat\\.txt$", full.names = TRUE)
stats_files <- list.files(stats_dir, pattern = "\\.stats\\.txt$", full.names = TRUE)
idx_files <- list.files(idx_dir, pattern = "\\.idxstats\\.txt$", full.names = TRUE)

flag_df <- purrr::map_dfr(flag_files, read_flagstat_one)
stats_df <- purrr::map_dfr(stats_files, read_stats_one)
idx_df <- purrr::map_dfr(idx_files, read_idxstats_one)

mapping_qc <- flag_df %>%
  full_join(stats_df, by = "sample") %>%
  full_join(idx_df, by = "sample") %>%
  mutate(
    mapped_pct_use = if_else(!is.na(mapped_primary_pct), mapped_primary_pct, mapped_pct),
    singleton_pct = if_else(
      is.na(singleton_pct) & !is.na(singleton_reads) & !is.na(total_reads) & total_reads > 0,
      100 * singleton_reads / total_reads,
      singleton_pct
    ),
    mate_diff_chr_mapq5_pct = if_else(
      !is.na(mate_diff_chr_mapq5_reads) & !is.na(total_reads) & total_reads > 0,
      100 * mate_diff_chr_mapq5_reads / total_reads,
      NA_real_
    ),
    duplicate_pct = if_else(
      !is.na(duplicate_reads) & !is.na(total_reads) & total_reads > 0,
      100 * duplicate_reads / total_reads,
      NA_real_
    ),
    mapped_chr_pct_from_idx = if_else(
      !is.na(mapped_chr_reads) & !is.na(unmapped_reads_idxstats) &
        (mapped_chr_reads + unmapped_reads_idxstats) > 0,
      100 * mapped_chr_reads / (mapped_chr_reads + unmapped_reads_idxstats),
      NA_real_
    ),
    flag_low_mapped = !is.na(mapped_pct_use) & mapped_pct_use < 95,
    flag_low_pp = !is.na(properly_paired_pct) & properly_paired_pct < 90,
    flag_high_singleton = !is.na(singleton_pct) & singleton_pct > 5,
    flag_high_diffchr = !is.na(mate_diff_chr_mapq5_pct) & mate_diff_chr_mapq5_pct > 1,
    flag_insert_size = !is.na(insert_size_average) & (insert_size_average < 100 | insert_size_average > 1000),
    mapping_qc_class = case_when(
      flag_low_mapped | flag_low_pp | flag_insert_size ~ "HOLD",
      flag_high_singleton | flag_high_diffchr ~ "WARN",
      TRUE ~ "PASS"
    )
  ) %>%
  arrange(readr::parse_number(sample), sample)

write_tsv(mapping_qc, file.path(map_dir, "mapping_qc_summary.tsv"))
mapping_qc %>%
  filter(mapping_qc_class != "PASS") %>%
  write_tsv(file.path(map_dir, "mapping_qc_bad.tsv"))

###########################
### STEP 02 mosdepth QC ###
###########################

# =========================
# Directory
# =========================
mos_dir <- MOS_DIR
mos_out_file <- file.path(dirname(mos_dir), "mosdepth_master_qc.tsv")

# =========================
# 1) Read one mosdepth summary
# =========================
read_mos_summary_one <- function(file) {
  sample <- sub("\\.mosdepth\\.summary\\.txt$", "", basename(file))
  
  x <- read_tsv(
    file,
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    filter(!str_detect(chrom, "_region$")) %>%
    mutate(
      length = as.numeric(length),
      bases = as.numeric(bases),
      mean = as.numeric(mean),
      min = as.numeric(min),
      max = as.numeric(max)
    ) %>%
    filter(!is.na(length), !is.na(bases), !is.na(mean)) %>%
    mutate(
      is_nuclear = length >= 1e5
    )
  
  if (nrow(x) == 0) {
    return(tibble(
      sample = sample,
      n_contigs_all = NA_integer_,
      genome_total_length_all = NA_real_,
      genome_total_bases_all = NA_real_,
      genome_mean_depth_all = NA_real_,
      chr_mean_depth_mean_all = NA_real_,
      chr_mean_depth_median_all = NA_real_,
      chr_mean_depth_min_all = NA_real_,
      chr_mean_depth_max_all = NA_real_,
      chr_mean_depth_sd_all = NA_real_,
      chr_mean_depth_cv_all = NA_real_,
      n_contigs_nuclear = NA_integer_,
      genome_total_length_nuclear = NA_real_,
      genome_total_bases_nuclear = NA_real_,
      genome_mean_depth_nuclear = NA_real_,
      chr_mean_depth_mean_nuclear = NA_real_,
      chr_mean_depth_median_nuclear = NA_real_,
      chr_mean_depth_min_nuclear = NA_real_,
      chr_mean_depth_max_nuclear = NA_real_,
      chr_mean_depth_sd_nuclear = NA_real_,
      chr_mean_depth_cv_nuclear = NA_real_
    ))
  }
  
  chr_mean_depth_mean_all <- mean(x$mean, na.rm = TRUE)
  chr_mean_depth_sd_all <- sd(x$mean, na.rm = TRUE)
  
  x_nuc <- x %>% filter(is_nuclear)
  
  if (nrow(x_nuc) == 0) {
    n_contigs_nuclear <- NA_integer_
    genome_total_length_nuclear <- NA_real_
    genome_total_bases_nuclear <- NA_real_
    genome_mean_depth_nuclear <- NA_real_
    chr_mean_depth_mean_nuclear <- NA_real_
    chr_mean_depth_median_nuclear <- NA_real_
    chr_mean_depth_min_nuclear <- NA_real_
    chr_mean_depth_max_nuclear <- NA_real_
    chr_mean_depth_sd_nuclear <- NA_real_
    chr_mean_depth_cv_nuclear <- NA_real_
  } else {
    chr_mean_depth_mean_nuclear <- mean(x_nuc$mean, na.rm = TRUE)
    chr_mean_depth_sd_nuclear <- sd(x_nuc$mean, na.rm = TRUE)
    
    n_contigs_nuclear <- nrow(x_nuc)
    genome_total_length_nuclear <- sum(x_nuc$length, na.rm = TRUE)
    genome_total_bases_nuclear <- sum(x_nuc$bases, na.rm = TRUE)
    genome_mean_depth_nuclear <- genome_total_bases_nuclear / genome_total_length_nuclear
    chr_mean_depth_median_nuclear <- median(x_nuc$mean, na.rm = TRUE)
    chr_mean_depth_min_nuclear <- min(x_nuc$mean, na.rm = TRUE)
    chr_mean_depth_max_nuclear <- max(x_nuc$mean, na.rm = TRUE)
    chr_mean_depth_cv_nuclear <- ifelse(
      !is.na(chr_mean_depth_mean_nuclear) && chr_mean_depth_mean_nuclear > 0,
      chr_mean_depth_sd_nuclear / chr_mean_depth_mean_nuclear,
      NA_real_
    )
  }
  
  tibble(
    sample = sample,
    n_contigs_all = nrow(x),
    genome_total_length_all = sum(x$length, na.rm = TRUE),
    genome_total_bases_all = sum(x$bases, na.rm = TRUE),
    genome_mean_depth_all = sum(x$bases, na.rm = TRUE) / sum(x$length, na.rm = TRUE),
    chr_mean_depth_mean_all = chr_mean_depth_mean_all,
    chr_mean_depth_median_all = median(x$mean, na.rm = TRUE),
    chr_mean_depth_min_all = min(x$mean, na.rm = TRUE),
    chr_mean_depth_max_all = max(x$mean, na.rm = TRUE),
    chr_mean_depth_sd_all = chr_mean_depth_sd_all,
    chr_mean_depth_cv_all = ifelse(
      !is.na(chr_mean_depth_mean_all) && chr_mean_depth_mean_all > 0,
      chr_mean_depth_sd_all / chr_mean_depth_mean_all,
      NA_real_
    ),
    n_contigs_nuclear = n_contigs_nuclear,
    genome_total_length_nuclear = genome_total_length_nuclear,
    genome_total_bases_nuclear = genome_total_bases_nuclear,
    genome_mean_depth_nuclear = genome_mean_depth_nuclear,
    chr_mean_depth_mean_nuclear = chr_mean_depth_mean_nuclear,
    chr_mean_depth_median_nuclear = chr_mean_depth_median_nuclear,
    chr_mean_depth_min_nuclear = chr_mean_depth_min_nuclear,
    chr_mean_depth_max_nuclear = chr_mean_depth_max_nuclear,
    chr_mean_depth_sd_nuclear = chr_mean_depth_sd_nuclear,
    chr_mean_depth_cv_nuclear = chr_mean_depth_cv_nuclear
  )
}

# =========================
# 2) Read one thresholds.bed.gz
# =========================
read_thresholds_one <- function(file) {
  sample <- sub("\\.thresholds\\.bed\\.gz$", "", basename(file))
  
  x <- read_tsv(
    file,
    comment = "#",
    col_names = c("chrom", "start", "end", "region", "X1", "X5", "X10", "X20", "X50", "X100"),
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    mutate(
      start = as.numeric(start),
      end = as.numeric(end),
      width = end - start,
      X1 = as.numeric(X1),
      X5 = as.numeric(X5),
      X10 = as.numeric(X10),
      X20 = as.numeric(X20),
      X50 = as.numeric(X50),
      X100 = as.numeric(X100)
    ) %>%
    filter(!is.na(width), width > 0)
  
  if (nrow(x) == 0) {
    return(tibble(
      sample = sample,
      threshold_total_bases = NA_real_,
      bases_ge_1x = NA_real_,
      bases_ge_5x = NA_real_,
      bases_ge_10x = NA_real_,
      bases_ge_20x = NA_real_,
      bases_ge_50x = NA_real_,
      bases_ge_100x = NA_real_,
      pct_bases_ge_1x = NA_real_,
      pct_bases_ge_5x = NA_real_,
      pct_bases_ge_10x = NA_real_,
      pct_bases_ge_20x = NA_real_,
      pct_bases_ge_50x = NA_real_,
      pct_bases_ge_100x = NA_real_
    ))
  }
  
  total_bases <- sum(x$width, na.rm = TRUE)
  
  tibble(
    sample = sample,
    threshold_total_bases = total_bases,
    bases_ge_1x = sum(x$X1, na.rm = TRUE),
    bases_ge_5x = sum(x$X5, na.rm = TRUE),
    bases_ge_10x = sum(x$X10, na.rm = TRUE),
    bases_ge_20x = sum(x$X20, na.rm = TRUE),
    bases_ge_50x = sum(x$X50, na.rm = TRUE),
    bases_ge_100x = sum(x$X100, na.rm = TRUE),
    pct_bases_ge_1x = 100 * sum(x$X1, na.rm = TRUE) / total_bases,
    pct_bases_ge_5x = 100 * sum(x$X5, na.rm = TRUE) / total_bases,
    pct_bases_ge_10x = 100 * sum(x$X10, na.rm = TRUE) / total_bases,
    pct_bases_ge_20x = 100 * sum(x$X20, na.rm = TRUE) / total_bases,
    pct_bases_ge_50x = 100 * sum(x$X50, na.rm = TRUE) / total_bases,
    pct_bases_ge_100x = 100 * sum(x$X100, na.rm = TRUE) / total_bases
  )
}

# =========================
# Read all mosdepth QC
# =========================
summary_files <- list.files(mos_dir, pattern = "\\.mosdepth\\.summary\\.txt$", full.names = TRUE)
threshold_files <- list.files(mos_dir, pattern = "\\.thresholds\\.bed\\.gz$", full.names = TRUE)

summary_df <- purrr::map_dfr(summary_files, read_mos_summary_one)
threshold_df <- purrr::map_dfr(threshold_files, read_thresholds_one)

mosdepth_qc <- summary_df %>%
  full_join(threshold_df, by = "sample") %>%
  mutate(
    flag_low_depth = !is.na(genome_mean_depth_nuclear) & genome_mean_depth_nuclear < 20,
    flag_depth_warn = !is.na(genome_mean_depth_nuclear) &
      genome_mean_depth_nuclear >= 20 & genome_mean_depth_nuclear < 30,
    flag_low_10x = !is.na(pct_bases_ge_10x) & pct_bases_ge_10x < 95,
    flag_low_20x = !is.na(pct_bases_ge_20x) & pct_bases_ge_20x < 90,
    mos_qc_class = case_when(
      flag_low_depth ~ "HOLD",
      flag_depth_warn | flag_low_10x | flag_low_20x ~ "WARN",
      TRUE ~ "PASS"
    ),
    mos_qc_note = pmap_chr(
      list(flag_low_depth, flag_depth_warn, flag_low_10x, flag_low_20x),
      function(a, b, c, d) {
        notes <- c(
          if (isTRUE(a)) "low_depth_lt20x" else NULL,
          if (isTRUE(b)) "depth_20_to_30x" else NULL,
          if (isTRUE(c)) "low_pct_ge_10x" else NULL,
          if (isTRUE(d)) "low_pct_ge_20x" else NULL
        )
        if (length(notes) == 0) NA_character_ else paste(notes, collapse = "; ")
      }
    )
  ) %>%
  arrange(readr::parse_number(sample), sample)

write_tsv(mosdepth_qc, mos_out_file)

##############################
### STEP 03 Master sample  ###
##############################
master_out <- COMBINED_DIR
master_sample_sheet <- mapping_qc %>%
  full_join(mosdepth_qc, by = "sample") %>%
  mutate(
    master_qc_class = case_when(
      mapping_qc_class == "HOLD" | mos_qc_class == "HOLD" ~ "HOLD",
      mapping_qc_class == "WARN" | mos_qc_class == "WARN" ~ "WARN",
      TRUE ~ "PASS"
    ),
    master_qc_note = pmap_chr(
      list(mapping_qc_class, mos_qc_class, flag_low_mapped, flag_low_pp, flag_insert_size,
           flag_high_singleton, flag_high_diffchr, flag_low_depth, flag_depth_warn,
           flag_low_10x, flag_low_20x),
      function(mc, qc2, a, b, c, d, e, f, g, h, i) {
        notes <- c(
          if (isTRUE(a)) "low_mapped" else NULL,
          if (isTRUE(b)) "low_properly_paired" else NULL,
          if (isTRUE(c)) "abnormal_insert_size" else NULL,
          if (isTRUE(d)) "high_singleton" else NULL,
          if (isTRUE(e)) "high_diffchr" else NULL,
          if (isTRUE(f)) "low_depth_lt20x" else NULL,
          if (isTRUE(g)) "depth_20_to_30x" else NULL,
          if (isTRUE(h)) "low_pct_ge_10x" else NULL,
          if (isTRUE(i)) "low_pct_ge_20x" else NULL
        )
        if (length(notes) == 0) NA_character_ else paste(unique(notes), collapse = "; ")
      }
    )
  ) %>%
  arrange(readr::parse_number(sample), sample)

write_tsv(master_sample_sheet, file.path(master_out, "hetero_master_sample_sheet.tsv"))

########################
### Quick inspection ###
########################

mapping_qc %>% count(mapping_qc_class) %>% print()
mosdepth_qc %>% count(mos_qc_class) %>% print()
master_sample_sheet %>% count(master_qc_class) %>% print()

master_sample_sheet %>%
  select(
    sample,
    mapping_qc_class,
    mos_qc_class,
    master_qc_class,
    mapped_pct_use,
    properly_paired_pct,
    insert_size_average,
    genome_mean_depth_nuclear,
    pct_bases_ge_10x,
    pct_bases_ge_20x,
    chr_mean_depth_cv_nuclear,
    master_qc_note
  ) %>%
  print(n = 30)