# Figure 3 plot code
#
# Inputs: ../SourceData. Default PDF output: ../Figure3_panels.
# Set FIG3_EXPORT=false to calculate without writing PDFs.
# FIG3_OUTPUT_DIR can redirect PDF output; no input files are modified.

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(DESeq2)
  library(ComplexHeatmap)
  library(circlize)
  # DESeq2 pulls in AnnotationDbi/S4Vectors/IRanges, which mask several dplyr and
  # tidyr verbs (count, select, rename, slice, filter, first). Re-attach the tidy
  # packages so those verbs resolve to the intended implementations.
  for (p in c("purrr", "tibble", "readr", "tidyr", "dplyr")) {
    try(detach(paste0("package:", p), character.only = TRUE, unload = FALSE), silent = TRUE)
  }
  library(purrr); library(tibble); library(readr); library(tidyr); library(dplyr)
})

get_script_dir <- function() {
  sourced <- vapply(sys.frames(), function(fr) {
    if (is.null(fr$ofile)) "" else as.character(fr$ofile)[1]
  }, character(1))
  sourced <- sourced[nzchar(sourced)]
  if (length(sourced)) {
    return(dirname(normalizePath(tail(sourced, 1), mustWork = TRUE)))
  }
  args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(args)) {
    return(dirname(normalizePath(sub("^--file=", "", args[1]), mustWork = TRUE)))
  }
  candidates <- file.path(getwd(), c(".", "Codes", "Revise02/Codes"))
  found <- candidates[file.exists(file.path(candidates, "fig3_submit.R"))]
  if (length(found) != 1L) {
    stop("Run fig3_submit.R with source() or Rscript so its directory can be resolved.")
  }
  normalizePath(found, mustWork = TRUE)
}
SCRIPT_DIR <- get_script_dir()
SOURCE_DIR <- normalizePath(file.path(SCRIPT_DIR, "..", "SourceData"), mustWork = TRUE)
OUTPUT_DIR <- Sys.getenv("FIG3_OUTPUT_DIR", file.path(SCRIPT_DIR, "..", "Figure3_panels"))
DERIVED_DIR <- Sys.getenv("FIG3_DERIVED_DIR", file.path(SCRIPT_DIR, "..", "Derived"))
EXPORT_PANELS <- tolower(Sys.getenv("FIG3_EXPORT", "true")) %in% c("1", "true", "yes")

kb_xlsx <- file.path(SOURCE_DIR, "kbpap_SourceData.xlsx")
smg_xlsx <- file.path(SOURCE_DIR, "smg_SourceData.xlsx")
pap_parent_csv <- file.path(SOURCE_DIR, "parent_pap_feature_results.csv")
pap_clone_csv <- file.path(SOURCE_DIR, "selected_pap_feature_results.csv")
master_tsv <- file.path(SOURCE_DIR, "master_wgs_table.tsv")
count_tsv <- file.path(SOURCE_DIR, "gene_counts.txt")
gtf_path <- file.path(SOURCE_DIR, "refs", "GCF_000149245.1_CNA3_genomic.gtf")
gaf_files <- file.path(SOURCE_DIR, "refs", c("FungiDB-68_CneoformansH99_GO.gaf.gz",
                                             "FungiDB-68_CneoformansH99_Curated_GO.gaf.gz"))
input_files <- c(kb_xlsx, smg_xlsx, pap_parent_csv, pap_clone_csv, master_tsv,
                 count_tsv, gtf_path, gaf_files)
if (any(!file.exists(input_files))) {
  stop("Missing SourceData input(s): ",
       paste(basename(input_files[!file.exists(input_files)]), collapse = ", "))
}

require_columns <- function(d, columns, label) {
  missing <- setdiff(columns, names(d))
  if (length(missing)) stop(label, " is missing: ", paste(missing, collapse = ", "))
}
require_unique <- function(d, columns, label) {
  if (anyDuplicated(d[, columns, drop = FALSE])) {
    stop(label, " contains duplicate keys: ", paste(columns, collapse = ", "))
  }
}

SMG_TP <- "72h"   # primary SMG timepoint
MIC_TP <- "48h"   # timepoint for plate-MIC50 read-out
CONC   <- c(128, 64, 32, 16, 8, 4, 2, 1, 0.5)

POS_MIN <- 0.15   # blank-subtracted no-drug control floor; below = replicate not evaluable

safe_med <- function(x) { x <- x[is.finite(x)]; if (length(x)) stats::median(x) else NA_real_ }

mic50_one <- function(conc, od, pos) {
  if (!is.finite(pos) || pos < POS_MIN) return(NA_real_)
  keep <- is.finite(conc) & is.finite(od)
  conc <- conc[keep]; od <- od[keep]
  if (!length(conc)) return(NA_real_)
  o <- order(conc); conc <- conc[o]; od <- od[o]
  h <- which(od <= 0.5 * pos)
  if (length(h)) conc[min(h)] else max(CONC, na.rm = TRUE) * 2
}

tp_by_side <- function(tp, side) {
  if (length(tp) == 1L) rep(unname(tp), length(side)) else unname(tp[as.character(side)])
}


G_HR   <- "HR-candidate"
G_CTRL <- "Control"
G_R    <- "Resistant"
G_H99  <- "H99"
group_levels <- c(G_HR, G_CTRL, G_R, G_H99)
group_colors <- setNames(c("#E69F00", "#009E73", "#D55E00", "#000000"), group_levels)  # Okabe-Ito
side_levels  <- c("parent", "clone")
side_labels <- c(parent = "Parent", clone = "Selected")
theme_f4 <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"), legend.position = "top")
canon <- function(s){ s <- as.character(s); ifelse(nchar(s) > 3, substr(s, 1, nchar(s) - 2), s) }
trapz <- function(x, y){ o <- order(x); x <- x[o]; y <- y[o]; sum(diff(x) * (head(y,-1)+tail(y,-1))/2) }


## 1. KB measurements and isolate metadata
ci <- read_excel(kb_xlsx, sheet = "Clinical_Isolates") %>%
  mutate(strain = as.character(ID))
require_columns(ci, c("strain", "Group", "ori_delta",
                      paste0("re_", 1:5, "_delta")), "Clinical_Isolates")
meta <- ci %>%
  transmute(
    strain,
    group = case_when(
      str_to_lower(as.character(Group)) == "hr"                  ~ G_HR,
      str_to_lower(as.character(Group)) == "control"             ~ G_CTRL,
      str_to_lower(as.character(Group)) == "r"                   ~ G_R,
      str_to_lower(as.character(Group)) %in% c("standard", "h99") ~ G_H99,
      TRUE ~ as.character(Group)
    )
  ) %>%
  mutate(group = factor(group, levels = group_levels))

kb_parent <- ci %>%
  select(strain, ori_delta, re_1_delta, re_2_delta, re_3_delta, re_4_delta, re_5_delta) %>%
  pivot_longer(-strain, names_to = "rep", values_to = "delta") %>% mutate(side = "parent")
x1 <- read_excel(kb_xlsx, sheet = "xxx01") %>% mutate(strain = as.character(ID))
require_columns(x1, c("strain", paste0("re_", 1:3, "_delta")), "xxx01")
kb_clone <- x1 %>%
  select(strain, re_1_delta, re_2_delta, re_3_delta) %>%
  pivot_longer(-strain, names_to = "rep", values_to = "delta") %>% mutate(side = "clone")
kb_long <- bind_rows(kb_parent, kb_clone) %>% filter(!is.na(delta))
require_unique(meta, "strain", "Clinical_Isolates")
require_unique(x1, "strain", "xxx01")
if (!setequal(meta$strain, x1$strain)) {
  stop("Parent and Selected KB isolate IDs do not match.")
}
if (anyNA(meta$group)) stop("Unrecognized group in Clinical_Isolates.")


## 2. PAP from image-analysis project measurements
pap_well_map <- tribble(
  ~plate, ~well, ~flc,
  "L", "TL", 0, "L", "TM", 0.5, "L", "TR", 1,
  "L", "BL", 2, "L", "BM", 4, "L", "BR", 8,
  "H", "TL", 0, "H", "TM", 8, "H", "TR", 16,
  "H", "BL", 32, "H", "BM", 64, "H", "BR", 128,
  "S", "TL", 0, "S", "TM", 8, "S", "TR", 12,
  "S", "BL", 16, "S", "BM", 24, "S", "BR", 32
)

read_pap_features <- function(path, is_clone) {
  d <- read_csv(path, show_col_types = FALSE)
  require_columns(d, c("image", "well", "raw_mean"), basename(path))
  d <- d %>%
    mutate(
      image_base = str_remove(str_remove(as.character(image), "\\.[A-Za-z0-9]+$"), "_crop$"),
      well = as.character(well)
    ) %>%
    filter(!str_detect(image_base, "^000blank")) %>%
    mutate(
      plate = if (is_clone) "S" else str_extract(image_base, "[HL]$"),
      sample_id = if (is_clone) image_base else str_remove(image_base, "[HL]$")
    ) %>%
    left_join(pap_well_map, by = c("plate", "well"))
  if (anyNA(d$flc) || anyNA(d$sample_id) || any(!is.finite(d$raw_mean))) {
    stop("Invalid sample/well identifier or missing raw_mean in ", basename(path))
  }
  d %>%
    group_by(sample_id, flc) %>%
    summarise(raw_mean = mean(raw_mean, na.rm = TRUE), .groups = "drop") %>%
    arrange(sample_id, flc)
}
pap_parent_avg <- read_pap_features(pap_parent_csv, FALSE)
pap_clone_avg <- read_pap_features(pap_clone_csv, TRUE)

# Relative growth (raw/FLC0) makes the two plate normalizations comparable;
# AUC is taken on the common grid {0, 8, 16, 32}.
pap_one <- function(d, is_clone){
  d <- d %>% mutate(strain = if (is_clone) canon(sample_id) else as.character(sample_id))
  long <- d %>% group_by(strain) %>%
    mutate(base = raw_mean[flc == 0][1]) %>% filter(!is.na(base), base != 0) %>%
    mutate(rel = raw_mean / base) %>% select(strain, flc, rel) %>% ungroup()
  auc <- long %>% group_by(strain) %>%
    summarise(pap_auc_native = trapz(flc[flc <= 32], rel[flc <= 32]),
              pap_auc_common = { ok <- all(c(0,8,16,32) %in% flc)
              if (ok) trapz(c(0,8,16,32), rel[match(c(0,8,16,32), flc)]) else NA_real_ },
              .groups = "drop")
  list(long = long, auc = auc)
}
pp <- pap_one(pap_parent_avg, FALSE); pc <- pap_one(pap_clone_avg, TRUE)
pap_long <- bind_rows(mutate(pp$long, side = "parent"), mutate(pc$long, side = "clone"))
pap_auc  <- bind_rows(mutate(pp$auc,  side = "parent"), mutate(pc$auc,  side = "clone"))
require_unique(pap_auc, c("strain", "side"), "PAP AUC")
if (!setequal(pp$auc$strain, meta$strain) || !setequal(pc$auc$strain, meta$strain)) {
  stop("PAP measurements do not cover the KB isolate set on both sides.")
}
if (any(!is.finite(pap_auc$pap_auc_common))) {
  stop("PAP AUC requires a nonzero FLC0 reference and the full 0/8/16/32 grid.")
}


## 3. SMG from raw OD readings
# Each time block: separator, nine FLC wells, Positive, Negative, Negative.
# Subtract the mean of the two negative controls. Preserve negative corrected OD.
read_block <- function(raw, tp, i){
  lab <- 3 + (i - 1) * 13
  conc_col <- (lab + 1):(lab + 9); pos_col <- lab + 10; neg_col <- c(lab + 11, lab + 12)
  blank <- rowMeans(raw[, neg_col], na.rm = TRUE)
  tibble(ID = rep(as.character(raw[[1]]), length(CONC)),
         rep = rep(raw[[2]], length(CONC)),
         timepoint = tp,
         conc = rep(CONC, each = nrow(raw)),
         od  = as.vector(as.matrix(raw[, conc_col])) - rep(blank, length(CONC)),
         pos = rep(raw[[pos_col]] - blank, length(CONC)))
}
read_smg_sheet <- function(sheet, tps, side, is_clone = FALSE){
  raw <- as.data.frame(read_excel(smg_xlsx, sheet, .name_repair = "unique_quiet"))
  d <- map_dfr(seq_along(tps), ~read_block(raw, tps[.x], .x))
  d$side <- side
  if (is_clone) d$ID <- canon(d$ID)
  d
}
smg_all <- bind_rows(
  read_smg_sheet("ori_1", c("24h","48h","72h"),        "parent"),
  read_smg_sheet("ori_2", c("24h","48h","72h","96h"),  "parent"),
  read_smg_sheet("xxx01", c("24h","48h","72h","96h"),  "clone", is_clone = TRUE)) %>%
  mutate(mic_tp = tp_by_side(MIC_TP, side), smg_tp = tp_by_side(SMG_TP, side))

mic_rep_tbl <- smg_all %>%
  filter(timepoint == mic_tp | timepoint == smg_tp) %>%
  group_by(ID, side, rep, mic_tp, smg_tp) %>%
  summarise(
    mic_primary  = { tp <- first(mic_tp)
                     mic50_one(conc[timepoint == tp], od[timepoint == tp], first(pos[timepoint == tp])) },
    mic_fallback = { tp <- first(smg_tp)
                     mic50_one(conc[timepoint == tp], od[timepoint == tp], first(pos[timepoint == tp])) },
    .groups = "drop") %>%
  mutate(MIC_rep     = dplyr::coalesce(mic_primary, mic_fallback),
         mic_tp_used = dplyr::if_else(is.na(mic_primary), smg_tp, mic_tp))

mic_plate <- mic_rep_tbl %>%
  group_by(ID, side) %>%
  summarise(MIC_plate      = safe_med(MIC_rep),
            n_mic_rep      = sum(is.finite(MIC_rep)),
            n_mic_fallback = sum(mic_tp_used == smg_tp & is.finite(MIC_rep)),
            .groups = "drop")

mic_use <- mic_plate %>% transmute(ID, side, MIC_use = MIC_plate)

mic_dropped <- mic_plate %>% filter(!is.finite(MIC_plate))
if (nrow(mic_dropped)) message("[v2] no plate MIC (dropped from SMG): ",
  paste(mic_dropped$ID, mic_dropped$side, sep = "/", collapse = ", "))
message("[v2] MIC replicates using the 72h fallback: ", sum(mic_plate$n_mic_fallback))

smg_rep <- smg_all %>% filter(timepoint == smg_tp) %>%
  left_join(mic_use, by = c("ID", "side")) %>%
  group_by(ID, side, rep, MIC_use) %>%
  summarise(SMG = { m <- first(MIC_use); p <- first(pos)
                    if (!is.finite(m) || !is.finite(p) || p < POS_MIN) NA_real_
                    else { ab <- od[conc > m & is.finite(od)]
                           if (length(ab)) mean(ab) / p else NA_real_ } },
            .groups = "drop") %>%
  rename(strain = ID) %>% mutate(side = factor(side, levels = side_levels))
if (!setequal(unique(smg_all$ID), meta$strain)) stop("SMG and KB isolate IDs differ.")


## 4. Per-isolate summaries and ordering
# Panel A uses medians; the displacement in Panel D uses the mean KB readings.
# SMG is the mean of evaluable repeats, separately for each side.
ms <- function(df, val, s) {
  df %>% filter(side == s) %>% group_by(strain) %>%
    summarise(v = mean(.data[[val]], na.rm = TRUE), .groups = "drop") %>% deframe()
}
master_out <- meta %>%
  mutate(
    kb_delta_parent = ms(kb_long, "delta", "parent")[strain],
    kb_delta_clone = ms(kb_long, "delta", "clone")[strain],
    pap_auc_parent = ms(pap_auc, "pap_auc_common", "parent")[strain],
    pap_auc_clone = ms(pap_auc, "pap_auc_common", "clone")[strain],
    smg_parent = ms(smg_rep, "SMG", "parent")[strain],
    smg_clone = ms(smg_rep, "SMG", "clone")[strain],
    d_kb = kb_delta_clone - kb_delta_parent,
    d_pap = pap_auc_clone - pap_auc_parent,
    d_smg = smg_clone - smg_parent
  ) %>% arrange(group, strain)

strain_order <- master_out %>%
  mutate(
    .zk = d_kb / sd(d_kb, na.rm = TRUE),
    .zp = d_pap / sd(d_pap, na.rm = TRUE),
    .shift = sqrt(coalesce(.zk, 0)^2 + coalesce(.zp, 0)^2)
  ) %>% arrange(group, desc(.shift)) %>% pull(strain)
add_grp <- function(df) df %>% left_join(meta %>% select(strain, group), by = "strain") %>%
  mutate(group = factor(group, levels = group_levels),
         strain = factor(strain, levels = strain_order),
         side = factor(side, levels = side_levels))


## 5. Panel A: KB delta by group
kb_pm  <- add_grp(kb_long) %>%
  group_by(group, strain, side) %>% summarise(delta = median(delta, na.rm = TRUE), .groups = "drop")
kb_pair <- kb_pm %>% group_by(strain) %>% filter(all(side_levels %in% side)) %>% ungroup() %>%
  filter(group %in% c(G_HR, G_CTRL)) %>% arrange(group, strain, side)
kb_rng <- range(kb_pm$delta, na.rm = TRUE); kb_py <- kb_rng[2] + diff(kb_rng) * 0.05
F4A2 <- ggplot(kb_pm, aes(group, delta)) +
  geom_boxplot(aes(fill = group, alpha = side), color = "grey35",
               width = .6, outlier.shape = NA, position = position_dodge(.7), linewidth = .3) +
  geom_point(aes(color = group, shape = side),
             position = position_jitterdodge(.1, dodge.width = .7, seed = 1), size = 1.8, alpha = .9) +
  ggpubr::geom_pwc(data = kb_pair, aes(group = side), method = "wilcox_test",
           method.args = list(paired = TRUE), label = "p.format", p.adjust.method = "none",
           y.position = kb_py, tip.length = 0.02, size = 0.4, label.size = 3.0,
           vjust = -0.2, hide.ns = FALSE) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_color_manual(values = group_colors, guide = "none") +
  scale_alpha_manual(values = c(parent = .35, clone = .85), breaks = side_levels,
                     labels = side_labels[side_levels], name = "Group") +
  scale_shape_manual(values = c(parent = 1, clone = 16), breaks = side_levels,
                     labels = side_labels[side_levels], name = "Group") +
  coord_cartesian(ylim = c(kb_rng[1] - diff(kb_rng) * 0.05, kb_rng[2] + diff(kb_rng) * 0.22), clip = "off") +
  labs(x = NULL, y = "KB delta (mm)") +
  theme_f4 + theme(plot.margin = margin(5.5, 10, 5.5, 5.5))
panelA <- F4A2


## 6. Panel B: paired PAP AUC (F4B)
pap_auc_g <- add_grp(pap_auc)
F4B <- ggplot(pap_auc_g, aes(side, pap_auc_common, group = strain, color = group)) +
  geom_line(linewidth = .5, alpha = .7) + geom_point(size = 2) +
  geom_text_repel(data = filter(pap_auc_g, side == "clone"), aes(label = strain),
                  size = 2.4, nudge_x = .15, show.legend = FALSE, seed = 1,
                  max.overlaps = Inf, max.time = 2) +
  scale_color_manual(values = group_colors) + facet_wrap(~ group, nrow = 1) +
  scale_x_discrete(labels = side_labels) +
  labs(x = NULL, y = "PAP-AUC (0-32, relative growth)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))
panelB <- F4B


## 7. Panel C: paired SMG (F4D)
# Per-isolate mean of evaluable repeats; HR/Control paired Wilcoxon tests only.
smg_pm  <- add_grp(smg_rep) %>%
  group_by(group, strain, side) %>%
  summarise(SMG = if (any(is.finite(SMG))) mean(SMG, na.rm = TRUE) else NA_real_, .groups = "drop") %>%
  filter(is.finite(SMG))
smg_pair <- smg_pm %>% group_by(strain) %>% filter(all(side_levels %in% side)) %>% ungroup() %>%
  filter(group %in% c(G_HR, G_CTRL)) %>% arrange(group, strain, side)
smg_rng <- range(smg_pm$SMG, na.rm = TRUE); smg_py <- smg_rng[2] + diff(smg_rng) * 0.05
F4D <- ggplot(smg_pm, aes(group, SMG)) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = .4) +
  geom_boxplot(aes(fill = group, alpha = side), color = "grey35",
               width = .6, outlier.shape = NA, position = position_dodge(.7), linewidth = .3) +
  geom_point(aes(color = group, shape = side),
             position = position_jitterdodge(.1, dodge.width = .7, seed = 1), size = 1.8, alpha = .9) +
  ggpubr::geom_pwc(data = smg_pair, aes(group = side), method = "wilcox_test",
           method.args = list(paired = TRUE), label = "p.format", p.adjust.method = "none",
           y.position = smg_py, tip.length = 0.02, size = 0.4, label.size = 3.0,
           vjust = -0.2, hide.ns = FALSE) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_color_manual(values = group_colors, guide = "none") +
  scale_alpha_manual(values = c(parent = .35, clone = .85), breaks = side_levels,
                     labels = side_labels, name = "Group") +
  scale_shape_manual(values = c(parent = 1, clone = 16), breaks = side_levels,
                     labels = side_labels, name = "Group") +
  coord_cartesian(ylim = c(min(0, smg_rng[1]) - diff(smg_rng) * 0.05,
                           smg_rng[2] + diff(smg_rng) * 0.22), clip = "off") +
  labs(x = NULL, y = expression(SMG[72*h])) +
  theme_f4 + theme(plot.margin = margin(5.5, 10, 5.5, 5.5))
panelC <- F4D


## 8. Panel D: phenotype-shift magnitude (F4E)
disp <- master_out %>% drop_na(d_kb, d_pap) %>%
  mutate(z_kb = d_kb / sd(d_kb, na.rm = TRUE), z_pap = d_pap / sd(d_pap, na.rm = TRUE),
         shift = sqrt(z_kb^2 + z_pap^2)) %>%
  arrange(shift) %>% mutate(strain = factor(strain, levels = strain))
F4E <- ggplot(disp, aes(shift, strain, color = group)) +
  geom_segment(aes(x = 0, xend = shift, yend = strain), linewidth = .5, alpha = .55) +
  geom_point(size = 2.6) +
  scale_color_manual(values = group_colors, name = "group") +
  labs(x = "phenotype-shift magnitude\n(standardized KB + PAP displacement)", y = NULL) +
  theme_f4
panelD <- F4E


## 9. Panel E: WGS heterozygosity and k-mer ploidy
mw <- read_tsv(master_tsv, show_col_types = FALSE)
require_columns(mw, c("sample", "strain_id", "batch", "figure", "group",
                      "gate_n_het_pf", "gate_kmer_ploidy"), basename(master_tsv))
ren_group <- function(g) {
  factor(case_when(
    str_to_lower(as.character(g)) == "hr" ~ G_HR,
    str_to_lower(as.character(g)) == "control" ~ G_CTRL,
    str_to_lower(as.character(g)) == "r" ~ G_R,
    str_to_lower(as.character(g)) %in% c("standard", "h99") ~ G_H99,
    TRUE ~ as.character(g)
  ), levels = group_levels)
}
ploidy_lv <- c("haploid", "heterozygous_diploid", "het_low_confidence", "unknown")
ploidy_sh <- c(haploid = 21, heterozygous_diploid = 23,
               het_low_confidence = 4, unknown = 4)
clone4 <- mw %>% filter(figure == "Fig4") %>%
  transmute(
    strain = as.character(strain_id), group = ren_group(group),
    n_het = gate_n_het_pf,
    ploidy = factor(gate_kmer_ploidy, levels = ploidy_lv)
  )
parent4 <- mw %>% filter(batch == "parent", strain_id %in% clone4$strain) %>%
  mutate(strain = as.character(strain_id)) %>%
  filter(strain != "169" | sample == "Cn169b") %>%
  transmute(strain, parent_sample = sample, n_het_parent = gate_n_het_pf)
require_unique(clone4, "strain", "Selected WGS")
require_unique(parent4, "strain", "Parent WGS")
if (!setequal(clone4$strain, meta$strain) || !setequal(parent4$strain, meta$strain)) {
  stop("WGS measurements do not cover the phenotype isolate set on both sides.")
}
fig4w <- clone4 %>% left_join(parent4, by = "strain")
if (any(!is.finite(fig4w$n_het)) || any(!is.finite(fig4w$n_het_parent)) ||
    any(fig4w$n_het <= 0) || any(fig4w$n_het_parent <= 0) || anyNA(fig4w$ploidy)) {
  stop("Missing or invalid WGS heterozygosity/ploidy input.")
}
strain_ord <- as.character(strain_order)
fig4w <- fig4w %>% mutate(strain = factor(strain, levels = rev(strain_ord)))
het_rng <- range(c(fig4w$n_het_parent, fig4w$n_het), na.rm = TRUE)
band <- tibble(x = 10^seq(log10(het_rng[1] * .7), log10(het_rng[2] * 1.3), length.out = 60)) %>%
  mutate(lo = x / 2, hi = x * 2)
F4G <- ggplot(fig4w, aes(n_het_parent, n_het)) +
  geom_ribbon(data = band, aes(x = x, ymin = lo, ymax = hi), inherit.aes = FALSE,
              fill = "grey50", alpha = .12) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
  geom_point(aes(colour = group, fill = group, shape = ploidy), size = 2.0, stroke = 1) +
  geom_text_repel(aes(label = strain, colour = group), size = 2.3,
                  max.overlaps = Inf, show.legend = FALSE, seed = 1) +
  scale_x_log10(labels = label_number()) + scale_y_log10(labels = label_number()) +
  scale_colour_manual(values = group_colors, name = "group") +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_shape_manual(values = ploidy_sh, drop = TRUE, name = "k-mer ploidy") +
  coord_fixed() +
  annotate("text", x = het_rng[2], y = het_rng[2] * .45,
           label = "on diagonal = no change", hjust = 1, size = 2.8, colour = "grey45") +
  labs(x = "parent het sites (popFreq, log)", y = "Selected het sites (popFreq, log)") +
  theme_f4 + theme(legend.position = "right", legend.box = "vertical")
panelE <- F4G


## 10. RNA-seq configuration
# Test set: clinical isolates with a paired Selected/Monoclonal transcriptome.
# Cn160 is excluded because only 41-45% of its reads map to H99
clinical_strains <- c(
  "101", "106", "126", "128", "162", "208", "213", "215", "217",  # HR
  "147", "150", "169", "174", "225", "226"                        # Control
)
positive_control_strains <- "H99"
display_strains <- c(clinical_strains, positive_control_strains)

contrast_batches <- c("01", "02")    # Selected (01) vs Monoclonal (02); ref = 02
drop_samples     <- character(0)
min_count        <- 10               # keep a gene with >= min_count in >= min_n samples
min_n            <- 3
padj_cut         <- 0.05
lfc_cut          <- 1

# Authoritative symbol/CNAG map from the FungiDB GO annotation. The count matrix,
# the GTF and the annotation share the H99 CNA3 identifier space, so this is a
# pure lookup with no re-mapping.
.read_gaf_sym <- function(f){
  ln <- readLines(gzfile(f)); ln <- ln[!startsWith(ln, "!")]
  x  <- strsplit(ln, "\t", fixed = TRUE)
  tibble(gene_id = vapply(x, `[`, "", 2L), symbol = vapply(x, `[`, "", 3L)) %>%
    filter(startsWith(gene_id, "CNAG"), symbol != "", symbol != gene_id)
}
cnag2sym <- map_dfr(gaf_files, .read_gaf_sym) %>%
  distinct(gene_id, .keep_all = TRUE)
if (!nrow(cnag2sym)) stop("No symbol/CNAG pairs parsed from the GO annotation files.")

# Three-tier target panel.
target_panel <- tibble::tribble(
  ~symbol,  ~tier,
  "ERG11","dosage",  "AFR1","dosage",  "ERG2","dosage",
  "CRZ1","output",   "PMC1","output",  "VCX1","output",  "PMR1","output",
  "CHS5","output",   "CHS6","output",  "CHS7","output",
  "CDA2","output",   "KRE6","output",  "FKS1","output",
  "HSP78","output",  "TPS1","output",  "TPS2","output",  "NTH1","output",
  "CNA1","context",  "CNB1","context", "HSP90","context",
  "TOR1","context",  "YPK1","context", "PKC1","context",
  "MPK1","context",  "HOG1","context"
)
.manual <- c(HSP90 = "CNAG_06150", HSP78 = "CNAG_03347", CHS6  = "CNAG_06487",
             CRZ1  = "CNAG_00156", CDA2  = "CNAG_01230")

target_panel <- target_panel %>%
  dplyr::left_join(cnag2sym %>% dplyr::filter(!symbol %in% names(.manual)), by = "symbol")
dup <- target_panel %>% dplyr::count(symbol) %>% dplyr::filter(n > 1)
if (nrow(dup)) stop("Multi-mapped annotation symbol: ", paste(dup$symbol, collapse = ", "))
target_panel <- target_panel %>%
  dplyr::mutate(gene_id = dplyr::coalesce(unname(.manual[symbol]), gene_id),
                tier = factor(tier, levels = c("dosage","output","context")))
stopifnot(!anyNA(target_panel$gene_id))
.expected <- c(HSP78 = "CNAG_03347", HSP90 = "CNAG_06150", CHS6 = "CNAG_06487")
stopifnot(nrow(target_panel) == 25L,
          !anyDuplicated(target_panel$symbol), !anyDuplicated(target_panel$gene_id),
          identical(unname(target_panel$gene_id[match(names(.expected), target_panel$symbol)]),
                    unname(.expected)),
          !("HSP104" %in% target_panel$symbol), !("CHS4" %in% target_panel$symbol))

target_genes <- setNames(target_panel$gene_id, target_panel$symbol)
target_genes[is.na(target_genes)] <- .manual[names(target_genes)[is.na(target_genes)]]
stopifnot(!anyNA(target_genes))


## 11. RNA-seq differential expression (Selected vs Monoclonal)
# gene_id -> product description, used to annotate results and confirm targets.
gtf_lines <- readLines(gtf_path)
gtf_lines <- gtf_lines[!startsWith(gtf_lines, "#")]
.gid  <- sub('.*gene_id "([^"]+)".*', "\\1", gtf_lines)
.prod <- ifelse(grepl('product "', gtf_lines),
                sub('.*product "([^"]+)".*', "\\1", gtf_lines), NA_character_)
gene_product <- tibble(gene_id = .gid, product = .prod) %>%
  filter(!is.na(product)) %>% distinct(gene_id, .keep_all = TRUE)

fc  <- read.delim(count_tsv, comment.char = "#", check.names = FALSE)
require_columns(fc, "Geneid", basename(count_tsv))
mat <- as.matrix(fc[, -(1:6)]); rownames(mat) <- fc$Geneid
colnames(mat) <- sub(".*/(Cn[^/]+)\\.Aligned.*", "\\1", colnames(mat))

samp <- tibble(sample = colnames(mat)) %>%
  tidyr::extract(sample, c("strain", "batch"),
                 "^Cn(\\d{3}|H99)(\\d{2})$", remove = FALSE)

keep <- samp %>%
  filter(strain %in% clinical_strains,
         batch  %in% contrast_batches,
         !sample %in% drop_samples)
if (!nrow(keep)) stop("No RNA-seq samples match the Figure 3 contrast.")

mat_use <- mat[, keep$sample]
col <- keep %>%
  mutate(strain = factor(strain, levels = intersect(clinical_strains, strain)),
         batch  = relevel(factor(batch), ref = contrast_batches[2])) %>%
  as.data.frame()
rownames(col) <- col$sample

message(sprintf("RNA-seq Figure 3: %d samples | %d strains | per-batch n = %s",
                nrow(col), nlevels(col$strain),
                paste(names(table(col$batch)), table(col$batch), sep = ":", collapse = " ")))
.n_per_strain <- table(keep$strain)
unpaired <- names(.n_per_strain)[.n_per_strain < 2]
if (length(unpaired)) message("  NOTE unpaired strain(s) (only one batch present): ",
                              paste(unpaired, collapse = ", "))

dds <- DESeqDataSetFromMatrix(round(mat_use), col, design = ~ strain + batch)
.keep <- rowSums(counts(dds) >= min_count) >= min_n
dds <- dds[.keep | rownames(dds) %in% target_genes, ]
dds <- DESeq(dds)

coef_name <- grep("^batch", resultsNames(dds), value = TRUE)   # "batch_01_vs_02"
if (length(coef_name) != 1L) stop("Could not resolve a single batch coefficient.")
res <- results(dds, name = coef_name, alpha = padj_cut)
# Shrink LFC for ranking and plotting; apeglm preferred, normal as a fallback.
res_s <- tryCatch(
  lfcShrink(dds, coef = coef_name, type = "apeglm"),
  error = function(e) { message("apeglm unavailable -> type='normal'")
    lfcShrink(dds, coef = coef_name, type = "normal") })

res_tab <- as_tibble(res_s, rownames = "gene_id") %>%
  left_join(cnag2sym,     by = "gene_id") %>%
  left_join(gene_product, by = "gene_id") %>%
  mutate(symbol = ifelse(is.na(symbol) & gene_id %in% target_genes,
                         names(target_genes)[match(gene_id, target_genes)], symbol)) %>%
  arrange(padj)

n_up   <- sum(res_tab$padj < padj_cut & res_tab$log2FoldChange >  lfc_cut, na.rm = TRUE)
n_down <- sum(res_tab$padj < padj_cut & res_tab$log2FoldChange < -lfc_cut, na.rm = TRUE)
message(sprintf("DE genes (padj<%.2f & |LFC|>%g): %d up / %d down",
                padj_cut, lfc_cut, n_up, n_down))

# Display dataset: clinical isolates plus the H99 positive control. Size factors only.
keep_display <- samp %>%
  filter(strain %in% display_strains, batch %in% contrast_batches, !sample %in% drop_samples)
col_display <- keep_display %>%
  mutate(strain = factor(strain, levels = display_strains),
         batch  = relevel(factor(batch), ref = contrast_batches[2])) %>%
  as.data.frame()
rownames(col_display) <- col_display$sample
dds_display <- DESeqDataSetFromMatrix(
  countData = round(mat[, keep_display$sample, drop = FALSE]),
  colData   = col_display, design = ~ strain + batch)
.keep_display <- rowSums(counts(dds_display) >= min_count) >= min_n
dds_display <- dds_display[.keep_display | rownames(dds_display) %in% target_genes, ]
dds_display <- estimateSizeFactors(dds_display)


## 12. Panel F: PCA of variance-stabilized counts
set.seed(123)
vsd <- vst(dds, blind = FALSE)
batch_labs <- c("02" = "Monoclonal", "01" = "Selected", "03" = "Passage")
batch_cols <- c("Monoclonal" = "#4F91C7", "Selected" = "#D66A78", "Passage" = "#8A6BB8")

pca <- plotPCA(vsd, intgroup = c("batch", "strain"), returnData = TRUE) %>%
  mutate(batch = factor(as.character(batch), levels = c("02", "01", "03")),
         batch_stage = factor(batch_labs[as.character(batch)],
                              levels = c("Monoclonal", "Selected", "Passage")),
         strain = as.character(strain))
pv  <- round(100 * attr(pca, "percentVar"))
pca_lab <- pca %>% filter(!is.na(batch_stage)) %>%
  group_by(batch_stage) %>%
  summarise(PC1 = median(PC1, na.rm = TRUE), PC2 = median(PC2, na.rm = TRUE), .groups = "drop") %>%
  mutate(label = as.character(batch_stage))
p_pca <- ggplot(pca, aes(PC1, PC2, colour = batch_stage)) +
  stat_ellipse(aes(group = batch_stage, fill = batch_stage), geom = "polygon",
               type = "norm", level = 0.80, alpha = 0.12,
               linewidth = 0.4, show.legend = FALSE) +
  stat_ellipse(aes(group = batch_stage), type = "norm",
               level = 0.80, linewidth = 0.5, show.legend = FALSE) +
  geom_point(size = 2.5) +
  ggrepel::geom_text_repel(aes(label = strain), size = 2.6, show.legend = FALSE,
                           max.overlaps = Inf, min.segment.length = 0, seed = 123) +
  ggrepel::geom_label_repel(data = pca_lab, aes(label = label, fill = batch_stage),
                            colour = "white", size = 3.0, label.size = 0.25,
                            alpha = 0.75, show.legend = FALSE, seed = 123) +
  scale_colour_manual(values = batch_cols, drop = TRUE) +
  scale_fill_manual(values = batch_cols, drop = TRUE) +
  labs(x = sprintf("PC1 (%d%%)", pv[1]), y = sprintf("PC2 (%d%%)", pv[2]),
       colour = "batch") +
  guides(fill = "none") +
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.9),
    axis.ticks = element_line(linewidth = 0.6),
    axis.title.x = element_text(size = 12, face = "plain"),
    axis.title.y = element_text(size = 12, face = "plain")
  )
panelF <- p_pca


## 13. Panel G: volcano plot
set.seed(123)
vol_plot <- res_tab %>%
  mutate(neglog10_padj = -log10(pmax(padj, .Machine$double.xmin)),
         sig_class = case_when(
           padj < padj_cut & log2FoldChange <= -lfc_cut ~ "down",
           padj < padj_cut & log2FoldChange >=  lfc_cut ~ "up",
           TRUE ~ "stable"),
         sig_class = factor(sig_class, levels = c("down", "stable", "up")))
label_df <- vol_plot %>%
  filter(sig_class %in% c("up", "down"), !is.na(symbol), is.finite(neglog10_padj)) %>%
  group_by(sig_class) %>% slice_max(order_by = neglog10_padj, n = 12, with_ties = FALSE) %>%
  ungroup()
class_note <- sprintf(
  "Classification: up = padj < %.2g and log2FC >= %.2f; down = padj < %.2g and log2FC <= -%.2f; stable = others.",
  padj_cut, lfc_cut, padj_cut, lfc_cut)

p_vol <- ggplot(vol_plot, aes(neglog10_padj, log2FoldChange)) +
  geom_point(aes(colour = sig_class), size = 1.2, alpha = 0.70) +
  geom_vline(xintercept = -log10(padj_cut), linetype = "dotted", colour = "grey35", linewidth = 0.55) +
  geom_hline(yintercept = c(-lfc_cut, lfc_cut), linetype = "dotted", colour = "grey35", linewidth = 0.55) +
  ggrepel::geom_label_repel(
    data = label_df, aes(label = symbol, colour = sig_class),
    fill = "white", size = 3.1, label.size = 0.25,
    label.padding = unit(0.12, "lines"), box.padding = 0.45,
    point.padding = 0.25, segment.size = 0.3, segment.alpha = 0.6,
    min.segment.length = 0, max.overlaps = Inf, max.time = 2,
    seed = 123, show.legend = FALSE) +
  scale_colour_manual(
    values = c("down" = "#2B5C87", "stable" = "grey75", "up" = "#A63446"),
    guide = "none") +
  labs(x = "-log10 (padj)", y = "log2 (FoldChange)",
       title = "Selected vs Monoclonal", caption = class_note) +
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "grey25", fill = NA, linewidth = 0.9),
    panel.grid.major = element_line(colour = "grey88", linewidth = 0.35),
    panel.grid.minor = element_line(colour = "grey92", linewidth = 0.25),
    axis.title.x = element_text(size = 12, face = "plain"),
    axis.title.y = element_text(size = 12, face = "plain"),
    axis.text = element_text(size = 12, colour = "grey25"),
    axis.ticks = element_line(linewidth = 0.8),
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold"),
    plot.caption = element_text(size = 9, colour = "grey25", hjust = 0))
panelG <- p_vol


## 14. Panel H: within-strain target-gene heatmap (ht)
tg_tab <- res_tab %>% filter(gene_id %in% target_genes) %>%
  mutate(symbol = names(target_genes)[match(gene_id, target_genes)]) %>%
  select(symbol, gene_id, baseMean, log2FoldChange, lfcSE, pvalue, padj, product)
if (nrow(tg_tab) != nrow(target_panel)) {
  stop("Target panel genes missing from the differential-expression table: ",
       paste(setdiff(target_panel$symbol, tg_tab$symbol), collapse = ", "))
}

tier_labels <- c(dosage  = "Azole resistance (dosage)",
                 output  = "Stress / cell-wall output",
                 context = "Upstream kinase / sensor")
tier_fill <- c(dosage = "#EAD9C9", output = "#DCE6DC", context = "#D6DEE8")
tier_block_labels <- c(dosage  = "Azole\nresistance",
                       output  = "Stress /\ncell-wall\noutput",
                       context = "Upstream\nkinase /\nsensor")
star <- function(p) dplyr::case_when(is.na(p) ~ "", p < 1e-3 ~ "***",
                                     p < 1e-2 ~ "**", p < 0.05 ~ "*", TRUE ~ "")

grp_map <- tibble(strain = display_strains) %>%
  left_join(meta %>% transmute(strain = as.character(strain), group), by = "strain") %>%
  mutate(group = factor(group, levels = group_levels))

present <- intersect(unname(target_genes), rownames(dds_display))
nc <- counts(dds_display, normalized = TRUE)[present, , drop = FALSE]

fc_wide <- as_tibble(nc, rownames = "gene_id") %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "n") %>%
  left_join(keep_display, by = "sample") %>%
  transmute(gene_id, strain, batch, n) %>%
  pivot_wider(names_from = batch, values_from = n) %>%
  mutate(l2fc = log2((`01` + 1) / (`02` + 1))) %>%
  left_join(dplyr::select(target_panel, gene_id, symbol), by = "gene_id")

M <- fc_wide %>% dplyr::select(symbol, strain, l2fc) %>%
  pivot_wider(names_from = strain, values_from = l2fc) %>%
  tibble::column_to_rownames("symbol") %>% as.matrix()

row_meta <- tg_tab %>%
  transmute(symbol, tier = target_panel$tier[match(symbol, target_panel$symbol)],
            lfc = log2FoldChange, padj) %>%
  mutate(tier = factor(tier, levels = c("dosage", "output", "context"))) %>%
  arrange(tier, desc(lfc))
M <- M[row_meta$symbol, , drop = FALSE]
row_split <- factor(tier_labels[as.character(row_meta$tier)],
                    levels = tier_labels[c("dosage", "output", "context")])

col_grp <- droplevels(factor(grp_map$group[match(colnames(M), grp_map$strain)],
                             levels = group_levels))

lim     <- max(abs(M), na.rm = TRUE)
col_fun <- colorRamp2(c(-lim, 0, lim), c("#2B5C87", "#FBFEF9", "#A63446"))

top_ann <- HeatmapAnnotation(
  Group = col_grp,
  col = list(Group = group_colors[levels(col_grp)]),
  annotation_label = "Group",
  annotation_name_gp = gpar(fontsize = 8),
  simple_anno_size = unit(3.5, "mm"))

ra <- rowAnnotation(
  "log2FC" = anno_barplot(
    row_meta$lfc, bar_width = 0.7, axis_param = list(gp = gpar(fontsize = 6)),
    gp = gpar(fill = ifelse(row_meta$lfc > 0, "#A63446", "#2B5C87"), col = NA)),
  width = unit(1.6, "cm"), annotation_name_gp = gpar(fontsize = 8.8, fontface = "bold"))
la <- rowAnnotation(
  tier = anno_block(
    gp        = gpar(fill = tier_fill[c("dosage", "output", "context")],
                     col = "grey40", lwd = 0.8),
    labels    = tier_block_labels[c("dosage", "output", "context")],
    labels_gp = gpar(fontsize = 8, fontface = "bold", col = "grey15"),
    labels_rot = 0),
  width = unit(1.7, "cm"),
  show_annotation_name = FALSE)

ht <- Heatmap(
  M, name = "Selected\n   vs\nMonoclonal", col = col_fun, na_col = "grey90",
  cluster_rows = FALSE, row_order = row_meta$symbol, cluster_column_slices = FALSE,
  cluster_columns = TRUE, clustering_distance_columns = "euclidean",
  clustering_method_columns = "ward.D2", show_column_dend = TRUE,
  column_dend_height = unit(6, "mm"),
  row_split = row_split, column_split = col_grp, row_title = NULL,
  row_labels = paste0(row_meta$symbol, star(row_meta$padj)),
  rect_gp = gpar(col = "grey92", lwd = 0.4),
  row_title_rot = 0, row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  row_gap = unit(2, "mm"), column_gap = unit(2, "mm"),
  row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 7),
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  left_annotation = la, top_annotation = top_ann, right_annotation = ra,
  heatmap_legend_param = list(legend_height = unit(3, "cm")))
panelH <- ht


## 15. Derived table for downstream figures
dir.create(DERIVED_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(DERIVED_DIR)) stop("Cannot create derived directory: ", DERIVED_DIR)
write_csv(res_tab, file.path(DERIVED_DIR, "fig3_DE_batch01_vs_02.csv"))
write_csv(tg_tab, file.path(DERIVED_DIR, "fig3_targets_DE.csv"))
message("Figure 3 derived tables written to: ", normalizePath(DERIVED_DIR))


## 16. Calculated checks and PDF output
paired_stats <- function(d, value, endpoint) {
  w <- d %>% filter(group %in% c(G_HR, G_CTRL)) %>%
    transmute(group, strain, side, value = .data[[value]]) %>%
    pivot_wider(names_from = side, values_from = value) %>%
    filter(is.finite(parent), is.finite(clone))
  w %>% group_by(group) %>% summarise(
    endpoint = endpoint, n_pairs = n(),
    p = if (n() >= 3) wilcox.test(parent, clone, paired = TRUE)$p.value else NA_real_,
    .groups = "drop"
  ) %>% select(endpoint, group, n_pairs, p)
}
figure3_stats <- bind_rows(
  paired_stats(kb_pm, "delta", "A: KB delta"),
  paired_stats(smg_pm, "SMG", "C: SMG")
)
smg_n_pairs <- smg_pm %>% count(strain) %>% filter(n == 2L) %>% nrow()
message(
  "Figure 3: ", n_distinct(kb_pm$strain), " KB pairs; ",
  n_distinct(pap_auc_g$strain), " PAP pairs; ",
  nrow(smg_pm), " evaluable isolate/side SMG summaries (", smg_n_pairs, " complete pairs); ",
  nrow(disp), " displacement values; ", nrow(fig4w), " WGS pairs; ",
  nrow(col), " RNA-seq samples (", n_up, " up / ", n_down, " down)."
)
print(figure3_stats, n = Inf)

if (EXPORT_PANELS) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(OUTPUT_DIR)) stop("Cannot create output directory: ", OUTPUT_DIR)
  if (!capabilities("cairo")) stop("PDF output requires Cairo support in R.")
  panel_plots <- list(
    Fig3A = panelA, Fig3B = panelB, Fig3C = panelC, Fig3D = panelD,
    Fig3E = panelE, Fig3F = panelF, Fig3G = panelG
  )
  panel_sizes <- list(
    Fig3A = c(3.9, 3.8), Fig3B = c(5.3, 3.9), Fig3C = c(3.7, 3.5), Fig3D = c(4.2, 4.9),
    Fig3E = c(5.0, 4.2), Fig3F = c(5.1, 4.2), Fig3G = c(6.0, 5.1)
  )
  iwalk(panel_plots, function(p, name) {
    size <- panel_sizes[[name]]
    ggsave(
      file.path(OUTPUT_DIR, paste0(name, ".pdf")), p,
      width = size[1], height = size[2], units = "in", device = cairo_pdf
    )
  })
  # Panel H is a ComplexHeatmap object and is drawn straight to a device.
  cairo_pdf(file.path(OUTPUT_DIR, "Fig3H.pdf"), width = 6.4, height = 5)
  ComplexHeatmap::draw(panelH)
  dev.off()
  message("Figure 3 PDFs written to: ", normalizePath(OUTPUT_DIR))
}
