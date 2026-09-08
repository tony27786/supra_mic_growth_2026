# Figure 2 plot code
#
# Inputs: ../SourceData. Default PDF output: ../Figure2_panels.
# Set FIG2_EXPORT=false to calculate without writing PDFs.
# FIG2_OUTPUT_DIR can redirect PDF output; no input files are modified.

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
  found <- candidates[file.exists(file.path(candidates, "fig2_submit.R"))]
  if (length(found) != 1L) {
    stop("Run fig2_submit.R with source() or Rscript so its directory can be resolved.")
  }
  normalizePath(found, mustWork = TRUE)
}
SCRIPT_DIR <- get_script_dir()
SOURCE_DIR <- normalizePath(file.path(SCRIPT_DIR, "..", "SourceData"), mustWork = TRUE)
OUTPUT_DIR <- Sys.getenv("FIG2_OUTPUT_DIR", file.path(SCRIPT_DIR, "..", "Figure2_panels"))
EXPORT_PANELS <- tolower(Sys.getenv("FIG2_EXPORT", "true")) %in% c("1", "true", "yes")

kb_xlsx <- file.path(SOURCE_DIR, "kbpap_SourceData.xlsx")
smg_xlsx <- file.path(SOURCE_DIR, "smg_SourceData.xlsx")
pap_parent_csv <- file.path(SOURCE_DIR, "parent_pap_feature_results.csv")
pap_clone_csv <- file.path(SOURCE_DIR, "monoclonal_pap_feature_results.csv")
master_tsv <- file.path(SOURCE_DIR, "master_wgs_table.tsv")
input_files <- c(kb_xlsx, smg_xlsx, pap_parent_csv, pap_clone_csv, master_tsv)
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
CONC   <- c(128,64,32,16,8,4,2,1,0.5)

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
side_labels <- c(parent = "Parent", clone = "Monoclonal")
theme_f3 <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"), legend.position = "top")
canon <- function(s){ s <- as.character(s); ifelse(nchar(s) > 3, substr(s, 1, nchar(s) - 2), s) }
trapz <- function(x, y){ o <- order(x); x <- x[o]; y <- y[o]; sum(diff(x) * (head(y,-1)+tail(y,-1))/2) }


## 1. KB measurements and isolate metadata
ci <- read_excel(kb_xlsx, sheet = "Clinical_Isolates") %>%
  mutate(strain = as.character(ID))
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
x2 <- read_excel(kb_xlsx, sheet = "xxx02") %>% mutate(strain = as.character(ID))
kb_clone <- x2 %>%
  select(strain, re_1_delta, re_2_delta, re_3_delta) %>%
  pivot_longer(-strain, names_to = "rep", values_to = "delta") %>% mutate(side = "clone")
kb_long <- bind_rows(kb_parent, kb_clone) %>% filter(!is.na(delta))
require_unique(meta, "strain", "Clinical_Isolates")
require_unique(x2, "strain", "xxx02")
if (!setequal(meta$strain, x2$strain)) {
  stop("Parent and Monoclonal KB isolate IDs do not match.")
}
if (anyNA(meta$group)) stop("Unrecognized group in Clinical_Isolates.")



## 2. PAP from image-analysis project measurements
# Reproduce the raw_mean portion of PAPDual/PAPSingle:
# Parent H/L plates overlap at FLC 0 and 8; average raw_mean across both plates
# at each shared concentration BEFORE dividing by the averaged FLC0 reference.
# Monoclonal uses a single plate with wells 0, 8, 12, 16, 24, 32.
# Background images are excluded as in the original normalization = "0" calls.
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
  read_smg_sheet("xxx02", c("24h","48h","72h","96h"),  "clone", is_clone = TRUE)) %>%
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
# Panel A uses medians; the displacement in Panel D uses the original mean KB
# readings. SMG is the mean of evaluable repeats, separately for each side.
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


## 5. Panel A: outer diameter and KB delta
pull_outer <- function(d, reps, s)
  d %>%
  select(strain, all_of(paste0(reps, "_outer_flc"))) %>%
  pivot_longer(-strain, names_to = "rep", values_to = "outer") %>%
  mutate(rep = str_remove(rep, "_outer_flc$"), side = s)

kb_out <- bind_rows(
  pull_outer(ci, c("ori", paste0("re_", 1:5)), "parent"),
  pull_outer(x2, paste0("re_", 1:3),           "clone")
) %>%
  filter(!is.na(outer)) %>%
  left_join(meta %>% select(strain, group), by = "strain") %>%
  mutate(group = factor(group, levels = group_levels),
         side  = factor(side,  levels = side_levels))

kb_out_pm <- kb_out %>%
  group_by(group, strain, side) %>%
  summarise(outer = median(outer, na.rm = TRUE), .groups = "drop") %>%
  group_by(strain) %>% filter(all(side_levels %in% side)) %>% ungroup()
# Align rows by isolate for within-group paired tests.
kb_out_pm <- kb_out_pm %>% arrange(group, strain, side)
kb_out_rng <- range(kb_out_pm$outer, na.rm = TRUE)
kb_out_y   <- kb_out_rng[2] + diff(kb_out_rng) * .08

F3A0 <- ggplot(kb_out_pm, aes(group, outer)) +
  geom_boxplot(aes(fill = group, alpha = side), color = "grey35",
               width = .6, outlier.shape = NA, position = position_dodge(.7), linewidth = .3) +
  geom_point(aes(color = group, shape = side),
             position = position_jitterdodge(.1, dodge.width = .7, seed = 1), size = 1.8, alpha = .9) +
  ggpubr::geom_pwc(
    data = kb_out_pm %>% filter(group %in% c(G_HR, G_CTRL)),
    aes(group = side), method = "wilcox_test",
    method.args = list(paired = TRUE), label = "p.format",
    p.adjust.method = "none", y.position = kb_out_y,
    tip.length = 0.02, size = 0.4,
    label.size = 3.0, vjust = -0.2, hide.ns = FALSE) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_color_manual(values = group_colors, guide = "none") +
  scale_alpha_manual(values = c(parent = .35, clone = .85), breaks = side_levels, 
                     labels = side_labels[side_levels], name = "Group") +
  scale_shape_manual(values = c(parent = 1, clone = 16), breaks = side_levels, 
                     labels = side_labels[side_levels], name = "Group") +
  coord_cartesian(ylim = c(kb_out_rng[1] - diff(kb_out_rng) * 0.05,
      kb_out_rng[2] + diff(kb_out_rng) * 0.22), clip = "off") +
  labs(x = NULL, y = "Outer zone diameter (mm)") +
  theme_f3 + theme(plot.margin = margin(5.5, 10, 5.5, 5.5))
kb_pm  <- add_grp(kb_long) %>%
  group_by(group, strain, side) %>% summarise(delta = median(delta, na.rm = TRUE), .groups = "drop")
kb_pair <- kb_pm %>% group_by(strain) %>% filter(all(side_levels %in% side)) %>% ungroup() %>%
  filter(group %in% c(G_HR, G_CTRL)) %>% arrange(group, strain, side)
kb_rng <- range(kb_pm$delta, na.rm = TRUE); kb_py <- kb_rng[2] + diff(kb_rng) * 0.05
F3A2 <- ggplot(kb_pm, aes(group, delta)) +
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
  theme_f3 + theme(plot.margin = margin(5.5, 10, 5.5, 5.5))
panelA_outer <- F3A0
panelA_delta <- F3A2
panelA <- (F3A0 | F3A2) + plot_layout(guides = "collect") &
  theme(legend.position = "top")



## 6. Panel B: paired PAP AUC
pap_auc_g <- add_grp(pap_auc)
F3B <- ggplot(pap_auc_g, aes(side, pap_auc_common, group = strain, color = group)) +
  geom_line(linewidth = .5, alpha = .7) + geom_point(size = 2) +
  geom_text_repel(data = filter(pap_auc_g, side == "clone"), aes(label = strain),
                  size = 2.4, nudge_x = .3, direction = "y",
                  box.padding = .35, point.padding = .1, min.segment.length = 0,
                  max.overlaps = Inf, max.iter = 100000, max.time = 2,
                  show.legend = FALSE, seed = 1) +
  scale_color_manual(values = group_colors) + facet_wrap(~ group, nrow = 1) +
  scale_x_discrete(labels = c(parent = "parent", clone = "mono\n-clonal")) +
  labs(x = NULL, y = "PAP-AUC (0–32, relative growth)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))
panelB <- F3B



## 7. Panel C: paired SMG
# Per-isolate mean of evaluable repeats; HR/Control paired Wilcoxon tests only.
smg_pm  <- add_grp(smg_rep) %>%
  group_by(group, strain, side) %>%                     # [v2] mean, not median (sync with Fig 1)
  summarise(SMG = if (any(is.finite(SMG))) mean(SMG, na.rm = TRUE) else NA_real_, .groups = "drop") %>%
  filter(is.finite(SMG))                                # [v2] plan A: no MIC -> no SMG -> drop
smg_pair <- smg_pm %>% group_by(strain) %>% filter(all(side_levels %in% side)) %>% ungroup() %>%
  filter(group %in% c(G_HR, G_CTRL)) %>% arrange(group, strain, side)
smg_rng <- range(smg_pm$SMG, na.rm = TRUE); smg_py <- smg_rng[2] + diff(smg_rng) * 0.05
F3D <- ggplot(smg_pm, aes(group, SMG)) +
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
  scale_alpha_manual(values = c(parent = .35, clone = .85),
                     breaks = side_levels, labels = side_labels, name = "Group") +
  scale_shape_manual(values = c(parent = 1, clone = 16),
                     breaks = side_levels, labels = side_labels, name = "Group") +
  coord_cartesian(ylim = c(min(0, smg_rng[1]) - diff(smg_rng) * 0.05, smg_rng[2] + diff(smg_rng) * 0.22), clip = "off") +
  labs(x = NULL, y = expression(SMG[72*h])) +
  theme_f3 + theme(plot.margin = margin(5.5, 10, 5.5, 5.5))
panelC <- F3D



## 8. Panel D: phenotype-shift magnitude
disp <- master_out %>% drop_na(d_kb, d_pap) %>%
  mutate(z_kb = d_kb / sd(d_kb, na.rm = TRUE), z_pap = d_pap / sd(d_pap, na.rm = TRUE),
         shift = sqrt(z_kb^2 + z_pap^2)) %>%
  arrange(shift) %>% mutate(strain = factor(strain, levels = strain))
F3E <- ggplot(disp, aes(shift, strain, color = group)) +
  geom_segment(aes(x = 0, xend = shift, yend = strain), linewidth = .5, alpha = .55) +
  geom_point(size = 2.6) +
  scale_color_manual(values = group_colors, name = "group") +
  labs(x = "phenotype-shift magnitude\n(standardized KB + PAP displacement)", y = NULL) +
  theme_f3
panelD <- F3E



## 9. Panel E: WGS heterozygosity and k-mer ploidy (F3G)
# Fig3 is the original WGS table's designation for the Monoclonal samples.
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
clone3 <- mw %>% filter(figure == "Fig3") %>%
  transmute(
    strain = as.character(strain_id), group = ren_group(group),
    n_het = gate_n_het_pf,
    ploidy = factor(gate_kmer_ploidy, levels = ploidy_lv)
  )
parent3 <- mw %>% filter(batch == "parent", strain_id %in% clone3$strain) %>%
  mutate(strain = as.character(strain_id)) %>%
  filter(strain != "169" | sample == "Cn169b") %>%
  transmute(strain, parent_sample = sample, n_het_parent = gate_n_het_pf)
require_unique(clone3, "strain", "Monoclonal WGS")
require_unique(parent3, "strain", "Parent WGS")
if (!setequal(clone3$strain, meta$strain) || !setequal(parent3$strain, meta$strain)) {
  stop("WGS measurements do not cover the phenotype isolate set on both sides.")
}
fig3w <- clone3 %>% left_join(parent3, by = "strain")
if (any(!is.finite(fig3w$n_het)) || any(!is.finite(fig3w$n_het_parent)) ||
    any(fig3w$n_het <= 0) || any(fig3w$n_het_parent <= 0) || anyNA(fig3w$ploidy)) {
  stop("Missing or invalid WGS heterozygosity/ploidy input.")
}
strain_ord <- as.character(strain_order)
fig3w <- fig3w %>% mutate(strain = factor(strain, levels = rev(strain_ord)))
het_rng <- range(c(fig3w$n_het_parent, fig3w$n_het), na.rm = TRUE)
band <- tibble(x = 10^seq(log10(het_rng[1] * .7), log10(het_rng[2] * 1.3), length.out = 60)) %>%
  mutate(lo = x / 2, hi = x * 2)
F3G <- ggplot(fig3w, aes(n_het_parent, n_het)) +
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
  labs(x = "parent het sites (popFreq, log)", y = "clone het sites (popFreq, log)") +
  theme_f3 + theme(legend.position = "right", legend.box = "vertical")
panelE <- F3G



## 10. Calculated checks and PDF output
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
figure2_stats <- bind_rows(
  paired_stats(kb_out_pm, "outer", "A: outer diameter"),
  paired_stats(kb_pm, "delta", "A: KB delta"),
  paired_stats(smg_pm, "SMG", "C: SMG")
)
smg_n_pairs <- smg_pm %>% count(strain) %>% filter(n == 2L) %>% nrow()
message(
  "Figure 2: ", n_distinct(kb_pm$strain), " KB pairs; ",
  n_distinct(pap_auc_g$strain), " PAP pairs; ",
  nrow(smg_pm), " evaluable isolate/side SMG summaries (", smg_n_pairs, " complete pairs); ",
  nrow(disp), " displacement values; ", nrow(fig3w), " WGS pairs."
)
print(figure2_stats, n = Inf)

if (EXPORT_PANELS) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(OUTPUT_DIR)) stop("Cannot create output directory: ", OUTPUT_DIR)
  if (!capabilities("cairo")) stop("PDF output requires Cairo support in R.")
  panel_plots <- list(
    Fig2A_outer = panelA_outer, Fig2A_delta = panelA_delta, Fig2A = panelA,
    Fig2B = panelB, Fig2C = panelC, Fig2D = panelD, Fig2E = panelE
  )
  panel_sizes <- list(
    Fig2A_outer = c(3.9, 3.8), Fig2A_delta = c(3.9, 3.8), Fig2A = c(7.8, 3.8),
    Fig2B = c(5.3, 3.9), Fig2C = c(3.7, 3.5),
    Fig2D = c(4.2, 4.9), Fig2E = c(5.0, 4.2)
  )
  iwalk(panel_plots, function(p, name) {
    size <- panel_sizes[[name]]
    ggsave(
      file.path(OUTPUT_DIR, paste0(name, ".pdf")), p,
      width = size[1], height = size[2], units = "in", device = cairo_pdf
    )
  })
  message("Figure 2 PDFs written to: ", normalizePath(OUTPUT_DIR))
}
