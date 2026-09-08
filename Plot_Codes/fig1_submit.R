# Figure 1 plot code
#
# Inputs: ../SourceData. Default PDF output: ../Figure1_panels.
# Set FIG1_EXPORT=false to calculate without writing PDFs.
# FIG1_OUTPUT_DIR can redirect PDF output; no input files are modified.

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(ggrepel)

## -----------------------------------------------------------------------------
## Paths and constants
## -----------------------------------------------------------------------------
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)))
  }

  frames <- sys.frames()
  if (length(frames) && !is.null(frames[[1]]$ofile)) {
    return(dirname(normalizePath(frames[[1]]$ofile, mustWork = TRUE)))
  }

  normalizePath(getwd(), mustWork = TRUE)
}

SCRIPT_DIR <- get_script_dir()
SOURCE_DIR <- normalizePath(
  Sys.getenv("FIG1_SOURCE_DIR", file.path(SCRIPT_DIR, "..", "SourceData")),
  mustWork = TRUE
)
OUTPUT_DIR <- Sys.getenv("FIG1_OUTPUT_DIR", file.path(SCRIPT_DIR, "..", "Figure1_panels"))
EXPORT_PANELS <- tolower(Sys.getenv("FIG1_EXPORT", "true")) %in% c("1", "true", "yes")

cohort_path <- file.path(SOURCE_DIR, "KB_Cohort_SourceData.xlsx")
kbpap_path  <- file.path(SOURCE_DIR, "kbpap_SourceData.xlsx")
smg_path    <- file.path(SOURCE_DIR, "smg_SourceData.xlsx")

input_paths <- c(cohort_path, kbpap_path, smg_path)
if (any(!file.exists(input_paths))) {
  stop("Missing Figure 1 Source Data file(s): ",
       paste(basename(input_paths[!file.exists(input_paths)]), collapse = ", "))
}

INOCULUM    <- 1e6
H_CFU_FLOOR <- 800
M_CFU       <- 2
LOD_CFU     <- 0.5
PRIMARY_TIME <- "72"

CONC       <- c(128, 64, 32, 16, 8, 4, 2, 1, 0.5)
MIC_TP     <- "48h"
SMG_TP     <- "72h"
POS_MIN    <- 0.15

G_HR   <- "HR-candidate"
G_CTRL <- "Control"
G_R    <- "Resistant"
G_H99  <- "H99"
GROUP_LEVELS <- c(G_HR, G_CTRL, G_R, G_H99)
group_colors <- c(
  "HR-candidate" = "#E69F00",
  "Control"      = "#009E73",
  "Resistant"    = "#D55E00",
  "H99"          = "#000000"
)

safe_max <- function(x) if (length(x)) max(x) else NA_real_
safe_med <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) stats::median(x) else NA_real_
}

## -----------------------------------------------------------------------------
## Panel B: cohort FLC disk-diffusion distributions
## -----------------------------------------------------------------------------
cohort <- read_excel(cohort_path, sheet = "Fig1B")

required_cohort <- c(
  "Isolate_ID", "FLC25_inner_mm", "FLC25_outer_mm", "FLC25_delta_mm"
)
stopifnot(all(required_cohort %in% names(cohort)))
stopifnot(nrow(cohort) == 144L)
stopifnot(!anyDuplicated(cohort$Isolate_ID))
stopifnot(all(complete.cases(cohort[required_cohort])))
stopifnot(all(abs(
  cohort$FLC25_outer_mm - cohort$FLC25_inner_mm - cohort$FLC25_delta_mm
) < 1e-8))

kb_data <- cohort %>%
  select(FLC25_inner_mm, FLC25_outer_mm) %>%
  pivot_longer(everything(), names_to = "Measurement", values_to = "Value") %>%
  mutate(
    Measurement = recode(
      Measurement,
      FLC25_inner_mm = "FLC25 (Inner inhibition zone)",
      FLC25_outer_mm = "FLC25 (Outer inhibition zone)"
    ),
    Measurement = factor(
      Measurement,
      levels = c(
        "FLC25 (Inner inhibition zone)",
        "FLC25 (Outer inhibition zone)"
      )
    )
  )

normality_test <- kb_data %>%
  group_by(Measurement) %>%
  summarise(
    Sample_Size = n(),
    Mean = mean(Value),
    SD = sd(Value),
    Shapiro_W = unname(shapiro.test(Value)$statistic),
    P_value = shapiro.test(Value)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    Label = paste0(
      "N=", Sample_Size, "\n",
      "Mean=", round(Mean, 2), "\n",
      "p=", format.pval(P_value, digits = 3)
    )
  )

normal_curve <- kb_data %>%
  group_by(Measurement) %>%
  group_modify(~ {
    x <- seq(min(.x$Value), max(.x$Value), length.out = 200)
    tibble(Value = x, density = dnorm(x, mean(.x$Value), sd(.x$Value)))
  }) %>%
  ungroup()

delta_data <- cohort %>%
  transmute(
    Drug = factor("Fluconazole (FLC)", levels = "Fluconazole (FLC)"),
    Delta_Value = FLC25_delta_mm
  )

delta_stats_summary <- delta_data %>%
  summarise(
    Drug = first(Drug),
    N = n(),
    Mean = mean(Delta_Value),
    Max = max(Delta_Value),
    Label = paste0("N = ", N, "\nMean = ", round(Mean, 1), " mm")
  )

pB_hist <- ggplot(kb_data, aes(Value)) +
  geom_histogram(
    aes(y = after_stat(density)), bins = 20,
    fill = "#a8dadc", color = "white", alpha = 0.7
  ) +
  geom_density(color = "#e63946", linewidth = 1, fill = "#e63946", alpha = 0.1) +
  geom_line(
    data = normal_curve, aes(y = density),
    color = "gray40", linetype = "dashed", linewidth = 0.5
  ) +
  facet_wrap(~Measurement, scales = "free", ncol = 1) +
  geom_text(
    data = normality_test, aes(x = Inf, y = Inf, label = Label),
    hjust = 1.1, vjust = 1.1, size = 3.5, color = "gray20"
  ) +
  labs(x = "Diameter (mm)", y = "Density") +
  theme_gray(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold", color = "#1d3557"),
    panel.grid.minor = element_blank()
  )

pB_delta <- ggplot(delta_data, aes(Drug, Delta_Value, fill = Drug)) +
  geom_violin(alpha = 0.4, color = NA, trim = FALSE, width = 0.7) +
  geom_boxplot(width = 0.1, color = "gray30", alpha = 0.8, outlier.shape = NA) +
  geom_jitter(aes(color = Drug), width = 0.15, size = 1.5, alpha = 0.6) +
  geom_text(
    data = delta_stats_summary,
    aes(y = Max + 2, label = Label),
    color = "black", size = 3.5, fontface = "italic"
  ) +
  scale_fill_manual(values = c("Fluconazole (FLC)" = "#E15759")) +
  scale_color_manual(values = c("Fluconazole (FLC)" = "#9e3d3e")) +
  labs(x = NULL, y = "Delta zone diameter (outer - inner) [mm]") +
  theme_gray(base_size = 10) +
  theme(
    axis.text.x = element_text(face = "bold", color = "black"),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

panelB <- pB_hist + pB_delta + plot_layout(widths = c(2, 1))

## -----------------------------------------------------------------------------
## Shared data for Panels C-F and H
## -----------------------------------------------------------------------------
clin <- read_excel(kbpap_path, sheet = "Clinical_Isolates")

delta_cols <- c(
  "ori_delta", "re_1_delta", "re_2_delta",
  "re_3_delta", "re_4_delta", "re_5_delta"
)
rep_delta_cols <- delta_cols[-1]
pap_cols <- grep("^FLC_", names(clin), value = TRUE)

required_kbpap <- c(
  "ID", "Group", "Broth_FLC_ori", "Broth_add",
  delta_cols, pap_cols
)
stopifnot(all(required_kbpap %in% names(clin)))
stopifnot(nrow(clin) == 20L)
stopifnot(!anyDuplicated(as.character(clin$ID)))

map_group <- function(x) {
  case_when(
    str_to_lower(x) == "hr" ~ G_HR,
    str_to_lower(x) == "control" ~ G_CTRL,
    str_to_lower(x) == "r" ~ G_R,
    str_to_lower(x) %in% c("standard", "h99") ~ G_H99,
    TRUE ~ as.character(x)
  )
}

meta <- clin %>%
  transmute(
    ID = as.character(ID),
    Group = factor(map_group(Group), levels = GROUP_LEVELS),
    MIC = coalesce(Broth_FLC_ori, Broth_add)
  )

stopifnot(!anyNA(meta$MIC))
stopifnot(identical(
  as.integer(table(meta$Group)),
  c(10L, 6L, 3L, 1L)
))

grp_levels <- c(G_R, G_CTRL, G_HR)
grp_cols <- group_colors[grp_levels]

## -----------------------------------------------------------------------------
## Panel C: phenotypic grouping rationale
## -----------------------------------------------------------------------------
delta_long <- clin %>%
  filter(Group %in% c("HR", "control", "R")) %>%
  transmute(ID = as.character(ID), Group = map_group(Group), across(all_of(delta_cols))) %>%
  pivot_longer(all_of(delta_cols), names_to = "rep", values_to = "delta") %>%
  filter(!is.na(delta)) %>%
  mutate(Group = factor(Group, levels = grp_levels), x = 1)

delta_consensus <- delta_long %>%
  group_by(ID, Group) %>%
  summarise(delta_med = median(delta), .groups = "drop") %>%
  mutate(x = 1)

stopifnot(nrow(delta_long) == 114L, nrow(delta_consensus) == 19L)

grp_bands <- delta_consensus %>%
  group_by(Group) %>%
  summarise(
    ymin = min(delta_med), ymax = max(delta_med), n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    ymin_band = if_else(Group == G_R, -0.6, ymin - 0.3),
    ymax_band = if_else(Group == G_R, 0.9, ymax + 0.3),
    ymid = (ymin_band + ymax_band) / 2,
    label = if_else(
      Group == G_R,
      paste0(Group, "\nDelta = 0 (no zone)   n=", n),
      paste0(Group, "\nDelta ", round(ymin, 1), "-", round(ymax, 1), " mm   n=", n)
    )
  )

set.seed(1)
panelC <- ggplot() +
  geom_rect(
    data = grp_bands,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin_band, ymax = ymax_band, fill = Group),
    alpha = 0.10
  ) +
  geom_violin(
    data = delta_long, aes(x, delta, group = 1),
    width = 0.85, fill = "grey88", color = "grey70", alpha = 0.5, trim = FALSE
  ) +
  geom_jitter(
    data = delta_long, aes(x, delta, color = Group),
    width = 0.16, height = 0, size = 1.4, alpha = 0.35
  ) +
  geom_point(
    data = delta_consensus, aes(x, delta_med, color = Group),
    position = position_jitter(width = 0.05, height = 0, seed = 1),
    shape = 23, size = 3, stroke = 0.6, fill = "white"
  ) +
  geom_point(
    data = delta_consensus, aes(x, delta_med, color = Group),
    position = position_jitter(width = 0.05, height = 0, seed = 1),
    shape = 18, size = 2.6
  ) +
  geom_text(
    data = grp_bands, aes(x = 1.62, y = ymid, label = label, color = Group),
    hjust = 0, size = 4, fontface = "bold", lineheight = 0.95
  ) +
  scale_fill_manual(values = grp_cols) +
  scale_color_manual(values = grp_cols) +
  scale_x_continuous(limits = c(0.4, 2.9), breaks = NULL, name = NULL) +
  scale_y_continuous(breaks = seq(0, 25, 5)) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Phenotypic grouping rationale",
    y = "Delta zone diameter (outer - inner) [mm]"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

## -----------------------------------------------------------------------------
## Panel D: original versus repeated disk-diffusion measurements
## -----------------------------------------------------------------------------
e_ori <- clin %>%
  filter(Group %in% c("HR", "control", "R")) %>%
  transmute(
    ID = as.character(ID),
    Group = factor(map_group(Group), levels = grp_levels),
    delta = ori_delta
  ) %>%
  filter(!is.na(delta))

e_rep <- clin %>%
  filter(Group %in% c("HR", "control", "R")) %>%
  transmute(ID = as.character(ID), Group = map_group(Group), across(all_of(rep_delta_cols))) %>%
  pivot_longer(all_of(rep_delta_cols), names_to = "rep", values_to = "delta") %>%
  filter(!is.na(delta)) %>%
  mutate(Group = factor(Group, levels = grp_levels))

stopifnot(nrow(e_ori) == 19L, nrow(e_rep) == 95L)

comparisons <- list(
  c(G_R, G_CTRL),
  c(G_CTRL, G_HR),
  c(G_R, G_HR)
)

mk_delta_box <- function(d, subtitle) {
  ytop <- max(d$delta)
  ystep <- diff(range(d$delta))
  if (ystep == 0) ystep <- 1

  ggplot(d, aes(Group, delta, color = Group)) +
    geom_boxplot(width = 0.55, fill = NA, linewidth = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.16, height = 0, size = 1.6, alpha = 0.55) +
    stat_compare_means(
      comparisons = comparisons, method = "wilcox.test",
      size = 3, tip.length = 0.01,
      label.y = ytop + ystep * c(0.12, 0.28, 0.44)
    ) +
    stat_compare_means(
      method = "kruskal.test", size = 3, color = "gray45",
      label.x = 0.6, label.y = ytop + ystep * 0.62
    ) +
    scale_color_manual(values = grp_cols, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    coord_cartesian(clip = "off") +
    labs(
      subtitle = subtitle, x = NULL,
      y = "Delta zone diameter (outer - inner) [mm]"
    ) +
    theme_pubr(base_size = 11) +
    theme(
      legend.position = "none",
      plot.subtitle = element_text(face = "bold", color = "#1d3557"),
      plot.margin = margin(t = 16, r = 10, b = 6, l = 8)
    )
}

panelD <- mk_delta_box(e_ori, "Original (ori)") |
  mk_delta_box(e_rep, "Replicates (re 1-5)")

## -----------------------------------------------------------------------------
## Panels E-F: 72 h CFU PAP and discrete-colony reach
## -----------------------------------------------------------------------------
cfu_long <- clin %>%
  mutate(
    ID = as.character(ID),
    across(all_of(pap_cols), as.character)
  ) %>%
  select(ID, all_of(pap_cols)) %>%
  pivot_longer(all_of(pap_cols), names_to = "col", values_to = "raw") %>%
  mutate(
    plate = if_else(str_detect(col, "^FLC_M_"), "M", "H"),
    body = str_remove(col, "^FLC_(M_)?"),
    conc = as.numeric(str_extract(body, "^[0-9.]+")),
    time = str_extract(body, "[0-9]+$"),
    raw_num = suppressWarnings(as.numeric(raw)),
    cfu_class = case_when(
      raw == "H" ~ "H",
      raw == "M" ~ "M",
      !is.na(raw_num) & raw_num == 0 ~ "0",
      !is.na(raw_num) & raw_num > 0 ~ "N",
      TRUE ~ NA_character_
    ),
    cfu = case_when(
      cfu_class == "H" ~ H_CFU_FLOOR,
      cfu_class == "M" ~ M_CFU,
      cfu_class == "N" ~ raw_num,
      cfu_class == "0" ~ 0,
      TRUE ~ NA_real_
    ),
    survival = cfu / INOCULUM
  ) %>%
  left_join(meta, by = "ID") %>%
  mutate(fold_mic = conc / MIC)

curve_concs <- c(0, 0.5, 1, 2, 4, 8, 16, 32, 64, 128)
cfu_curve <- cfu_long %>%
  filter(
    time == PRIMARY_TIME,
    plate == "H",
    conc %in% curve_concs,
    !is.na(fold_mic), fold_mic > 0
  ) %>%
  mutate(
    survival_plot = if_else(cfu == 0, LOD_CFU / INOCULUM, survival),
    censored = cfu_class == "H"
  )

stopifnot(nrow(cfu_curve) == 180L)

main_folds <- c(1, 2, 4, 8, 16)
cfu_main_band <- cfu_curve %>%
  filter(Group %in% c(G_HR, G_CTRL), fold_mic %in% main_folds) %>%
  group_by(Group, fold_mic) %>%
  summarise(
    mean = mean(survival_plot),
    lo = quantile(survival_plot, 0.25),
    hi = quantile(survival_plot, 0.75),
    .groups = "drop"
  )

cfu_main_ref <- cfu_curve %>%
  filter(Group %in% c(G_R, G_H99), fold_mic %in% main_folds) %>%
  group_by(Group, fold_mic) %>%
  summarise(mean = mean(survival_plot), .groups = "drop")

full_log2_breaks <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128)
full_log2_labels <- ifelse(
  full_log2_breaks < 1,
  sub("0+$", "", sprintf("%.2f", full_log2_breaks)),
  as.character(full_log2_breaks)
)

panelE <- ggplot() +
  annotate(
    "rect", xmin = 4, xmax = 16, ymin = LOD_CFU / INOCULUM, ymax = Inf,
    fill = "grey95", alpha = 0.6
  ) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "grey55") +
  geom_ribbon(
    data = cfu_main_band,
    aes(fold_mic, ymin = lo, ymax = hi, fill = Group),
    alpha = 0.18
  ) +
  geom_line(
    data = cfu_main_band,
    aes(fold_mic, mean, color = Group),
    linewidth = 1.2
  ) +
  geom_point(
    data = cfu_main_band,
    aes(fold_mic, mean, color = Group),
    size = 2.4
  ) +
  geom_line(
    data = cfu_main_ref,
    aes(fold_mic, mean, color = Group),
    linetype = "dashed", linewidth = 0.8, alpha = 0.65
  ) +
  geom_point(
    data = cfu_main_ref,
    aes(fold_mic, mean, color = Group),
    shape = 1, size = 2, alpha = 0.65
  ) +
  scale_x_continuous(
    trans = "log2", breaks = full_log2_breaks, labels = full_log2_labels
  ) +
  scale_y_log10() +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors, guide = "none") +
  labs(
    x = "FLC / MIC (fold)",
    y = "Survival fraction (CFU / 1e6, log)",
    color = "Group"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

cfu_lhf_tbl <- cfu_long %>%
  filter(time == PRIMARY_TIME, !is.na(fold_mic), fold_mic > 0) %>%
  group_by(ID, Group, MIC) %>%
  summarise(
    lhf_N = safe_max(fold_mic[cfu_class == "N"]),
    .groups = "drop"
  )

lhf_plot_data <- cfu_lhf_tbl %>%
  filter(Group %in% c(G_HR, G_CTRL), !is.na(lhf_N))

auc_data <- lhf_plot_data %>%
  mutate(case = Group == G_HR)
auc_lhf <- with(
  auc_data,
  mean(outer(lhf_N[case], lhf_N[!case], `>`)) +
    0.5 * mean(outer(lhf_N[case], lhf_N[!case], `==`))
)

stopifnot(nrow(lhf_plot_data) == 15L)
stopifnot(abs(auc_lhf - 0.7222222) < 1e-6)
stopifnot(is.na(cfu_lhf_tbl$lhf_N[cfu_lhf_tbl$ID == "213"]))

panelF <- ggplot(lhf_plot_data, aes(Group, lhf_N, color = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.85) +
  scale_y_continuous(
    trans = "log2", breaks = c(1, 2, 4, 8, 16, 32),
    labels = paste0(c(1, 2, 4, 8, 16, 32), "x MIC")
  ) +
  scale_color_manual(values = group_colors) +
  labs(
    x = NULL,
    y = expression(LHF[N]~"(discrete-colony reach above MIC)"),
    title = sprintf("HR-candidate vs. Control AUROC ~ %.2f", auc_lhf)
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  )

## -----------------------------------------------------------------------------
## Panels G-H: SMG calculated directly from raw optical-density readings
## -----------------------------------------------------------------------------
read_block <- function(raw, timepoint, block_index, batch) {
  label_col <- 3 + (block_index - 1) * 13
  conc_cols <- (label_col + 1):(label_col + 9)
  pos_col <- label_col + 10
  neg_cols <- c(label_col + 11, label_col + 12)

  if (max(neg_cols) > ncol(raw)) {
    stop("Unexpected SMG sheet layout in ", batch, ".")
  }

  blank <- rowMeans(raw[, neg_cols, drop = FALSE], na.rm = TRUE)
  tibble(
    ID = rep(as.character(raw[[1]]), length(CONC)),
    rep = rep(raw[[2]], length(CONC)),
    batch = batch,
    timepoint = timepoint,
    conc = rep(CONC, each = nrow(raw)),
    od = as.vector(as.matrix(raw[, conc_cols, drop = FALSE])) - rep(blank, length(CONC)),
    pos = rep(raw[[pos_col]] - blank, length(CONC))
  )
}

read_smg_sheet <- function(path, sheet, timepoints, batch) {
  raw <- as.data.frame(read_excel(path, sheet = sheet, .name_repair = "minimal"))
  map_dfr(
    seq_along(timepoints),
    ~ read_block(raw, timepoints[[.x]], .x, batch)
  )
}

smg_long <- bind_rows(
  read_smg_sheet(smg_path, "ori_1", c("24h", "48h", "72h"), "batch1"),
  read_smg_sheet(smg_path, "ori_2", c("24h", "48h", "72h", "96h"), "batch2")
) %>%
  mutate(ID = if_else(str_detect(ID, "^H99"), G_H99, ID)) %>%
  left_join(meta %>% select(ID, Group), by = "ID") %>%
  mutate(
    Group = if_else(ID == G_H99, G_H99, as.character(Group)),
    Group = factor(Group, levels = GROUP_LEVELS)
  )

stopifnot(setequal(unique(smg_long$ID), meta$ID))

mic50_one <- function(conc, od, pos) {
  if (!is.finite(pos) || pos < POS_MIN) return(NA_real_)
  ord <- order(conc)
  conc <- conc[ord]
  od <- od[ord]
  inhibited <- which(od <= 0.5 * pos)
  if (length(inhibited)) conc[min(inhibited)] else max(CONC) * 2
}

mic_rep_tbl <- smg_long %>%
  filter(timepoint %in% c(MIC_TP, SMG_TP)) %>%
  group_by(ID, batch, rep) %>%
  summarise(
    mic_48 = mic50_one(
      conc[timepoint == MIC_TP],
      od[timepoint == MIC_TP],
      first(pos[timepoint == MIC_TP])
    ),
    mic_72 = mic50_one(
      conc[timepoint == SMG_TP],
      od[timepoint == SMG_TP],
      first(pos[timepoint == SMG_TP])
    ),
    .groups = "drop"
  ) %>%
  mutate(MIC_rep = if_else(is.na(mic_48), mic_72, mic_48))

mic_plate <- mic_rep_tbl %>%
  group_by(ID) %>%
  summarise(MIC_use = safe_med(MIC_rep), .groups = "drop")

smg_rep <- smg_long %>%
  filter(timepoint == SMG_TP) %>%
  left_join(mic_plate, by = "ID") %>%
  group_by(ID, batch, rep, MIC_use) %>%
  summarise(
    SMG = {
      above_mic <- od[conc > MIC_use]
      if (length(above_mic)) mean(above_mic) / first(pos) else NA_real_
    },
    .groups = "drop"
  )

smg_tbl <- smg_rep %>%
  group_by(ID) %>%
  summarise(
    smg = mean(SMG, na.rm = TRUE),
    smg_sd = sd(SMG, na.rm = TRUE),
    smg_n = sum(!is.na(SMG)),
    .groups = "drop"
  ) %>%
  left_join(meta %>% select(ID, Group), by = "ID") %>%
  arrange(Group, desc(smg))

stopifnot(nrow(smg_tbl) == 20L)
stopifnot(all(smg_tbl$smg_n[smg_tbl$ID != G_H99] == 4L))
stopifnot(smg_tbl$smg_n[smg_tbl$ID == G_H99] == 8L)

panelG_data <- smg_tbl %>% filter(!is.na(Group), !is.na(smg))
panelG_p <- wilcox.test(
  smg ~ Group,
  data = panelG_data %>% filter(Group %in% c(G_HR, G_CTRL))
)$p.value

panelG <- ggplot(panelG_data, aes(Group, smg, color = Group)) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.4) +
  geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.45) +
  geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.85) +
  geom_text_repel(aes(label = ID), size = 2.4, max.overlaps = Inf, show.legend = FALSE) +
  stat_compare_means(
    comparisons = list(c(G_HR, G_CTRL)), method = "wilcox.test",
    label = "p.format",
    label.y = max(panelG_data$smg) * 1.12,
    tip.length = 0.01, bracket.size = 0.35, size = 3.2
  ) +
  scale_color_manual(values = group_colors) +
  coord_cartesian(
    ylim = c(min(panelG_data$smg) * 1.1, max(panelG_data$smg) * 1.22),
    clip = "off"
  ) +
  labs(x = NULL, y = expression(SMG[72*h]), title = "FLC tolerance (SMG) by group") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    plot.margin = margin(5.5, 10, 5.5, 5.5)
  )

axes_data <- cfu_lhf_tbl %>%
  select(ID, Group, lhf_N) %>%
  left_join(smg_tbl %>% select(ID, smg), by = "ID") %>%
  filter(Group %in% c(G_HR, G_CTRL), !is.na(lhf_N), !is.na(smg))

axes_rho <- cor(axes_data$lhf_N, axes_data$smg, method = "spearman")
axes_med_x <- median(axes_data$lhf_N)
axes_med_y <- median(axes_data$smg)
axes_arch <- axes_data %>%
  filter(ID %in% c("215", "160")) %>%
  mutate(
    role = if_else(
      ID == "215", "reach-type (Cn215)", "tolerance-type (Cn160)"
    )
  )

stopifnot(nrow(axes_data) == 15L)
stopifnot(abs(axes_rho - 0.174378) < 1e-6)

panelH <- ggplot(axes_data, aes(lhf_N, smg)) +
  geom_hline(yintercept = axes_med_y, linetype = "dotted", color = "grey80") +
  geom_vline(xintercept = axes_med_x, linetype = "dotted", color = "grey80") +
  geom_point(aes(color = Group), size = 3, alpha = 0.9) +
  geom_point(
    data = axes_arch,
    shape = 21, size = 5.5, stroke = 1.2, color = "black", fill = NA
  ) +
  geom_text_repel(
    data = axes_arch, aes(label = role),
    size = 3, fontface = "bold", min.segment.length = 0,
    box.padding = 1, show.legend = FALSE
  ) +
  scale_x_continuous(trans = "log2", breaks = c(4, 8, 16)) +
  scale_color_manual(values = group_colors) +
  labs(
    x = "PAP discrete-subpopulation reach (lhf_N, fold-MIC)",
    y = expression(SMG[72*h]~"(tolerance)"),
    title = "Two orthogonal phenotype axes",
    subtitle = sprintf("Spearman rho = %.2f (n = %d, HR+Control)", axes_rho, nrow(axes_data))
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))


## -----------------------------------------------------------------------------
## PDF output
## -----------------------------------------------------------------------------
if (EXPORT_PANELS) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(OUTPUT_DIR)) stop("Cannot create output directory: ", OUTPUT_DIR)
  if (!capabilities("cairo")) stop("PDF output requires Cairo support in R.")
  panel_plots <- list(
    Fig1B = panelB, Fig1C = panelC, Fig1D = panelD, Fig1E = panelE,
    Fig1F = panelF, Fig1G = panelG, Fig1H = panelH
  )
  panel_sizes <- list(
    Fig1B = c(8.0, 6.6), Fig1C = c(5.8, 5.8), Fig1D = c(6.8, 4.0),
    Fig1E = c(6.7, 5.8), Fig1F = c(5.2, 4.5), Fig1G = c(4.3, 3.6),
    Fig1H = c(5.5, 3.9)
  )
  iwalk(panel_plots, function(p, name) {
    size <- panel_sizes[[name]]
    ggsave(
      file.path(OUTPUT_DIR, paste0(name, ".pdf")), p,
      width = size[1], height = size[2], units = "in", device = cairo_pdf
    )
  })
  message("Figure 1 PDFs written to: ", normalizePath(OUTPUT_DIR))
}

