# Figure 6 plot code
#
# Inputs: ../SourceData. Default PDF output: ../Figure6_panels.
# Set FIG6_EXPORT=false to calculate without writing PDFs.
# FIG6_OUTPUT_DIR can redirect PDF output; no input files are modified.

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  # ape and ggtree mask a few tidyverse verbs, so the tidy core is attached last
  # and wins on the search path.
  for (p in c("purrr", "tibble", "readr", "tidyr", "dplyr", "stringr")) {
    try(detach(paste0("package:", p), character.only = TRUE, unload = FALSE), silent = TRUE)
  }
  library(purrr); library(tibble); library(readr); library(tidyr)
  library(dplyr); library(stringr)
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
  dir_candidates <- file.path(getwd(), c(".", "Codes", "Revise02/Codes"))
  found <- dir_candidates[file.exists(file.path(dir_candidates, "fig6_submit.R"))]
  if (length(found) != 1L) {
    stop("Run fig6_submit.R with source() or Rscript so its directory can be resolved.")
  }
  normalizePath(found, mustWork = TRUE)
}
SCRIPT_DIR <- get_script_dir()
SOURCE_DIR <- normalizePath(file.path(SCRIPT_DIR, "..", "SourceData"), mustWork = TRUE)
OUTPUT_DIR <- Sys.getenv("FIG6_OUTPUT_DIR", file.path(SCRIPT_DIR, "..", "Figure6_panels"))
EXPORT_PANELS <- tolower(Sys.getenv("FIG6_EXPORT", "true")) %in% c("1", "true", "yes")

F_MASTER <- file.path(SOURCE_DIR, "master_wgs_table.tsv")
F_TRACKA <- file.path(SOURCE_DIR, "trackA_ploidy_summary.tsv")
F_DEPTH  <- file.path(SOURCE_DIR, "mosdepth_chr_normalized_depth_long.tsv")
F_HETCHR <- file.path(SOURCE_DIR, "het_chrom_summary.tsv")
F_KSNP   <- file.path(SOURCE_DIR, "tree.core_SNPs.ML.tre")
F_MLST   <- file.path(SOURCE_DIR, "Cn_all_merged_results.txt")
input_files <- c(F_MASTER, F_TRACKA, F_DEPTH, F_HETCHR, F_KSNP, F_MLST)
if (any(!file.exists(input_files))) {
  stop("Missing SourceData input(s): ",
       paste(basename(input_files[!file.exists(input_files)]), collapse = ", "))
}

require_columns <- function(d, columns, label) {
  missing <- setdiff(columns, names(d))
  if (length(missing)) stop(label, " is missing: ", paste(missing, collapse = ", "))
}


## -----------------------------------------------------------------------------
## Conventions
## -----------------------------------------------------------------------------
G_HR   <- "HR-candidate"
G_CTRL <- "Control"
G_R    <- "Resistant"
G_H99  <- "H99"
group_levels <- c(G_HR, G_CTRL, G_R, G_H99)
group_colors <- setNames(c("#E69F00", "#009E73", "#D55E00", "#000000"), group_levels)
GRP_COL <- group_colors
GRP_SHP <- setNames(c(16, 17, 18, 15), group_levels)

# Integrated WGS genome-state classes.
CLASS_COL <- c(monoclonal_clean = "#009E73", monoclonal_with_disomy = "#E69F00",
               heterozygous_diploid = "#0072B2", divergent_clone_artifact = "#D55E00",
               divergent_clone_artifact_marginal = "#CC79A7",
               ambiguous_needs_review = "#999999", reference_H99 = "#000000")
CLASS_LAB <- c(monoclonal_clean = "clean haploid", monoclonal_with_disomy = "haploid + disomy",
               heterozygous_diploid = "heterozygous diploid",
               divergent_clone_artifact = "divergent (reference artifact)",
               divergent_clone_artifact_marginal = "divergent (marginal)",
               ambiguous_needs_review = "ambiguous", reference_H99 = "H99 reference")

CHR_ACC <- sprintf("NC_0267%02d.1", 45:58)
CHR_NUM <- setNames(seq_along(CHR_ACC), CHR_ACC)
chr_levels <- paste0("Chr", CHR_NUM[CHR_ACC])
chr_factor <- function(acc) factor(paste0("Chr", CHR_NUM[acc]), levels = chr_levels)

# Thresholds drawn as guides in Panel A
GAIN_THRESHOLD <- 1.50   # whole-chromosome gain, cohort-relative depth ratio
HETZ_THRESHOLD <- 3      # reference-het robust z flagged as an outlier

parse_pct <- function(x) as.numeric(str_remove(x, "%"))
# CnXXXBB -> parent CnXXX; H99 clones map back to CnH99.
parent_of <- function(s) ifelse(str_detect(s, "^CnH99"), "CnH99",
                                paste0("Cn", str_sub(str_remove(s, "^Cn"), 1, -3)))

# The 19 clinical parents in phenotype-group order.
COHORT <- tribble(
  ~strain, ~group_code,
  "101", "HR",  "106", "HR",  "126", "HR",  "128", "HR",  "160", "HR",
  "162", "HR",  "208", "HR",  "213", "HR",  "215", "HR",  "217", "HR",
  "147", "control", "150", "control", "169", "control",
  "174", "control", "225", "control", "226", "control",
  "196", "R", "197", "R", "198", "R"
)
CODE_TO_GROUP <- c(HR = G_HR, control = G_CTRL, R = G_R)
grp_of <- setNames(unname(CODE_TO_GROUP[COHORT$group_code]), paste0("Cn", COHORT$strain))

# Tree display switches.
OUTGROUP      <- "CnH99"   # root both trees on H99 as the outgroup
BREAK_LONG    <- TRUE      # truncate an outlier long branch for display (self-disabling)
BREAK_IDS     <- "Cn160"   # tip(s) eligible for the display break
BREAK_TRIGGER <- 1.8       # break only if root-to-tip > this * next-furthest tip
BREAK_TARGET  <- 1.15      # place the broken tip at this * next-furthest root-to-tip
ST_NOVEL_FROM <- c("0")    # stringMLST ST codes rendered as "ST novel"
SUPPORT_CUTOFF <- 70       # only label nodes at or above this support
# Support labels sit just above the branch they belong to. If a label collides
# with its branch and reads as a partial number, raise SUPPORT_NUDGE_Y.
SUPPORT_NUDGE_Y <- 0.16


## -----------------------------------------------------------------------------
## 1. Load
## -----------------------------------------------------------------------------
master <- read_tsv(F_MASTER, show_col_types = FALSE)
require_columns(master, c("sample", "clone_type", "call_wgs_class", "screen_het_robust_z"),
                basename(F_MASTER))
trackA <- read_tsv(F_TRACKA, show_col_types = FALSE)
require_columns(trackA, c("sample", "gs_heterozygous_ab"), basename(F_TRACKA))
trackA <- trackA %>% mutate(kmer_het = parse_pct(gs_heterozygous_ab))
depth  <- read_tsv(F_DEPTH, show_col_types = FALSE)
require_columns(depth, c("sample", "chrom", "chr_depth_norm"), basename(F_DEPTH))
hetchr <- read_tsv(F_HETCHR, show_col_types = FALSE)
require_columns(hetchr, c("sample", "chrom", "het_per_kb_chrom"), basename(F_HETCHR))
kmer_of  <- setNames(trackA$kmer_het, trackA$sample)

# Cohort-relative chromosome ratio
cohort_med <- depth %>% filter(chrom %in% CHR_ACC) %>% group_by(chrom) %>%
  summarise(m = median(chr_depth_norm, na.rm = TRUE), .groups = "drop")
depth <- depth %>% left_join(cohort_med, by = "chrom") %>% mutate(rel = chr_depth_norm / m)


## -----------------------------------------------------------------------------
## 2. Panel A: cohort selection landscape
## -----------------------------------------------------------------------------
Co <- master %>%
  filter(clone_type == "parent", sample %in% names(grp_of)) %>%
  transmute(sample, group = factor(grp_of[sample], levels = group_levels),
            hetz = asinh(screen_het_robust_z)) %>%
  left_join(depth %>% filter(chrom %in% CHR_ACC) %>% group_by(sample) %>%
              summarise(gain = max(rel, na.rm = TRUE), .groups = "drop"), by = "sample")
if (nrow(Co) != length(grp_of)) {
  stop("Panel A expects one Parent row per clinical isolate; got ", nrow(Co),
       " of ", length(grp_of), ".")
}
if (any(!is.finite(Co$gain)) || any(!is.finite(Co$hetz))) {
  stop("Panel A has a Parent isolate with no depth or heterozygosity value.")
}

z_ticks <- c(-1, 0, 1, 5, 20, 68)
p_cohort <- ggplot(Co, aes(gain, hetz, color = group, shape = group)) +
  geom_vline(xintercept = GAIN_THRESHOLD, linetype = "dashed",
             colour = "grey80", linewidth = 0.3) +
  geom_hline(yintercept = asinh(HETZ_THRESHOLD), linetype = "dashed",
             colour = "grey80", linewidth = 0.3) +
  geom_point(size = 3) +
  geom_text_repel(
    data = filter(Co, sample %in% paste0("Cn", c("160", "150", "147", "128",
                                                 "215", "196", "197", "198"))),
    aes(label = sample), size = 2.8, colour = "grey25",
    min.segment.length = 0, max.overlaps = Inf, seed = 1) +
  annotate("text", x = 1.0, y = asinh(60), hjust = 0, size = 2.6, colour = "#D55E00",
           label = "'mixture' (reference)\n— k-mer: divergent artifact") +
  annotate("text", x = 1.9, y = asinh(-1.6), hjust = 1, size = 2.8, colour = "grey45",
           label = "clonal aneuploidy") +
  annotate("text", x = 1.0, y = asinh(-1.6), hjust = 0, size = 2.8, colour = "grey45",
           label = "clean monoclonal") +
  scale_color_manual(values = GRP_COL, name = NULL, drop = TRUE) +
  scale_shape_manual(values = GRP_SHP, name = NULL, drop = TRUE) +
  scale_y_continuous(breaks = asinh(z_ticks), labels = z_ticks) +
  labs(x = "Strongest chromosomal gain (cohort-relative depth ratio; Chr2 absorbed)",
       y = "Genotypic mixture (reference het robust z, asinh)") +
  theme_classic(base_size = 11) +
  theme(legend.position = c(0.87, 0.78))
panelA <- p_cohort


## -----------------------------------------------------------------------------
## 3. Panel C: k-mer heterozygosity change on cloning
## -----------------------------------------------------------------------------
# 01 is Selected, 02 Monoclonal and 03 Passage.
dA <- master %>% filter(clone_type == "clone") %>%
  transmute(sample, parent = parent_of(sample),
            cls = factor(recode(call_wgs_class,
                                divergent_clone_artifact_marginal = "divergent_clone_artifact"),
                         levels = names(CLASS_COL))) %>%
  mutate(kmer = kmer_of[sample], pkmer = kmer_of[parent], dkmer = kmer - pkmer) %>%
  filter(!is.na(kmer), !is.na(pkmer))
if (!nrow(dA)) stop("Panel C found no clone with both its own and its Parent k-mer value.")

pA <- ggplot(dA, aes(dkmer, kmer, color = cls)) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.4) +
  geom_point(size = 2.6, alpha = 0.9) +
  geom_text_repel(
    data = filter(dA, sample %in% c("Cn12802", "Cn12803", "Cn14701",
                                    "Cn14702", "Cn14703", "Cn21501")),
    aes(label = sample), size = 2.4, colour = "grey25",
    min.segment.length = 0, max.overlaps = Inf, seed = 1) +
  annotate("text", x = max(dA$dkmer) * 0.5, y = max(dA$kmer) * 0.95, hjust = 0,
           size = 2.7, colour = "#0072B2",
           label = "acquired het (haploid → het-diploid)\ne.g. Cn128 series") +
  annotate("text", x = min(dA$dkmer) * 0.9, y = 3.5, hjust = 0,
           size = 2.7, colour = "#0072B2",
           label = "lost het on passage\nCn14703") +
  annotate("text", x = 0, y = -0.4, hjust = 0.5, size = 2.5, colour = "grey55",
           label = "divergent / clean cluster at origin → self-excluded") +
  scale_color_manual(values = CLASS_COL, labels = CLASS_LAB,
                     name = "WGS genome state", drop = TRUE) +
  labs(x = "Δ k-mer heterozygosity  (clone − parent, %)",
       y = "k-mer heterozygosity (current state, %)") +
  theme_classic(base_size = 11) +
  theme(legend.key.size = unit(3.2, "mm"), legend.text = element_text(size = 7),
        legend.title = element_text(size = 8))
panelC <- pA


## -----------------------------------------------------------------------------
## 4. Panel D: genome-state archetype karyograms
## -----------------------------------------------------------------------------
karyogram <- function(samps, labs) {
  missing <- setdiff(samps, depth$sample)
  if (length(missing)) stop("Karyogram sample(s) absent from the depth table: ",
                            paste(missing, collapse = ", "))
  lev <- paste0(samps, "  ·  ", labs); labmap <- setNames(lev, samps)
  dep <- depth %>% filter(sample %in% samps, chrom %in% CHR_ACC) %>%
    transmute(chrom = chr_factor(chrom), value = chr_depth_norm,
              facet = factor(labmap[sample], levels = lev))
  het <- hetchr %>% filter(sample %in% samps, chrom %in% CHR_ACC) %>%
    transmute(chrom = chr_factor(chrom), value = het_per_kb_chrom,
              facet = factor(labmap[sample], levels = lev))
  if (!nrow(het)) stop("Karyogram found no per-chromosome heterozygosity rows.")
  ann <- tibble(sample = samps) %>%
    mutate(kmer = kmer_of[sample], facet = factor(labmap[sample], levels = lev),
           kmer_txt = sprintf("k-mer het = %.1f%%", kmer))
  htop <- max(het$value, na.rm = TRUE) * 1.1
  pdep <- ggplot(dep, aes(chrom, value)) +
    geom_col(data = filter(dep, chrom != "Chr2"), width = 0.72,
             fill = "#4C78A8", alpha = 0.9) +
    geom_col(data = filter(dep, chrom == "Chr2"), width = 0.72,
             fill = "#4C78A8", alpha = 0.30) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
    geom_hline(yintercept = 2, linetype = "dotted", colour = "grey80", linewidth = 0.3) +
    geom_text(data = ann, aes(x = 1, y = 2.35, label = kmer_txt), inherit.aes = FALSE,
              hjust = 0, size = 2.4, fontface = "italic", colour = "grey30") +
    facet_wrap(~ facet, ncol = 1, strip.position = "top") +
    coord_cartesian(ylim = c(0, 2.6)) +
    scale_x_discrete(limits = chr_levels, drop = FALSE) +
    labs(x = NULL, y = "depth ratio") +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 7.5, hjust = 0))
  phet <- ggplot(het, aes(chrom, value)) +
    geom_col(data = filter(het, chrom != "Chr2"), width = 0.72,
             fill = "#B2182B", alpha = 0.85) +
    geom_col(data = filter(het, chrom == "Chr2"), width = 0.72,
             fill = "#B2182B", alpha = 0.30) +
    facet_wrap(~ facet, ncol = 1, strip.position = "top") +
    coord_cartesian(ylim = c(0, htop)) +
    scale_x_discrete(limits = chr_levels, drop = FALSE) +
    labs(x = NULL, y = "reference het / kb",
         title = "reference het (Chr2 faded = artifact)") +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
          strip.background = element_blank(),
          strip.text = element_text(size = 7.5, hjust = 0, colour = "grey40"),
          plot.title = element_text(size = 8.5))
  pdep | phet
}
pB_arch <- karyogram(
  c("Cn101", "Cn10101", "Cn147", "Cn160"),
  c("clean haploid", "aneuploidy · Chr1 disomy", "heterozygous diploid",
    "divergent (artifact)"))
panelD <- pB_arch


## -----------------------------------------------------------------------------
## 5. Panel E: reference-free core-SNP tree
## -----------------------------------------------------------------------------
load_class <- function(path) {
  m <- read_tsv(path, show_col_types = FALSE)
  sc <- names(m)[which(tolower(names(m)) %in% c("sample", "sample_id"))[1]]
  cc <- names(m)[which(tolower(names(m)) %in%
                         c("call_wgs_class", "wgs_class", "class"))[1]]
  if (is.na(sc) || is.na(cc)) {
    stop("Could not find sample and genome-state class columns in ", basename(path), ".")
  }
  setNames(as.character(m[[cc]]), str_extract(m[[sc]], "Cn[0-9A-Za-z]+"))
}
load_st <- function(path) {
  m <- suppressWarnings(read_tsv(path, show_col_types = FALSE))
  sc <- names(m)[which(tolower(names(m)) %in% c("sample", "sample_id"))[1]]
  st <- names(m)[which(toupper(names(m)) == "ST")[1]]
  if (is.na(sc) || is.na(st)) {
    stop("Could not find Sample and ST columns in ", basename(path), ".")
  }
  v <- setNames(as.character(m[[st]]), str_extract(m[[sc]], "Cn[0-9A-Za-z]+"))
  v[!duplicated(names(v))]
}
tree_class_of <- load_class(F_MASTER)
st_of <- load_st(F_MLST)

build_labels <- function(ids) {
  vapply(ids, function(id) {
    cl <- if (id %in% names(tree_class_of)) tree_class_of[[id]] else NA_character_
    st <- if (id %in% names(st_of)) st_of[[id]] else NA_character_
    st_disp <- if (is.na(st)) NA_character_
    else if (st %in% ST_NOVEL_FROM) "ST novel"
    else paste0("ST", st)
    parts <- c(if (!is.na(cl)) cl, if (!is.na(st_disp)) st_disp)
    if (length(parts)) paste0(id, "  [", paste(parts, collapse = " · "), "]") else id
  }, character(1))
}

root_ladder <- function(tr) {
  if (!is.null(tr$edge.length)) tr$edge.length[tr$edge.length < 0] <- 0
  if (OUTGROUP %in% tr$tip.label) {
    tr <- root(tr, outgroup = OUTGROUP, resolve.root = TRUE)
  } else {
    stop("Outgroup ", OUTGROUP, " is not among the tree tips.")
  }
  ladderize(tr)
}

# One phylogram, with an optional self-disabling long-branch break. Cn160 is
# reference-divergent, so its terminal branch is truncated for display and marked
# with a slash
plot_tree_v2 <- function(tr, support = NULL, layout = "rectangular",
                         branch.length = "branch.length",
                         right_mar = 12, do_legend = TRUE,
                         tip_size = NULL, tip_point = NULL,
                         legend_title = NULL, legend_text = NULL,
                         legend_position = NULL,
                         support_cutoff = SUPPORT_CUTOFF,
                         show_break_value = FALSE) {
  Ntip <- length(tr$tip.label)
  ids  <- tr$tip.label
  grp <- ifelse(ids == OUTGROUP, G_H99, unname(grp_of[ids]))
  grp[is.na(grp)] <- "unknown"
  pal <- c(GRP_COL, unknown = "grey55")

  tip_df <- tibble(label = ids, id = ids, group = grp, disp = build_labels(ids))

  tr_disp <- tr
  brk <- list()
  do_break <- BREAK_LONG && layout %in% c("rectangular", "roundrect", "slanted")
  if (do_break && !is.null(tr_disp$edge.length)) {
    nd <- node.depth.edgelength(tr_disp)
    rd <- nd[seq_len(Ntip)]
    for (bid in intersect(BREAK_IDS, ids)) {
      bt <- match(bid, ids)
      rd_other <- max(rd[-bt], na.rm = TRUE)
      if (is.finite(rd_other) && rd[bt] > BREAK_TRIGGER * rd_other) {
        e  <- which(tr_disp$edge[, 2] == bt)
        pd <- nd[tr_disp$edge[e, 1]]
        pos_len <- tr_disp$edge.length[tr_disp$edge.length > 0]
        min_pos <- if (length(pos_len)) min(pos_len, na.rm = TRUE) else 1e-6
        newlen <- max(BREAK_TARGET * rd_other - pd, min_pos * 2)
        brk[[bid]] <- list(id = bid, edge = e, true = tr_disp$edge.length[e], tip = bt)
        tr_disp$edge.length[e] <- newlen
      }
    }
  }

  p0 <- ggtree(tr_disp, layout = layout, branch.length = branch.length, linewidth = 0.45)
  dd <- p0$data
  h <- max(dd$x, na.rm = TRUE)
  if (!is.finite(h) || h <= 0) h <- 1
  off <- 0.012 * h

  p <- p0 %<+% tip_df +
    geom_tiplab(aes(label = disp, color = group), align = TRUE, linetype = "dotted",
                linesize = 0.25, offset = off, size = tip_size, show.legend = FALSE) +
    geom_tippoint(aes(color = group), size = tip_point, show.legend = do_legend)

  if (!is.null(support) && any(nzchar(as.character(support)), na.rm = TRUE)) {
    sup <- suppressWarnings(as.numeric(support))
    if (any(!is.na(sup))) {
      sup <- if (max(sup, na.rm = TRUE) <= 1.5) sup * 100 else sup
      n <- min(length(sup), tr_disp$Nnode)
      if (n > 0) {
        sup_df <- tibble(node = (Ntip + 1):(Ntip + n), sup = sup[seq_len(n)]) %>%
          mutate(lab = ifelse(!is.na(sup) & sup >= support_cutoff,
                              as.character(round(sup)), NA_character_))
        sup_xy <- dd %>% left_join(sup_df, by = "node") %>% filter(!isTip, !is.na(lab))
        if (nrow(sup_xy)) {
          p <- p + geom_text(data = sup_xy, aes(x = x, y = y, label = lab),
                             inherit.aes = FALSE, nudge_x = 0.008 * h,
                             nudge_y = SUPPORT_NUDGE_Y, size = 1.65, color = "grey35")
        }
      }
    }
  }

  if (length(brk)) {
    brk_df <- bind_rows(lapply(brk, function(b) {
      e <- b$edge
      pa <- tr_disp$edge[e, 1]; ch <- tr_disp$edge[e, 2]
      x0 <- dd$x[match(pa, dd$node)]; x1 <- dd$x[match(ch, dd$node)]
      y  <- dd$y[match(ch, dd$node)]
      xm <- x0 + (x1 - x0) * 0.55
      dx <- max((x1 - x0) * 0.08, h * 0.004)
      dy <- 0.22
      tibble(xmin = xm - dx, xmax = xm + dx, ymin = y - dy, ymax = y + dy,
             x1a = xm - dx, y1a = y - dy, x2a = xm + dx * 0.3, y2a = y + dy,
             x1b = xm - dx * 0.3, y1b = y - dy, x2b = xm + dx, y2b = y + dy,
             tx = x1, ty = y + dy * 1.45,
             lab = sprintf("branch %.3g", b$true))
    }))
    p <- p +
      geom_rect(data = brk_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE, fill = "white", color = NA) +
      geom_segment(data = brk_df, aes(x = x1a, y = y1a, xend = x2a, yend = y2a),
                   inherit.aes = FALSE, linewidth = 0.45) +
      geom_segment(data = brk_df, aes(x = x1b, y = y1b, xend = x2b, yend = y2b),
                   inherit.aes = FALSE, linewidth = 0.45)
    if (show_break_value) {
      p <- p + geom_text(data = brk_df, aes(x = tx, y = ty, label = lab),
                         inherit.aes = FALSE, size = 1.65, color = "grey30",
                         fontface = "italic", hjust = 0)
    }
  }

  p +
    scale_color_manual(values = pal, breaks = group_levels, labels = group_levels,
                       drop = FALSE,
                       guide = if (do_legend) guide_legend(title = "Original group") else "none") +
    labs(x = NULL, y = NULL) +
    xlim(0, h + h * (0.32 + right_mar / 45)) +
    (if (layout %in% c("rectangular", "roundrect", "slanted")) theme_tree2() else theme_tree()) +
    theme(axis.text.x = element_text(size = 6.5),
          axis.ticks.x = element_line(linewidth = 0.25),
          axis.line.x = element_line(linewidth = 0.25),
          legend.position = legend_position,
          legend.justification = c(1, 1),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.title = element_text(size = legend_title),
          legend.text = element_text(size = legend_text),
          plot.margin = margin(5.5, 5.5 + right_mar * 2.5, 5.5, 5.5)) +
    coord_cartesian(clip = "off")
}

ksnp_tr <- read.tree(F_KSNP)
ksnp_tr$tip.label <- str_extract(ksnp_tr$tip.label, "Cn[0-9A-Za-z]+")
if (anyNA(ksnp_tr$tip.label)) stop("Some tree tip labels do not parse to a CnXXX id.")
expected_tips <- c(names(grp_of), OUTGROUP)
if (!setequal(ksnp_tr$tip.label, expected_tips)) {
  stop("The tree tips do not match the 19 clinical parents plus ", OUTGROUP, ". Missing: ",
       paste(setdiff(expected_tips, ksnp_tr$tip.label), collapse = ", "),
       "; unexpected: ",
       paste(setdiff(ksnp_tr$tip.label, expected_tips), collapse = ", "))
}
missing_st <- setdiff(ksnp_tr$tip.label, names(st_of))
if (length(missing_st)) {
  stop("No MLST sequence type for tree tip(s): ", paste(missing_st, collapse = ", "))
}
ksnp_tr <- root_ladder(ksnp_tr)
ksnp_support <- ksnp_tr$node.label   # aligned with node numbering after rooting

p_tree <- plot_tree_v2(ksnp_tr, support = ksnp_support,
                       layout = "rectangular", branch.length = "none",
                       right_mar = 7.5, tip_size = 2.6, tip_point = 1.5,
                       legend_title = 9, legend_text = 8,
                       legend_position = c(1, 0.98))
panelE <- p_tree


## -----------------------------------------------------------------------------
## 6. Calculated checks and PDF output
## -----------------------------------------------------------------------------
gain_flagged <- Co %>% filter(gain >= GAIN_THRESHOLD) %>% pull(sample)
hetz_flagged <- Co %>% filter(hetz >= asinh(HETZ_THRESHOLD)) %>% pull(sample)
st_table <- tibble(tip = ksnp_tr$tip.label,
                   st = unname(st_of[ksnp_tr$tip.label])) %>%
  mutate(st_disp = ifelse(st %in% ST_NOVEL_FROM, "ST novel", paste0("ST", st))) %>%
  count(st_disp, name = "n_tips") %>% arrange(desc(n_tips))

message(
  "Figure 6: ", nrow(Co), " Parent isolates in Panel A; ",
  nrow(dA), " clone states in Panel C; ",
  length(ksnp_tr$tip.label), " tips in Panel E."
)
message("Panel A: isolates at or above the gain threshold (", GAIN_THRESHOLD, "): ",
        if (length(gain_flagged)) paste(gain_flagged, collapse = ", ") else "none")
message("Panel A: isolates at or above the reference-het z threshold (", HETZ_THRESHOLD, "): ",
        if (length(hetz_flagged)) paste(hetz_flagged, collapse = ", ") else "none")
message("Panel E: sequence types across the tree tips")
print(st_table, n = Inf)

if (EXPORT_PANELS) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(OUTPUT_DIR)) stop("Cannot create output directory: ", OUTPUT_DIR)
  if (!capabilities("cairo")) stop("PDF output requires Cairo support in R.")
  panel_plots <- list(Fig6A = panelA, Fig6C = panelC, Fig6D = panelD, Fig6E = panelE)
  panel_sizes <- list(Fig6A = c(4.4, 3.7), Fig6C = c(5.0, 3.7),
                      Fig6D = c(6.3, 4.6), Fig6E = c(4.4, 3.8))
  iwalk(panel_plots, function(p, name) {
    size <- panel_sizes[[name]]
    ggsave(
      file.path(OUTPUT_DIR, paste0(name, ".pdf")), p,
      width = size[1], height = size[2], units = "in", device = cairo_pdf
    )
  })
  message("Figure 6 PDFs written to: ", normalizePath(OUTPUT_DIR))
}
