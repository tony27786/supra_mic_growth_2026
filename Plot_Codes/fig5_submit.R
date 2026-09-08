# Figure 5 plot code
#
# Inputs: ../SourceData. Default PDF output: ../Figure5_panels.
# Set FIG5_EXPORT=false to calculate without writing PDFs.
# FIG5_OUTPUT_DIR can redirect PDF output; no input files are modified.

suppressPackageStartupMessages({
  library(data.table)
  library(Biostrings)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  # data.table and Biostrings mask a number of tidyverse verbs, so the tidy core
  # is attached last and wins on the search path.
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
  found <- dir_candidates[file.exists(file.path(dir_candidates, "fig5_submit.R"))]
  if (length(found) != 1L) {
    stop("Run fig5_submit.R with source() or Rscript so its directory can be resolved.")
  }
  normalizePath(found, mustWork = TRUE)
}
SCRIPT_DIR <- get_script_dir()
SOURCE_DIR <- normalizePath(file.path(SCRIPT_DIR, "..", "SourceData"), mustWork = TRUE)
DERIVED_DIR <- Sys.getenv("FIG5_DERIVED_DIR", file.path(SCRIPT_DIR, "..", "Derived"))
OUTPUT_DIR <- Sys.getenv("FIG5_OUTPUT_DIR", file.path(SCRIPT_DIR, "..", "Figure5_panels"))
EXPORT_PANELS <- tolower(Sys.getenv("FIG5_EXPORT", "true")) %in% c("1", "true", "yes")

F_MASTER <- file.path(SOURCE_DIR, "master_wgs_table.tsv")
F_DEPTH  <- file.path(SOURCE_DIR, "mosdepth_chr_normalized_depth_long.tsv")
F_COUNTS <- file.path(SOURCE_DIR, "gene_counts.txt")
F_SMV    <- file.path(SOURCE_DIR, "smallvar_4state_long.tsv")
F_CALLABILITY <- file.path(SOURCE_DIR, "target_panel_callability.tsv")
F_GTF    <- file.path(SOURCE_DIR, "refs", "GCF_000149245.1_CNA3_genomic.gtf")
F_FNA    <- file.path(SOURCE_DIR, "refs", "GCF_000149245.1_CNA3_genomic.fna")
F_GAF    <- file.path(SOURCE_DIR, "refs", "FungiDB-68_CneoformansH99_GO.gaf.gz")
F_GAF_ALL <- file.path(SOURCE_DIR, "refs",
                       c("FungiDB-68_CneoformansH99_GO.gaf.gz",
                         "FungiDB-68_CneoformansH99_Curated_GO.gaf.gz"))
source_inputs <- c(F_MASTER, F_DEPTH, F_COUNTS, F_SMV, F_CALLABILITY,
                   F_GTF, F_FNA, F_GAF_ALL)
if (any(!file.exists(source_inputs))) {
  stop("Missing SourceData input(s): ",
       paste(basename(source_inputs[!file.exists(source_inputs)]), collapse = ", "))
}

F_DE3    <- file.path(DERIVED_DIR, "fig3_DE_batch01_vs_02.csv")
F_DE4    <- file.path(DERIVED_DIR, "fig4_DE_batch03_vs_01.csv")
F_PHENO  <- file.path(DERIVED_DIR, "fig4_master_out.csv")
derived_inputs <- c(F_DE3, F_DE4, F_PHENO)
if (any(!file.exists(derived_inputs))) {
  stop("Missing derived input(s): ",
       paste(basename(derived_inputs[!file.exists(derived_inputs)]), collapse = ", "),
       "\nRun fig3_submit.R and fig4_submit.R first; they write these into ",
       DERIVED_DIR, ".")
}

require_columns <- function(d, columns, label) {
  missing <- setdiff(columns, names(d))
  if (length(missing)) stop(label, " is missing: ", paste(missing, collapse = ", "))
}


## -----------------------------------------------------------------------------
## Shared conventions
## -----------------------------------------------------------------------------

# Integrated WGS genome-state classes carried over from the genome-state call.
CLASS_COL <- c(monoclonal_clean = "#009E73", monoclonal_with_disomy = "#E69F00",
               heterozygous_diploid = "#0072B2", divergent_clone_artifact = "#D55E00",
               divergent_clone_artifact_marginal = "#CC79A7",
               ambiguous_needs_review = "#999999", reference_H99 = "#000000")
CLASS_LAB <- c(monoclonal_clean = "clean haploid", monoclonal_with_disomy = "haploid + disomy",
               heterozygous_diploid = "heterozygous diploid",
               divergent_clone_artifact = "divergent (reference artifact)",
               divergent_clone_artifact_marginal = "divergent (marginal)",
               ambiguous_needs_review = "ambiguous", reference_H99 = "H99 reference")

# Chromosome accessions: NC_026745.1 .. NC_026758.1 = Chr1..Chr14. Chr2 is kept
CHR_ACC <- sprintf("NC_0267%02d.1", 45:58)
CHR_NUM <- setNames(seq_along(CHR_ACC), CHR_ACC)
CHR2 <- "NC_026746.1"

# One cohort definition serves both halves of this script: 19 clinical isolates in phenotype-group order, then the H99 reference.
COHORT <- tribble(
  ~strain, ~group_code,
  "101", "HR",  "106", "HR",  "126", "HR",  "128", "HR",  "160", "HR",
  "162", "HR",  "208", "HR",  "213", "HR",  "215", "HR",  "217", "HR",
  "147", "control", "150", "control", "169", "control",
  "174", "control", "225", "control", "226", "control",
  "196", "R", "197", "R", "198", "R"
)
STR_ORDER <- c(COHORT$strain, "H99")

# CnXXXBB -> strain XXX and state BB; a bare CnXXX is the Parent bulk.
parse_batch  <- function(s){ b <- str_remove(s, "^Cn")
  ifelse(str_length(b) > 3 & str_sub(b, -2) %in% c("01","02","03"), str_sub(b, -2), "parent") }
parse_strain <- function(s){ b <- str_remove(s, "^Cn")
  ifelse(str_length(b) > 3 & str_sub(b, -2) %in% c("01","02","03"), str_sub(b, 1, -3), b) }
BATCH_LAB <- c(parent = "Parent", `02` = "Monoclonal", `01` = "Selected", `03` = "Passage")

# Lollipop switches for Panel C.
FILTER_HYPOTHETICAL_PROTEIN <- TRUE
PD_GENES_PER_CHR            <- 6L


## -----------------------------------------------------------------------------
## 1. WGS depth and genome-state load
## -----------------------------------------------------------------------------
master <- read_tsv(F_MASTER, show_col_types = FALSE)
require_columns(master, c("sample", "call_wgs_class"), basename(F_MASTER))
depth  <- read_tsv(F_DEPTH, show_col_types = FALSE)
require_columns(depth, c("sample", "chrom", "chr_depth_norm"), basename(F_DEPTH))
class_of <- setNames(master$call_wgs_class, master$sample)

# Cohort-relative chromosome ratio: divide each chromosome by its cohort median across all sequenced samples.
cohort_med <- depth %>% filter(chrom %in% CHR_ACC) %>% group_by(chrom) %>%
  summarise(m = median(chr_depth_norm, na.rm = TRUE), .groups = "drop")
depth <- depth %>% left_join(cohort_med, by = "chrom") %>% mutate(rel = chr_depth_norm / m)


## -----------------------------------------------------------------------------
## 2. Panel A: copy-number landscape (p_cnv_hm)
## -----------------------------------------------------------------------------
cnv <- depth %>%
  filter(chrom %in% CHR_ACC) %>%
  mutate(strain = parse_strain(sample), batch = parse_batch(sample)) %>%
  filter(strain %in% STR_ORDER) %>%
  group_by(sample) %>%
  mutate(rel_hm = rel / median(rel[chrom != CHR2], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(strain = factor(strain, levels = STR_ORDER),
         batch  = factor(BATCH_LAB[batch], levels = BATCH_LAB),
         chr    = factor(paste0("Chr", CHR_NUM[chrom]),
                         levels = rev(paste0("Chr", CHR_NUM[CHR_ACC]))))
if (!nrow(cnv)) stop("No cohort samples found in ", basename(F_DEPTH), ".")

p_cnv_hm <- ggplot(cnv, aes(strain, chr, fill = rel_hm)) +
  geom_tile() +
  annotate("rect", xmin = 0.5, xmax = length(STR_ORDER) + 0.5,
           ymin = which(levels(cnv$chr) == "Chr1") - 0.5,
           ymax = which(levels(cnv$chr) == "Chr1") + 0.5,
           fill = NA, colour = "grey10", linewidth = 0.5) +
  annotate("rect", xmin = 0.5, xmax = length(STR_ORDER) + 0.5,
           ymin = which(levels(cnv$chr) == "Chr4") - 0.5,
           ymax = which(levels(cnv$chr) == "Chr4") + 0.5,
           fill = NA, colour = "grey10", linewidth = 0.5) +
  facet_wrap(~ batch, nrow = 1) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 1,
                       guide = guide_colourbar(direction = "horizontal"),
                       limits = c(0.75, 2), oob = squish, name = "copy ratio") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.5, face = "bold"),
        axis.text.y = element_text(size = 6.5, face = "bold"), panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 8.5), legend.position = "top",
        legend.key.height = unit(4, "mm"))
panelA <- p_cnv_hm


## -----------------------------------------------------------------------------
## 3. RNA-seq setup for Panels B and C
## -----------------------------------------------------------------------------
# Median-of-ratios size factors, then gene-centred log2. The Cn160 libraries are dropped because only 41-45% of their reads map to H99.
RNA_EXCLUDE <- c("Cn16001", "Cn16002", "Cn16003")
gc_raw <- read_tsv(F_COUNTS, comment = "#", show_col_types = FALSE)
require_columns(gc_raw, c("Geneid", "Chr"), basename(F_COUNTS))
meta_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
samp_cols <- setdiff(names(gc_raw), meta_cols)
clone_ids <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", basename(samp_cols))
names(gc_raw)[match(samp_cols, names(gc_raw))] <- clone_ids
clone_ids <- setdiff(clone_ids, RNA_EXCLUDE)
if (!length(clone_ids)) stop("No usable RNA-seq libraries in ", basename(F_COUNTS), ".")
gene_chr  <- setNames(sub(";.*", "", gc_raw$Chr), gc_raw$Geneid)
cnts <- as.matrix(gc_raw[, clone_ids]); rownames(cnts) <- gc_raw$Geneid
gkeep  <- rowSums(cnts == 0) == 0
logref <- rowMeans(log(cnts[gkeep, , drop = FALSE]))
sf     <- apply(cnts[gkeep, , drop = FALSE], 2, function(col) median(exp(log(col) - logref)))
norm   <- sweep(cnts, 2, sf, "/")
l2     <- log2(norm + 1)
l2c    <- sweep(l2, 1, apply(l2, 1, median), "-")


## -----------------------------------------------------------------------------
## 4. Panel B: DNA dosage against RNA dosage (pC)
## -----------------------------------------------------------------------------
cr <- depth %>% filter(sample %in% clone_ids, chrom %in% CHR_ACC) %>%
  select(sample, acc = chrom, copyratio = chr_depth_norm)
shift <- bind_rows(lapply(CHR_ACC, function(acc) {
  g <- intersect(names(gene_chr)[gene_chr == acc], rownames(l2c))
  if (!length(g)) return(NULL)
  data.frame(acc = acc, sample = colnames(l2c),
             exprshift = apply(l2c[g, , drop = FALSE], 2, median), row.names = NULL)
}))
Cdat <- inner_join(cr, shift, by = c("sample", "acc")) %>%
  mutate(cls = factor(dplyr::recode(class_of[sample],
                                    divergent_clone_artifact_marginal = "divergent_clone_artifact"),
                      levels = names(CLASS_COL)),
         is_chr2 = acc == CHR2)
if (!nrow(Cdat)) stop("No sample-chromosome pairs with both WGS depth and RNA-seq.")
rC_all <- cor(Cdat$copyratio, Cdat$exprshift, use = "complete.obs")
rC_no2 <- with(filter(Cdat, !is_chr2), cor(copyratio, exprshift, use = "complete.obs"))
dose <- data.frame(copyratio = seq(0.7, 2.0, length.out = 100)) %>%
  mutate(exprshift = log2(copyratio))

pC <- ggplot(Cdat, aes(copyratio, exprshift)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 1, colour = "grey85", linewidth = 0.3) +
  geom_line(data = dose, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_point(aes(color = cls), data = filter(Cdat, !is_chr2), size = 1.8, alpha = 0.75) +
  geom_point(data = filter(Cdat, is_chr2), shape = 4, size = 2, colour = "grey35", stroke = 0.6) +
  annotate("text", x = 1.22, y = -0.55, size = 2.5, colour = "grey35",
           label = "Chr2 (×): depth-inflated,\nexpression-flat → artifact") +
  scale_color_manual(values = CLASS_COL, labels = CLASS_LAB,
                     name = "WGS genome state", drop = TRUE) +
  labs(x = "WGS chromosome copy ratio (normalized depth)",
       y = "RNA chromosome expression shift (median gene-centered log2)") +
  theme_classic(base_size = 11) +
  theme(legend.key.size = unit(3.2, "mm"), legend.text = element_text(size = 7),
        legend.title = element_text(size = 8))
panelB <- pC


## -----------------------------------------------------------------------------
## 5. Panel C: copy-number-associated genes on Chr1 and Chr4 (pD_lolli)
## -----------------------------------------------------------------------------
CN <- depth %>% filter(chrom %in% CHR_ACC, sample %in% clone_ids) %>%
  select(sample, chrom, chr_depth_norm) %>%
  pivot_wider(names_from = sample, values_from = chr_depth_norm)
CNm <- as.matrix(CN[, clone_ids]); rownames(CNm) <- CN$chrom
chr_spread <- apply(CNm, 1, function(x) diff(range(x, na.rm = TRUE)))
ANEU_ACC   <- names(chr_spread)[chr_spread > 0.4]

r_dosage <- vapply(rownames(l2), function(g) {
  acc <- gene_chr[g]
  if (is.na(acc) || !acc %in% rownames(CNm)) return(NA_real_)
  cn <- CNm[acc, clone_ids]
  if (sd(cn, na.rm = TRUE) < 1e-6) return(NA_real_)
  suppressWarnings(cor(l2[g, clone_ids], cn, use = "complete.obs"))
}, numeric(1))

# de3 = Selected vs Monoclonal (Figure 3G); de4 = Passage vs Selected (Figure 4G).
de3 <- read_csv(F_DE3, show_col_types = FALSE) %>% dplyr::rename(gene = gene_id)
de4 <- read_csv(F_DE4, show_col_types = FALSE) %>% dplyr::rename(gene = gene_id)
require_columns(de3, c("gene", "log2FoldChange", "padj", "symbol", "product"), basename(F_DE3))
require_columns(de4, c("gene", "log2FoldChange", "padj"), basename(F_DE4))
rev_set <- inner_join(de3, de4, by = "gene", suffix = c("_3", "_4")) %>%
  filter(padj_3 < 0.05, padj_4 < 0.05,
         sign(log2FoldChange_3) != sign(log2FoldChange_4),
         abs(log2FoldChange_3) > 0.5) %>% pull(gene)

Dg <- tibble(gene = rownames(l2), r_dosage = r_dosage) %>%
  mutate(acc = gene_chr[gene], chr = CHR_NUM[acc], aneu = acc %in% ANEU_ACC) %>%
  left_join(select(de3, gene, log2FoldChange, padj, symbol, product), by = "gene") %>%
  mutate(reversible = gene %in% rev_set) %>%
  filter(!is.na(r_dosage), !is.na(log2FoldChange))

AZOLE <- tibble(gene = c("CNAG_00040", "CNAG_00730", "CNAG_00854"),
                azlab = c("ERG11*", "AFR1*", "ERG2*"))
Dg <- Dg %>% mutate(cand = aneu & abs(r_dosage) > 0.5 & !is.na(padj) & padj < 0.05 &
                      abs(log2FoldChange) > 0.5 & reversible)
glab <- function(sym, gene, prod) {
  s <- ifelse(!is.na(sym), sym, str_remove(gene, "CNAG_"))
  paste0(s, "  ·  ", str_trunc(ifelse(is.na(prod), "", prod), 34))
}

Dg_lolli <- Dg %>%
  mutate(is_hypothetical = str_detect(coalesce(product, ""),
                                      regex("hypothetical\\s+protein", ignore_case = TRUE))) %>%
  filter(!FILTER_HYPOTHETICAL_PROTEIN | !is_hypothetical)

chr1_controls <- Dg_lolli %>%
  filter(gene %in% AZOLE$gene) %>%
  arrange(match(gene, AZOLE$gene)) %>%
  slice_head(n = PD_GENES_PER_CHR)
chr1_n_fill <- max(PD_GENES_PER_CHR - nrow(chr1_controls), 0L)
chr1_fill <- Dg_lolli %>%
  filter(cand, chr == 1, !gene %in% chr1_controls$gene) %>%
  arrange(desc(r_dosage), padj, gene) %>%
  slice_head(n = chr1_n_fill)
chr1_pick <- bind_rows(chr1_controls, chr1_fill) %>%
  slice_head(n = PD_GENES_PER_CHR) %>%
  mutate(grp = "Chr1")

chr4_pick <- Dg_lolli %>%
  filter(cand, chr == 4) %>%
  arrange(desc(r_dosage), padj, gene) %>%
  slice_head(n = PD_GENES_PER_CHR) %>%
  mutate(grp = "Chr4")

if (nrow(chr1_pick) < PD_GENES_PER_CHR || nrow(chr4_pick) < PD_GENES_PER_CHR) {
  warning("Panel C lollipop has fewer than ", PD_GENES_PER_CHR,
          " entries in at least one chromosome after filtering.")
}

pick <- bind_rows(chr1_pick, chr4_pick) %>%
  mutate(grp = factor(grp, levels = c("Chr1", "Chr4")),
         lab = glab(symbol, gene, product)) %>%
  arrange(grp, r_dosage) %>% mutate(gene = factor(gene, levels = gene))
lab_map <- setNames(pick$lab, pick$gene)

pD_lolli <- ggplot(pick, aes(r_dosage, gene)) +
  geom_segment(aes(x = 0, xend = r_dosage, yend = gene), colour = "grey80", linewidth = 0.5) +
  geom_vline(xintercept = 0.5, linetype = "dotted", colour = "grey60", linewidth = 0.4) +
  geom_point(aes(fill = log2FoldChange, size = -log10(padj)), shape = 21, colour = "grey20") +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_y_discrete(labels = lab_map) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-2, 2), oob = squish,
                       name = "log2FC\n(Selected vs\nMonoclonal)") +
  scale_size_continuous(range = c(1.6, 5), name = "-log10 padj") +
  labs(x = "Expression-copy-number correlation (Pearson r)", y = NULL) +
  theme_bw(base_size = 10) +
  theme(strip.placement = "outside", strip.background = element_blank(),
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 7.2),
        panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 7))
panelC <- pD_lolli


## -----------------------------------------------------------------------------
## 6. Panel E setup: frozen small-variant conventions
## -----------------------------------------------------------------------------
VAF_ABSENT   <- 0.05
VAF_SELECTED <- 0.8
VAF_RETAINED <- 0.5
DP_MIN       <- 10L
EXPECTED_N_CANDIDATES        <- 13L
EXPECTED_N_RETAINED_EVENTS   <- 9L
EXPECTED_N_RETAINED_CARRIERS <- 6L

SV_WORKING <- data.table(
  strain = STR_ORDER,
  group_code = c(COHORT$group_code, "H99")
)
CLINICAL_LEVELS <- SV_WORKING[group_code != "H99", strain]
STR_LEVELS <- SV_WORKING$strain
GROUP_COLORS <- c(group_hr = "#E69F00", group_control = "#009E73",
                  group_resistant = "#D55E00", group_h99 = "#333333")

# Figure 5 phenotype classes are fixed independently of the small-variant calls.
LIMITED_COMPOSITE <- c("106", "126", "169")
KB_SPECIFIC       <- "217"
KB_ONLY_SET       <- c("106", "126", KB_SPECIFIC)
PHENOTYPE_SET     <- c(LIMITED_COMPOSITE, KB_SPECIFIC)

## 6a. Phenotype cross-check against the Figure 4 master table
pheno <- fread(F_PHENO)
pheno[, strain := as.character(strain)]
needed_pheno_cols <- c("strain", "group", "d_sel_kb", "d_rev_kb", "d_sel_pap", "d_rev_pap",
                       "kb_parent", "kb_clone", "kb_passage",
                       "pap_parent", "pap_clone", "pap_passage")
if (!all(needed_pheno_cols %in% names(pheno))) {
  stop("The Figure 4 master table is missing required columns: ",
       paste(setdiff(needed_pheno_cols, names(pheno)), collapse = ", "))
}
pheno[, group := fcase(grepl("^HR", group, ignore.case = TRUE), "HR",
                       grepl("^Control", group, ignore.case = TRUE), "Control",
                       grepl("^R", group, ignore.case = TRUE), "R",
                       default = "H99")]

s_kb  <- sd(c(pheno$d_sel_kb, pheno$d_rev_kb), na.rm = TRUE)
s_pap <- sd(c(pheno$d_sel_pap, pheno$d_rev_pap), na.rm = TRUE)
pheno[, `:=`(z_sel_kb = d_sel_kb / s_kb, z_rev_kb = d_rev_kb / s_kb,
             z_sel_pap = d_sel_pap / s_pap, z_rev_pap = d_rev_pap / s_pap)]
pheno[, sel_mag := sqrt(z_sel_kb^2 + z_sel_pap^2)]
pheno[, recovery := -(z_rev_kb * z_sel_kb + z_rev_pap * z_sel_pap) / pmax(sel_mag, 1e-06)]

lowest_three <- pheno[group %chin% c("HR", "Control")][order(recovery), head(strain, 3)]
if (!setequal(lowest_three, LIMITED_COMPOSITE)) {
  stop("The recovery ranking drifted. Expected lowest composite recovery: ",
       paste(LIMITED_COMPOSITE, collapse = ", "), "; observed: ",
       paste(lowest_three, collapse = ", "))
}
cn217 <- pheno[strain == KB_SPECIFIC]
if (nrow(cn217) != 1L || !(cn217$kb_passage <= cn217$kb_clone &&
                           cn217$pap_passage < cn217$pap_parent)) {
  stop("Cn217 no longer satisfies the frozen KB-specific persistence context.")
}
h99_pheno <- pheno[strain == "H99"]
if (nrow(h99_pheno) != 1L || !(abs(h99_pheno$kb_passage - h99_pheno$kb_parent) < 1 &&
                               h99_pheno$pap_passage < h99_pheno$pap_clone)) {
  stop("H99 no longer satisfies the frozen drug-free recovery-control context.")
}


## 6b. Four-state small-variant calls
message("Reading the four-state small-variant table (large file, this takes a while)...")
smv <- fread(F_SMV, sep = "\t", showProgress = TRUE,
             select = c("strain", "chrom", "pos", "ref", "alt", "type",
                        "qual", "filter", "state", "dp", "vaf"))
smv[, `:=`(strain = as.character(strain), pos = as.integer(pos),
           dp = as.integer(dp), vaf = suppressWarnings(as.numeric(vaf)))]

# Candidate calls shown in the figure must pass the caller's site filter.
smv <- smv[filter == "PASS"]

wide <- dcast(smv, strain + chrom + pos + ref + alt + type + qual + filter ~ state,
              value.var = c("vaf", "dp"))
rm(smv)
invisible(gc())

rename_map <- c(vaf_parent = "v_par", vaf_01 = "v_01", vaf_02 = "v_02", vaf_03 = "v_03",
                dp_parent = "d_par", dp_01 = "d_01", dp_02 = "d_02", dp_03 = "d_03")
for (old in names(rename_map)) {
  if (!old %in% names(wide)) stop("Expected state column is missing: ", old)
  setnames(wide, old, rename_map[[old]])
}

scored_four_state <- wide[!is.na(v_par) & !is.na(v_01) & !is.na(v_02) & !is.na(v_03) &
                          d_par >= DP_MIN & d_01 >= DP_MIN & d_02 >= DP_MIN & d_03 >= DP_MIN]

candidates <- scored_four_state[v_par < VAF_ABSENT & v_02 < VAF_ABSENT & v_01 >= VAF_SELECTED]
candidates[, retained := v_03 >= VAF_RETAINED]
candidates[, site_id := paste(chrom, pos, ref, alt, sep = ":")]

if (nrow(candidates) != EXPECTED_N_CANDIDATES) {
  stop("Candidate count drifted: expected ", EXPECTED_N_CANDIDATES,
       ", observed ", nrow(candidates))
}
if (sum(candidates$retained) != EXPECTED_N_RETAINED_EVENTS) {
  stop("Retained-event count drifted: expected ", EXPECTED_N_RETAINED_EVENTS,
       ", observed ", sum(candidates$retained))
}
retained_carriers <- sort(unique(candidates[retained == TRUE, strain]))
if (length(retained_carriers) != EXPECTED_N_RETAINED_CARRIERS) {
  stop("Retained-carrier count drifted: expected ", EXPECTED_N_RETAINED_CARRIERS,
       ", observed ", length(retained_carriers))
}


## 6c. Reference annotation
gtf <- read_tsv(F_GTF, comment = "#", col_names = FALSE, show_col_types = FALSE,
                col_types = cols(.default = col_character()))
if (ncol(gtf) < 9L) stop("The GTF did not parse into nine columns.")
gtf <- gtf[, 1:9]
setDT(gtf)
setnames(gtf, c("chrom", "source", "feature", "start", "end",
                "score", "strand", "phase", "attribute"))
gtf[, `:=`(start = as.integer(start), end = as.integer(end),
           gene = str_match(attribute, "gene_id \"([^\"]+)\"")[, 2])]
gtf <- gtf[!is.na(gene)]
gtf_gene <- gtf[feature == "gene", .(chrom, start, end, strand, gene)]
gtf_cds  <- gtf[feature == "CDS", .(chrom, start, end, strand, gene, attribute)]

gene_at <- function(chrom_value, pos_value, features) {
  hit <- features[chrom == chrom_value & start <= pos_value & end >= pos_value, gene]
  if (length(hit)) hit[1] else NA_character_
}
candidates[, gene_context := mapply(gene_at, chrom, pos,
                                    MoreArgs = list(features = gtf_gene), USE.NAMES = FALSE)]
candidates[, gene_cds := mapply(gene_at, chrom, pos,
                                MoreArgs = list(features = gtf_cds), USE.NAMES = FALSE)]
candidates[, gene := fifelse(!is.na(gene_cds), gene_cds, gene_context)]

gaf_lines <- readLines(gzfile(F_GAF))
gaf_lines <- gaf_lines[!startsWith(gaf_lines, "!")]
gaf <- fread(text = paste(gaf_lines, collapse = "\n"), sep = "\t",
             header = FALSE, fill = TRUE, quote = "", showProgress = FALSE)
if (ncol(gaf) < 10L) stop("The GO annotation did not parse into the expected columns.")
gaf_map <- gaf[!is.na(V2) & V2 != "", {
  symbols <- unique(V3[!is.na(V3) & V3 != ""])
  preferred <- symbols[symbols != V2[1]]
  products <- unique(V10[!is.na(V10) & V10 != ""])
  list(symbol = if (length(preferred)) preferred[1] else V2[1],
       product = if (length(products)) products[1] else NA_character_)
}, by = .(gene = V2)]

candidates <- merge(candidates, gaf_map, by = "gene", all.x = TRUE, sort = FALSE)
candidates[is.na(symbol) & !is.na(gene), symbol := gene]
candidates[is.na(symbol), symbol := "non-CDS"]

candidates[, is_indel := type != "SNP" | nchar(ref) != nchar(alt)]
candidates[, event_class := fcase(is_indel & !is.na(gene_cds), "coding_indel",
                                  is_indel, "noncoding_indel",
                                  !is.na(gene_cds), "coding_snv",
                                  default = "noncoding_snv")]
candidates[, chr := paste0("Chr", CHR_NUM[chrom])]
candidates[is.na(CHR_NUM[chrom]), chr := chrom]
candidates[, site_label := paste0(symbol, "\n", chr, ":", pos, " ", ref, ">", alt)]


## 6d. Bounded 25-gene negative result and callability
.read_gaf_sym <- function(f){
  ln <- readLines(gzfile(f)); ln <- ln[!startsWith(ln, "!")]
  x  <- strsplit(ln, "\t", fixed = TRUE)
  tibble(gene_id = vapply(x, `[`, "", 2L), symbol = vapply(x, `[`, "", 3L)) %>%
    filter(startsWith(gene_id, "CNAG"), symbol != "", symbol != gene_id)
}
cnag2sym <- map_dfr(F_GAF_ALL, .read_gaf_sym) %>% distinct(gene_id, .keep_all = TRUE)
if (!nrow(cnag2sym)) stop("No symbol/CNAG pairs parsed from the GO annotation files.")

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
  dplyr::mutate(cnag = dplyr::coalesce(unname(.manual[symbol]), gene_id),
                tier = factor(tier, levels = c("dosage", "output", "context")))
.expected <- c(HSP78 = "CNAG_03347", HSP90 = "CNAG_06150", CHS6 = "CNAG_06487")
stopifnot(nrow(target_panel) == 25L, !anyNA(target_panel$cnag),
          !anyDuplicated(target_panel$symbol), !anyDuplicated(target_panel$cnag),
          identical(unname(target_panel$cnag[match(names(.expected), target_panel$symbol)]),
                    unname(.expected)),
          !("HSP104" %in% target_panel$symbol), !("CHS4" %in% target_panel$symbol))

n_panel_coding_candidates <- candidates[!is.na(gene_cds) &
                                        gene_cds %chin% target_panel$cnag, .N]
if (n_panel_coding_candidates != 0L) {
  stop("A qualifying coding candidate now overlaps the frozen 25-gene panel.")
}

callability <- fread(F_CALLABILITY)
if (!"common_frac_cds" %in% names(callability)) {
  stop("The callability table lacks common_frac_cds.")
}
callability_full <- sum(callability$common_frac_cds == 1, na.rm = TRUE)
callability_min  <- min(callability$common_frac_cds, na.rm = TRUE)


## 6e. Parent ERG11 p.G484S
ERG11 <- target_panel$cnag[target_panel$symbol == "ERG11"]
if (length(ERG11) != 1L) stop("ERG11 did not resolve to one frozen CNAG ID.")

genome <- readDNAStringSet(F_FNA)
names(genome) <- sub(" .*", "", names(genome))

erg11_exons <- gtf_cds[gene == ERG11][order(start)]
if (!nrow(erg11_exons)) stop("The ERG11 CDS was not found in the GTF.")
erg11_chrom  <- erg11_exons$chrom[1]
erg11_strand <- erg11_exons$strand[1]

erg11_seq <- do.call(c, lapply(seq_len(nrow(erg11_exons)), function(i)
  subseq(genome[[erg11_exons$chrom[i]]], erg11_exons$start[i], erg11_exons$end[i])))
erg11_coords <- unlist(Map(`:`, erg11_exons$start, erg11_exons$end))
if (erg11_strand == "-") {
  erg11_seq <- reverseComplement(erg11_seq)
  erg11_coords <- rev(erg11_coords)
}
erg11_index <- setNames(seq_along(erg11_coords), as.character(erg11_coords))

annotate_erg11_snv <- function(pos_value, ref_value, alt_value) {
  idx <- unname(erg11_index[as.character(pos_value)])
  if (is.na(idx) || nchar(ref_value) != 1L || nchar(alt_value) != 1L) return(NA_character_)
  if (erg11_strand == "-") {
    ref_value <- as.character(reverseComplement(DNAString(ref_value)))
    alt_value <- as.character(reverseComplement(DNAString(alt_value)))
  }
  codon_index  <- (idx - 1L) %/% 3L
  codon_offset <- (idx - 1L) %% 3L
  codon_start  <- codon_index * 3L + 1L
  if (codon_start + 2L > length(erg11_seq)) return(NA_character_)
  ref_codon <- as.character(subseq(erg11_seq, codon_start, codon_start + 2L))
  if (substr(ref_codon, codon_offset + 1L, codon_offset + 1L) != ref_value) {
    stop("Reference FASTA mismatch while annotating ERG11 at position ", pos_value)
  }
  alt_codon <- ref_codon
  substr(alt_codon, codon_offset + 1L, codon_offset + 1L) <- alt_value
  aa_ref <- as.character(translate(DNAString(ref_codon), no.init.codon = TRUE))
  aa_alt <- as.character(translate(DNAString(alt_codon), no.init.codon = TRUE))
  paste0("p.", aa_ref, codon_index + 1L, aa_alt)
}

erg11_parent <- wide[chrom == erg11_chrom & d_par >= DP_MIN & !is.na(v_par) &
                     v_par >= VAF_RETAINED & pos >= min(erg11_exons$start) &
                     pos <= max(erg11_exons$end) & nchar(ref) == 1L & nchar(alt) == 1L]
erg11_parent[, protein_change := mapply(annotate_erg11_snv, pos, ref, alt, USE.NAMES = FALSE)]
erg11_parent <- erg11_parent[protein_change == "p.G484S"]

g484s_carriers <- sort(unique(erg11_parent$strain))
if (!identical(g484s_carriers, c("196", "197", "198"))) {
  stop("The ERG11 p.G484S carrier set drifted: ", paste(g484s_carriers, collapse = ", "))
}


## 6f. Strain-level exploratory association
carrier_table <- function(case_set, carriers = retained_carriers, label = "cases") {
  ctrl_set <- setdiff(CLINICAL_LEVELS, case_set)
  matrix(c(sum(case_set %chin% carriers), sum(!case_set %chin% carriers),
           sum(ctrl_set %chin% carriers), sum(!ctrl_set %chin% carriers)),
         nrow = 2, byrow = TRUE,
         dimnames = list(phenotype = c(label, "remaining"),
                         retained_candidate = c("yes", "no")))
}
association_table <- carrier_table(LIMITED_COMPOSITE, label = "limited_KB_PAP_recovery")
expected_association_table <- matrix(c(3, 0, 3, 13), 2, byrow = TRUE)
if (!all(dim(association_table) == dim(expected_association_table)) ||
    !all(unname(association_table) == expected_association_table)) {
  stop("The carrier-by-phenotype table drifted from 3/0 versus 3/13. Observed: ",
       paste(as.vector(t(association_table)), collapse = ", "))
}
fisher_primary <- fisher.test(association_table, alternative = "two.sided")
kb_only_table  <- carrier_table(KB_ONLY_SET, label = "limited_KB_only_recovery")
fisher_kb_only <- fisher.test(kb_only_table, alternative = "two.sided")$p.value

retention_sensitivity <- rbindlist(lapply(c(0.3, 0.5, 0.8), function(threshold) {
  carriers_t <- unique(candidates[v_03 >= threshold, strain])
  tab_t <- carrier_table(LIMITED_COMPOSITE, carriers = carriers_t)
  data.table(retained_vaf_threshold = threshold,
             phenotype_yes = tab_t[1, 1], phenotype_no = tab_t[1, 2],
             remaining_yes = tab_t[2, 1], remaining_no = tab_t[2, 2],
             fisher_two_sided_p = fisher.test(tab_t)$p.value)
}))


## -----------------------------------------------------------------------------
## 7. Panel E: small-variant OncoPrint
## -----------------------------------------------------------------------------
candidate_strain_order <- c(LIMITED_COMPOSITE, KB_SPECIFIC,
                            setdiff(retained_carriers, PHENOTYPE_SET),
                            setdiff(unique(candidates$strain), retained_carriers))
candidates[, strain_rank := match(strain, candidate_strain_order)]
setorder(candidates, -retained, strain_rank, chrom, pos)

event_levels <- unique(candidates$site_id)
event_label_map <- setNames(candidates[match(event_levels, site_id), site_label], event_levels)

BLOCK_LEVELS <- c("Cohort", "ERG11", "Selected-clone-associated candidates", "Recovery")
column_meta <- rbindlist(list(
  data.table(col_id = "cohort_group", block = "Cohort", label = "Group"),
  data.table(col_id = "ERG11_p.G484S", block = "ERG11", label = "ERG11\np.G484S"),
  data.table(col_id = event_levels, block = "Selected-clone-associated candidates",
             label = unname(event_label_map[event_levels])),
  data.table(col_id = "fig5_phenotype", block = "Recovery", label = "Phenotype")))
column_meta[, block := factor(block, levels = BLOCK_LEVELS)]
column_levels <- column_meta$col_id

background <- CJ(strain = STR_LEVELS, col_id = column_levels, unique = TRUE)
background <- merge(background, column_meta, by = "col_id", all.x = TRUE, sort = FALSE)

group_cells <- copy(SV_WORKING)
group_cells[, `:=`(
  col_id = "cohort_group",
  fill_class = fcase(group_code == "HR", "group_hr",
                     group_code == "control", "group_control",
                     group_code == "R", "group_resistant",
                     default = "group_h99"),
  cell_text = fcase(group_code == "HR", "HR",
                    group_code == "control", "C",
                    group_code == "R", "R",
                    default = "H99"),
  text_color = fcase(group_code == "HR", "#1A1A1A", default = "#FFFFFF"),
  retained_outline = FALSE)]
group_cells <- group_cells[, .(strain, col_id, fill_class, cell_text, text_color, retained_outline)]

erg11_cells <- erg11_parent[, .(strain, col_id = paste0("ERG11_", protein_change),
                                fill_class = "erg11_g484s", cell_text = "",
                                text_color = "#1A1A1A", retained_outline = FALSE)]

candidate_cells <- candidates[, .(strain, col_id = site_id, fill_class = event_class,
                                  cell_text = "", text_color = "#1A1A1A",
                                  retained_outline = retained)]

context_cells <- rbindlist(list(
  data.table(strain = LIMITED_COMPOSITE, col_id = "fig5_phenotype",
             fill_class = "phenotype_limited", cell_text = "low",
             text_color = "#1A1A1A", retained_outline = FALSE),
  data.table(strain = KB_SPECIFIC, col_id = "fig5_phenotype",
             fill_class = "phenotype_kb", cell_text = "KB",
             text_color = "#1A1A1A", retained_outline = FALSE)))

active_cells <- rbindlist(list(group_cells, erg11_cells, candidate_cells, context_cells),
                          use.names = TRUE, fill = TRUE)
plot_dt <- merge(background, active_cells, by = c("strain", "col_id"), all.x = TRUE, sort = FALSE)
plot_dt[, `:=`(strain_f = factor(strain, levels = rev(STR_LEVELS)),
               col_f = factor(col_id, levels = column_levels),
               block_f = factor(as.character(block), levels = BLOCK_LEVELS))]
plot_dt[, outline_status := fifelse(retained_outline == TRUE,
                                    "Retained in passage 03", "Not retained in passage 03")]

FILL_COLORS <- c(GROUP_COLORS, erg11_g484s = "#0072B2",
                 coding_snv = "#E69F00", noncoding_snv = "#CC79A7",
                 coding_indel = "#D55E00", noncoding_indel = "#56B4E9",
                 phenotype_limited = "#CAB2D6", phenotype_kb = "#B2DF8A",
                 phenotype_h99 = "#A6CEE3")
LEGEND_BREAKS <- c("erg11_g484s", "coding_snv", "noncoding_snv", "coding_indel",
                   "noncoding_indel", "phenotype_limited", "phenotype_kb", "phenotype_h99")
LEGEND_LABELS <- c("ERG11 p.G484S", "coding SNV", "non-coding SNV", "coding indel",
                   "non-coding indel", "low KB+PAP-AUC recovery",
                   "KB-specific persistence", "H99 recovered control")
x_label_map <- setNames(column_meta$label, column_meta$col_id)

pI <- ggplot(plot_dt, aes(col_f, strain_f)) +
  geom_tile(fill = "#F4F4F4", colour = "#FFFFFF", linewidth = 0.65, width = 0.88, height = 0.8) +
  geom_tile(data = plot_dt[!is.na(fill_class)], aes(fill = fill_class),
            colour = "#FFFFFF", linewidth = 0.65, width = 0.88, height = 0.8) +
  geom_tile(data = plot_dt[block == "Selected-clone-associated candidates" & !is.na(fill_class)],
            aes(colour = outline_status), fill = NA,
            linewidth = 0.85, width = 0.88, height = 0.8) +
  geom_text(data = plot_dt[!is.na(cell_text) & cell_text != "" & text_color == "#1A1A1A"],
            aes(label = cell_text), colour = "#1A1A1A", size = 2.15, fontface = "bold") +
  geom_text(data = plot_dt[!is.na(cell_text) & cell_text != "" & text_color == "#FFFFFF"],
            aes(label = cell_text), colour = "#FFFFFF", size = 2.15, fontface = "bold") +
  facet_grid(. ~ block_f, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = FILL_COLORS, breaks = LEGEND_BREAKS,
                    labels = LEGEND_LABELS, name = NULL, drop = TRUE) +
  scale_colour_manual(values = c(`Retained in passage 03` = "#1A1A1A",
                                 `Not retained in passage 03` = "#FFFFFF"),
                      breaks = "Retained in passage 03", name = NULL) +
  scale_x_discrete(labels = x_label_map) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE,
                             override.aes = list(colour = "#FFFFFF"), order = 1),
         colour = guide_legend(override.aes = list(fill = "#F4F4F4", colour = "#1A1A1A",
                                                   linewidth = 0.6), order = 2)) +
  theme_minimal(base_size = 9.5) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 6.5,
                                   lineheight = 0.92, colour = "#252525"),
        axis.text.y = element_text(size = 7.8, colour = "#252525", margin = margin(r = 2)),
        axis.ticks = element_blank(), panel.grid = element_blank(),
        panel.spacing.x = unit(4.5, "pt"),
        strip.background = element_rect(fill = "#EEF1F4", colour = NA),
        strip.text = element_text(face = "bold", colour = "#263746", size = 7.6,
                                  margin = margin(t = 4, r = 3, b = 4, l = 3)),
        legend.position = "bottom", legend.justification = "left",
        legend.box.just = "left", legend.key.width = unit(12, "pt"),
        legend.key.height = unit(9, "pt"), legend.text = element_text(size = 7.2),
        plot.margin = margin(t = 8, r = 10, b = 6, l = 7))
panelE <- pI


## -----------------------------------------------------------------------------
## 8. Calculated checks, audit table and PDF output
## -----------------------------------------------------------------------------
dir.create(DERIVED_DIR, recursive = TRUE, showWarnings = FALSE)
candidate_audit <- candidates[, .(
  strain, group = SV_WORKING$group_code[match(strain, SV_WORKING$strain)],
  chrom, chr, pos, ref, alt, type, qual, filter, gene, symbol, product,
  gene_cds, gene_context, event_class,
  vaf_parent = v_par, dp_parent = d_par, vaf_02 = v_02, dp_02 = d_02,
  vaf_01 = v_01, dp_01 = d_01, vaf_03 = v_03, dp_03 = d_03,
  passage_retained = retained,
  figure5_context = fcase(strain %chin% LIMITED_COMPOSITE, "limited KB+PAP-AUC recovery",
                          strain == KB_SPECIFIC, "KB-specific persistence",
                          default = "remaining isolate"))]
fwrite(candidate_audit, file.path(DERIVED_DIR, "fig5_candidate_variants.tsv"),
       sep = "\t", quote = FALSE, na = "NA")

message(
  "Figure 5: ", n_distinct(cnv$sample), " samples in the CNV heatmap; ",
  nrow(Cdat), " sample-chromosome pairs in Panel B; ",
  nrow(pick), " genes in Panel C (Chr1 ", nrow(chr1_pick),
  ", Chr4 ", nrow(chr4_pick), "); ",
  nrow(candidates), " small-variant candidates, ", sum(candidates$retained),
  " retained in ", length(retained_carriers), " isolates."
)
message(sprintf("Panel B: Pearson r = %.2f overall, %.2f excluding Chr2.", rC_all, rC_no2))
message(sprintf("Panel E: primary two-sided Fisher exact P = %.4f (KB-only sensitivity P = %.4f).",
                fisher_primary$p.value, fisher_kb_only))
message("Panel E: ERG11 p.G484S parent carriers: ", paste(g484s_carriers, collapse = ", "))
message(sprintf("Panel E: four-state callability %d/%d strain-gene combinations fully callable; minimum common CDS fraction %.3f.",
                callability_full, nrow(callability), callability_min))
message("Panel E: retention-threshold sensitivity")
print(retention_sensitivity)
print(association_table)

if (EXPORT_PANELS) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(OUTPUT_DIR)) stop("Cannot create output directory: ", OUTPUT_DIR)
  if (!capabilities("cairo")) stop("PDF output requires Cairo support in R.")
  panel_plots <- list(Fig5A = panelA, Fig5B = panelB, Fig5C = panelC, Fig5E = panelE)
  panel_sizes <- list(Fig5A = c(9.6, 3.3), Fig5B = c(5.0, 3.7),
                      Fig5C = c(5.2, 4.0), Fig5E = c(9.6, 6.65))
  iwalk(panel_plots, function(p, name) {
    size <- panel_sizes[[name]]
    ggsave(
      file.path(OUTPUT_DIR, paste0(name, ".pdf")), p,
      width = size[1], height = size[2], units = "in", device = cairo_pdf
    )
  })
  message("Figure 5 PDFs written to: ", normalizePath(OUTPUT_DIR))
}
