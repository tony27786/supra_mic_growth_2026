# =============================================================================
# analysis_config.R -- shared batch paths for hetero_qc.R and fig6.R.
#
# source("./wgs_analysis_config.R") runs at the top of the other R scripts.
# Keep all of them (this, hetero_cnv.R, hetero_qc.R, wgs_aftershell.R,
# merge_batches.sh) in one folder and make that folder the R working directory.
#
# To analyze a different merged batch, edit BATCH_TAG only.
# Keep BATCH_TAG in sync with COMBINED_TAG in merge_batches.sh.
# =============================================================================

# ---- EDIT PER BATCH --------------------------------------------------------
BATCH_TAG <- Sys.getenv("BATCH_TAG", "260623")    # must equal COMBINED_TAG in merge_batches.sh
WORK_ROOT <- Sys.getenv("WORK_ROOT", "~/working") # scratch root; same value the SLURM stages used
# ----------------------------------------------------------------------------

COMBINED_DIR <- file.path(WORK_ROOT, paste0("cn_combined_", BATCH_TAG))

# Derived inputs consumed by the two R scripts (do not edit) -----------------
MAP_DIR          <- COMBINED_DIR                                 # holds qc_metrics_map/
MOS_DIR          <- file.path(COMBINED_DIR, "mosdepth")
DEPTH_DIR        <- COMBINED_DIR                                 # Step 4 writes its long tables here; Step 5 reads them
HET_FILE         <- file.path(COMBINED_DIR, "het_extracted", "all_het_snps.tsv")
VCF_FILTERED_DIR <- file.path(COMBINED_DIR, "vcf", "filtered")   # used only for Step 5 out-dir derivation

# Strain grouping (fig6.R): col1 = id (generic base, e.g. 101), col2 = group (HR/control/R).
GROUP_CSV <- "./group.csv"

# Track A reference-free ploidy summary, written by trackA_kmer_ploidy.slurm.
# Point this at that job's output directory; the tag is the SV/kmer run tag,
# which is deliberately independent of BATCH_TAG.
KMER_SUMMARY_TSV <- Sys.getenv(
  "KMER_SUMMARY_TSV",
  file.path(WORK_ROOT, "trackA_kmer", "trackA_ploidy_summary.tsv")
)

# Single-colony (monoclonal) detection from sample name. Current naming:
#   parent = Cn<strain>           (Cn160, Cn201, Cn169b)
#   clone  = Cn<strain><batch>    strain = 3 chars (3-digit OR H99), batch 01/02/03
#            (Cn16001, CnH9902) -- ALL single-colony-derived => monoclonal.
# This regex MUST match CLONE_RE in fig6_aftershell.R. (Old CnNNNm scheme retired.)
is_mono <- function(s) stringr::str_detect(s, "^Cn(\\d{3}|H99)0[12]$")
# Strictest no-selection floor instead = xxx01 only (parent single CFU, no drug):
#   is_mono <- function(s) stringr::str_detect(s, "^Cn(\\d{3}|H99)01$")

# H99 is display-only: never a het/depth baseline anchor (project convention).
is_h99 <- function(s) stringr::str_detect(toupper(s), "H99")

# Het mixture baseline source, applied to BOTH the genome-level (het_robust_z)
# and per-chromosome (chrom_het_baseline) references in fig6.R.
#   "mono" = single-colony clones (xxx01/02/03) only -- clean monoclonal reference;
#   "all"  = every sample.  (H99 + its clones are always excluded, see is_h99.)
HET_BASELINE <- "mono"

# Rows of `df` (must have a `sample` column) to use as the het baseline.
het_baseline_rows <- function(df) {
  df <- dplyr::filter(df, !is_h99(sample))          # H99 never in a baseline
  if (identical(HET_BASELINE, "mono")) {
    sub <- dplyr::filter(df, is_mono(sample))
    if (dplyr::n_distinct(sub$sample) >= 2) return(sub)
    warning("HET_BASELINE='mono' but <2 mono (m) samples; baseline falls back to all (non-H99).")
  }
  df
}

message("[analysis_config] BATCH_TAG = ", BATCH_TAG)
message("[analysis_config] COMBINED_DIR = ", COMBINED_DIR)