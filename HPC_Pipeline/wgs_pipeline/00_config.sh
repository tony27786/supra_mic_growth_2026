#!/usr/bin/env bash
# =============================================================================
# 00_config.sh  --  Shared configuration for the Fig6 WGS upstream pipeline.
#
# This file is SOURCED by every stage script. It only DEFINES variables and
# helper functions; it must not enable strict mode or have side effects.
#
# Stages keyed on a single RUN_TAG so all input/output directories line up:
#   00_fastp      -> $WORK_ROOT/cn_replenish_fastpout_${RUN_TAG}  (== CLEAN_DIR)
#   01_map        -> $WORK_ROOT/cn_replenish_map_${RUN_TAG}
#   02_mosdepth   -> $WORK_ROOT/cn_replenish_mosdepth_${RUN_TAG}
#   03_vcf        -> $WORK_ROOT/cn_re_per_sample_vcf_${RUN_TAG}/vcf[/filtered]
#   04_extract_het-> $WORK_ROOT/cn_re_per_sample_vcf_${RUN_TAG}/het_extracted/all_het_snps.tsv
#
# RUN_TAG is frozen ONCE by submit_all.sh (written to run_tag.txt and exported),
# so the chain stays consistent even if jobs queue across day boundaries.
#
# Mappability mask (NEW): built ONCE on REF by 00_make_mappability_mask.slurm,
# reused by stage 04 (-T) and by mask_het_table.sh. See the "Mappability mask"
# block below. Leave MASK_BED empty to disable masking everywhere.
# =============================================================================

# =============================================================================
# SITE SETTINGS -- set these before submitting anything.
#
# Every value below is site specific. Export them in your shell, or edit the
# defaults here. Nothing else in the pipeline hard-codes a path.
#
#   export WORK_ROOT=/scratch/$USER/hetero      # all run directories go here
#   export RAW_DIR=/data/hetero/fastq           # {sample}_R1.fq.gz / _R2.fq.gz
#   export REF_FA=/refs/GCF_000149245.1_CNA3_genomic.fna
#   export SLURM_PARTITION=<your partition>
#   export CONDA_LAUNCHER=$(command -v conda)
#
# The published run used a SLURM cluster with lustre scratch; the directory
# layout below is preserved so the stage-to-stage wiring is unchanged.
# =============================================================================

# ---- Identity / SLURM ----
# Optional accounting string passed to `sbatch --comment`. Left empty by
# default; stage scripts omit the flag entirely when it is empty.
SBATCH_COMMENT="${SBATCH_COMMENT:-}"
# SLURM partition. There is no portable default -- set it for your site.
SLURM_PARTITION="${SLURM_PARTITION:-}"

# ---- Conda ----
# Path to the conda launcher. Defaults to whatever is on PATH.
CONDA_LAUNCHER="${CONDA_LAUNCHER:-$(command -v conda || echo conda)}"
CONDA_ENV_FASTP="${CONDA_ENV_FASTP:-fastp}"             # fastp
CONDA_ENV_BWASAM="${CONDA_ENV_BWASAM:-bwasam_snake}"    # bwa / samtools / bcftools / snakemake / genmap / bedtools
CONDA_ENV_MOSDEPTH="${CONDA_ENV_MOSDEPTH:-mosdepth}"    # mosdepth

# ---- Roots ----------------------------------------------------------------
# WORK_ROOT: where all run directories live (scratch space on this cluster).
WORK_ROOT="${WORK_ROOT:?set WORK_ROOT to a writable scratch directory}"

# Reference FASTA (bwa/samtools index is built next to it on first run by 01_map).
# NCBI RefSeq GCF_000149245.1 -- C. neoformans var. grubii H99, assembly CNA3.
REF_DIR="${REF_DIR:-$WORK_ROOT/h99ref}"
REF_FA="${REF_FA:-$REF_DIR/GCF_000149245.1_CNA3_genomic.fna}"

# Raw sequencing reads for stage 00 (fastp). Flat layout, one pair per sample:
#   {sample}_R1.fq.gz  and  {sample}_R2.fq.gz
# Sample names must match the `sample` column of ../manifest/wgs_samples.tsv.
RAW_DIR="${RAW_DIR:?set RAW_DIR to the directory holding the FASTQ pairs}"

# ---- Mappability mask (built once on REF; see 00_make_mappability_mask.slurm)
# Uniquely-mappable, single-copy regions of REF. Counting het only inside this
# BED strips the paralogous / mis-mapping het that otherwise dominates the count
# (the ~5000/sample H99 baseline). Same .bed feeds stage 04 (`bcftools query
# -T`) and the standalone mask_het_table.sh (`bedtools intersect -b`).
#   - Set MASK_BED to "" to disable masking (stage 04 reverts to original behavior).
#   - It is intentionally NOT keyed on RUN_TAG: one mask per reference, reused.
MASK_K="${MASK_K:-150}"                              # GenMap k-mer length (match read length, ~151 here)
MASK_E="${MASK_E:-1}"                                # mismatches allowed when judging uniqueness
MASK_MIN_MAPPABILITY="${MASK_MIN_MAPPABILITY:-1}"    # keep positions with mappability >= this (1 = strictly unique)
GENMAP_DIR="${GENMAP_DIR:-$(dirname "$REF_FA")/genmap}"            # index + bedgraph live next to the reference
GENMAP_INDEX="${GENMAP_INDEX:-$GENMAP_DIR/index}"
GENMAP_PREFIX="${GENMAP_PREFIX:-$GENMAP_DIR/h99_K${MASK_K}_E${MASK_E}}"
MASK_BED="${MASK_BED:-$GENMAP_DIR/h99_mappable_K${MASK_K}_E${MASK_E}.bed}"

# ---- Run tag (frozen once by submit_all.sh) -------------------------------
# Priority: exported RUN_TAG env  >  run_tag.txt in pipeline dir  >  today
PIPELINE_DIR="${PIPELINE_DIR:-$(pwd)}"
if [[ -z "${RUN_TAG:-}" ]]; then
  if [[ -f "${PIPELINE_DIR}/run_tag.txt" ]]; then
    RUN_TAG="$(tr -d '[:space:]' < "${PIPELINE_DIR}/run_tag.txt")"
  else
    RUN_TAG="$(date +%y%m%d)"
  fi
fi

# ---- Derived run directories (all keyed on RUN_TAG) -----------------------
# fastp output == the cleaned-reads dir that stage 01 (mapping) reads.
CLEAN_DIR="${CLEAN_DIR:-$WORK_ROOT/cn_all_fastpout_${RUN_TAG}}"
MAP_OUT="$WORK_ROOT/cn_all_map_${RUN_TAG}"
MOS_OUT="$WORK_ROOT/cn_all_mosdepth_${RUN_TAG}"
VCF_OUT="$WORK_ROOT/cn_per_sample_vcf_${RUN_TAG}"
VCF_FILTERED_DIR="$VCF_OUT/vcf/filtered"            # consumed by stage 04
HET_OUT="$VCF_OUT/het_extracted"

# ---- Stage 00: fastp parallelism ------------------------------------------
FASTP_MAX_JOBS="${FASTP_MAX_JOBS:-5}"

# ---- Stage 01: mapping parallelism ----------------------------------------
MAP_MAX_JOBS="${MAP_MAX_JOBS:-6}"
MAP_SORT_THREADS="${MAP_SORT_THREADS:-2}"
MAP_INDEX_THREADS="${MAP_INDEX_THREADS:-2}"

# ---- Stage 02: mosdepth ----------------------------------------------------
MOS_MAX_JOBS="${MOS_MAX_JOBS:-16}"
MOS_WINDOW_SIZE="${MOS_WINDOW_SIZE:-1000}"
MOS_FAST_MODE="${MOS_FAST_MODE:-1}"
MOS_NO_PER_BASE="${MOS_NO_PER_BASE:-1}"
# Coverage thresholds: default ON (matches the original invocation). They are
# NOT used by the CNV/aneuploidy analysis (that reads mosdepth.summary.txt
# per-chrom depth + regions.bed.gz bins) and are only approximate under
# --fast-mode, so disabling is harmless. To disable: export MOS_THRESHOLDS=""
MOS_THRESHOLDS="${MOS_THRESHOLDS-1,5,10,20,50,100}"

# ---- Stage 03: variant calling (mirrored into config.yaml) -----------------
VC_MIN_MAPQ="${VC_MIN_MAPQ:-20}"
VC_MIN_BQ="${VC_MIN_BQ:-20}"
VC_PLOIDY="${VC_PLOIDY:-2}"          # 2 is REQUIRED: ploidy 1 emits only
                                      # homozygous GTs, which would discard the
                                      # het sites / AF that disomy + mixture
                                      # detection depend on.
VC_CALL_THREADS="${VC_CALL_THREADS:-1}"
VC_QUAL_MIN="${VC_QUAL_MIN:-30}"
VC_DEPTH_MIN="${VC_DEPTH_MIN:-10}"
VC_MIN_CONTIG_LEN="${VC_MIN_CONTIG_LEN:-1000}"
# bcftools mpileup -d (max reads/site). 250 == bcftools default (unchanged).
# For ~200x WGS you may raise it so DP/AD reflect true depth, e.g.
#   export VC_MAX_DEPTH=2000     (or 0 for unlimited)
VC_MAX_DEPTH="${VC_MAX_DEPTH:-250}"
VC_SNAKE_JOBS="${VC_SNAKE_JOBS:-30}"

# ---- Stage 04: het extraction ---------------------------------------------
HET_MAX_JOBS="${HET_MAX_JOBS:-16}"

# ============================================================================
# Stage 05: 4-state small-variant forensics (05_smallvar_4state.reviewed_v2)
# ============================================================================
# Per strain we co-pileup 4 states:
#   parent Cn{S}   02 Cn{S}02   01 Cn{S}01   03 Cn{S}03
# The stage takes a parent BAM dir AND a clone BAM dir, because in the original
# study the Parent bulks and the clone states were mapped in separate runs. If
# you mapped all 80 samples in one run, point both at that run's output:
#   export SV_PARENT_BAM_DIR=$WORK_ROOT/cn_all_map_${RUN_TAG}
#   export SV_CLONE_BAM_DIR=$SV_PARENT_BAM_DIR
# Expected content: 19 strains x 4 states = 76 BAMs, plus CnH99 x 4 = 80.
# NB a parent dir may also hold Cn169b (a second Cn169 sample); stage 05 uses Cn169.
SV_PARENT_BAM_DIR="${SV_PARENT_BAM_DIR:?set SV_PARENT_BAM_DIR to the directory holding the Parent BAMs}"
SV_CLONE_BAM_DIR="${SV_CLONE_BAM_DIR:-$SV_PARENT_BAM_DIR}"

# Frozen 25-gene panel, union-CDS, BED4 with col4 = CNAG (25/25 genes, 271
# intervals, 68,842 bp). Built locally by the repo; upload it next to the scripts:
#   scp hetero_fig6_wf/fig6_panel25_cds.bed <cluster>:$PIPELINE_DIR/
# Required: without it the per-gene 4-state callability track cannot be built,
# and a "no variant in the target panel" conclusion is not supportable.
SV_TARGET_BED="${SV_TARGET_BED:-$PIPELINE_DIR/fig6_panel25_cds.bed}"

# FRESH tag per forensics run -- deliberately NOT the pipeline RUN_TAG, so the
# outputs never mix with an earlier run. Stage 05 refuses to overwrite an
# existing dir unless FORCE=1, so re-running the same day needs an explicit tag.
SV_TAG="${SV_TAG:-smv$(date +%y%m%d)}"

# Calling / callability parameters.
SV_MAX_DEPTH="${SV_MAX_DEPTH:-4000}"     # mpileup -d; 250 (pipeline default) is too low for ~200x
SV_DP_CALLABLE="${SV_DP_CALLABLE:-10}"   # min depth for a base to count as callable
SV_QUAL_SOFT="${SV_QUAL_SOFT:-20}"       # QUAL below this -> soft FILTER=LowQual
SV_MAX_JOBS="${SV_MAX_JOBS:-20}"         # strains in parallel (FIFO semaphore)
SV_INCLUDE_H99="${SV_INCLUDE_H99:-1}"    # 1 = also process CnH99 as a control
SV_ALLOW_NO_MASK="${SV_ALLOW_NO_MASK:-0}"   # 1 = allow running without MASK_BED (not advised)
# Optional consequence annotation. bcftools csq needs a VALIDATED Ensembl-format
# GFF3 -- the RefSeq .gtf will NOT work. Leave empty to annotate locally in step 2.
SV_GFF3="${SV_GFF3:-}"

# ---- Stage 05: 4-state small-variant forensics (05_smallvar_4state.*) ------
# Joint per-strain calling across parent / 02 / 01 / 03 for the Fig6 "Panel I"
# analysis. Everything here is overridable via `sbatch --export=ALL,VAR=...`.
#
# BAM_DIR: MUST contain all four states of every strain:
#     Cn{S}.bam  Cn{S}01.bam  Cn{S}02.bam  Cn{S}03.bam   (+ .bai)
# The initial and replenish batches map into DIFFERENT dirs, so if your 80 BAMs
# are split, symlink them into one dir and point BAM_DIR there. The stage-05
# preflight writes input_manifest.tsv listing exists/quickcheck/index/@RG-SM per
# state, and refuses to run any strain that is not complete -- check that file
# first if strains get excluded.
BAM_DIR="${BAM_DIR:-$MAP_OUT}"

# SMALLVAR_TAG: keyed like RUN_TAG (env > smallvar_tag.txt > date), but it must
# be FRESH per analysis run: stage 05 refuses to write into an existing output
# dir unless FORCE=1, so frozen results are never silently overwritten/mixed.
if [[ -z "${SMALLVAR_TAG:-}" ]]; then
  if [[ -f "${PIPELINE_DIR}/smallvar_tag.txt" ]]; then
    SMALLVAR_TAG="$(tr -d '[:space:]' < "${PIPELINE_DIR}/smallvar_tag.txt")"
  else
    SMALLVAR_TAG="smv$(date +%y%m%d)"
  fi
fi
SMALLVAR_OUT="$WORK_ROOT/fig6_smallvar_${SMALLVAR_TAG}"

# Frozen 25-gene panel (generated locally, copy both files next to this config).
# TARGET_BED is BED4 with column 4 = CNAG id; it gates the callability track that
# a target-negative conclusion depends on. fig6_panel25.tsv is the human-readable
# freeze (symbol/CNAG/tier/alias) -- not read by the script, keep for the methods.
TARGET_BED="${TARGET_BED:-$PIPELINE_DIR/fig6_panel25_cds.bed}"
TARGET_TSV="${TARGET_TSV:-$PIPELINE_DIR/fig6_panel25.tsv}"

# Consequence annotation: bcftools csq needs a VALIDATED Ensembl-format GFF3.
# The RefSeq .gtf will NOT work. Leave empty to skip and annotate downstream in R.
GFF3="${GFF3:-}"

SV_MAX_DEPTH="${SV_MAX_DEPTH:-4000}"        # mpileup -d; 250 default is too low for ~200x
SV_MIN_MAPQ="${SV_MIN_MAPQ:-$VC_MIN_MAPQ}"  # 20
SV_MIN_BQ="${SV_MIN_BQ:-$VC_MIN_BQ}"        # 20
SV_DP_CALLABLE="${SV_DP_CALLABLE:-10}"      # min DP for a base to count as callable
SV_QUAL_SOFT="${SV_QUAL_SOFT:-20}"          # QUAL below this -> soft FILTER=LowQual
SV_MAX_JOBS="${SV_MAX_JOBS:-16}"            # strains in parallel (FIFO semaphore)
INCLUDE_H99="${INCLUDE_H99:-0}"             # 1 = also process H99's four states
ALLOW_NO_MASK="${ALLOW_NO_MASK:-0}"         # 1 = allow running without MASK_BED (not advised)
FORCE="${FORCE:-1}"                         # 1 = allow writing into an existing SMALLVAR_OUT

# ============================================================================
# Helpers
# ============================================================================

# Activate a conda env while preserving the caller's 'set -u' state.
activate_conda() {
  local env_name="$1"
  local had_u=0
  case "$-" in *u*) had_u=1 ;; esac
  local conda_base
  conda_base="$("$CONDA_LAUNCHER" info --base)"
  # shellcheck source=/dev/null
  source "${conda_base}/etc/profile.d/conda.sh"
  set +u
  conda activate "$env_name"
  (( had_u )) && set -u
  echo "[$(date '+%F %T')] conda env = ${env_name} (CONDA_PREFIX=${CONDA_PREFIX:-NA})"
}

# Apply the standard strict mode + single-thread BLAS used by all stages.
apply_strict_mode() {
  set -euo pipefail
  trap 'code=$?; echo -e "\n[ERROR] Exit $code at line $LINENO" >&2' ERR
  IFS=$'\n\t'
  shopt -s nullglob
  export PYTHONUNBUFFERED=1
  export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
}

# Print the resolved configuration banner.
print_config_banner() {
  echo "[$(date '+%F %T')] >>> RUN_TAG     = $RUN_TAG"
  echo "[$(date '+%F %T')] >>> WORK_ROOT   = $WORK_ROOT"
  echo "[$(date '+%F %T')] >>> REF_FA      = $REF_FA"
  echo "[$(date '+%F %T')] >>> RAW_DIR     = $RAW_DIR"
  echo "[$(date '+%F %T')] >>> CLEAN_DIR   = $CLEAN_DIR"
  echo "[$(date '+%F %T')] >>> MAP_OUT     = $MAP_OUT"
  echo "[$(date '+%F %T')] >>> MOS_OUT     = $MOS_OUT"
  echo "[$(date '+%F %T')] >>> VCF_OUT     = $VCF_OUT"
  echo "[$(date '+%F %T')] >>> HET_OUT     = $HET_OUT"
  echo "[$(date '+%F %T')] >>> MASK_BED    = ${MASK_BED:-<none>}"
}

# Stage-05 banner (small-variant forensics); called by 05_smallvar_*.slurm.
print_smallvar_banner() {
  echo "[$(date '+%F %T')] >>> SV_TAG            = ${SV_TAG:-<unset>}"
  echo "[$(date '+%F %T')] >>> SV_PARENT_BAM_DIR = ${SV_PARENT_BAM_DIR:-<unset>}"
  echo "[$(date '+%F %T')] >>> SV_CLONE_BAM_DIR  = ${SV_CLONE_BAM_DIR:-<unset>}"
  echo "[$(date '+%F %T')] >>> SV_TARGET_BED     = ${SV_TARGET_BED:-<unset>}"
  echo "[$(date '+%F %T')] >>> SV_GFF3           = ${SV_GFF3:-<none, annotate in step 2>}"
  echo "[$(date '+%F %T')] >>> depth<=${SV_MAX_DEPTH} callable_DP>=${SV_DP_CALLABLE} jobs=${SV_MAX_JOBS} H99=${SV_INCLUDE_H99}"
}
