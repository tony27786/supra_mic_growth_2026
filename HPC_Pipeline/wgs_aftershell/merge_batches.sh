#!/usr/bin/env bash
# =============================================================================
# merge_batches.sh -- combine per-sample upstream outputs from multiple
# sequencing batches into ONE directory that hetero_qc.R and fig6.R read from.
#
# Merged content:
#   mosdepth   : *.mosdepth.summary.txt, *.regions.bed.gz(.csi),
#                *.thresholds.bed.gz(.csi), *.mosdepth.{global,region}.dist.txt
#   map QC     : qc_metrics_map/{flagstat,stats,idxstats}/*
#   het table  : all_het_snps.tsv  (headers de-duplicated, rows concatenated)
#
# Per-sample files are SYMLINKED (default), so both source batches stay intact
# and there is no extra disk use. Set MERGE_MODE=copy to copy instead.
# Re-running is safe (idempotent).
#
# FUTURE BATCH: append its dirs to the three arrays below and bump COMBINED_TAG.
# Keep COMBINED_TAG in sync with BATCH_TAG in analysis_config.R.
#
# Usage:
#   bash merge_batches.sh
#   COMBINED_TAG=260815 MERGE_MODE=copy bash merge_batches.sh
# =============================================================================
set -euo pipefail

WORK_ROOT="${WORK_ROOT:-$HOME/working}"
COMBINED_TAG="${1:-${COMBINED_TAG:?set COMBINED_TAG, and keep BATCH_TAG in wgs_analysis_config.R equal to it}}"
MERGE_MODE="${MERGE_MODE:-symlink}"          # symlink | copy
COMBINED_DIR="$WORK_ROOT/cn_combined_${COMBINED_TAG}"

# ---- source batches --------------------------------------------------------
# One entry per pipeline run to combine. If all samples went through the WGS
# pipeline in a single run, give one entry per array, using that run's RUN_TAG:
#
#   RUN_TAG=<your run tag>
#   MOSDEPTH_DIRS=( "$WORK_ROOT/cn_all_mosdepth_${RUN_TAG}/mosdepth" )
#   MAPQC_DIRS=(    "$WORK_ROOT/cn_replenish_map_${RUN_TAG}/qc_metrics_map" )
#   HET_TABLES=(    "$WORK_ROOT/cn_re_per_sample_vcf_${RUN_TAG}/het_extracted/all_het_snps.tsv" )
#
# The original study merged two runs, because the Parent bulks and the clone
# states were sequenced in separate batches. Merging is still required for a
# single run: wgs_aftershell.R derives every path from COMBINED_DIR, and the
# heterozygosity robust z score is computed across the whole merged cohort.
RUN_TAG="${RUN_TAG:?set RUN_TAG (the tag frozen by wgs_pipeline/submit_all.sh)}"
MOSDEPTH_DIRS=(
  "$WORK_ROOT/cn_all_mosdepth_${RUN_TAG}/mosdepth"
)
MAPQC_DIRS=(
  "$WORK_ROOT/cn_replenish_map_${RUN_TAG}/qc_metrics_map"
)
HET_TABLES=(
  "$WORK_ROOT/cn_re_per_sample_vcf_${RUN_TAG}/het_extracted/all_het_snps.tsv"
)
# ----------------------------------------------------------------------------

echo ">>> COMBINED_DIR = $COMBINED_DIR   (mode: $MERGE_MODE)"
mkdir -p "$COMBINED_DIR/mosdepth" \
         "$COMBINED_DIR/qc_metrics_map/flagstat" \
         "$COMBINED_DIR/qc_metrics_map/stats" \
         "$COMBINED_DIR/qc_metrics_map/idxstats" \
         "$COMBINED_DIR/het_extracted"

COLLISIONS=0

place_one() {                 # place_one <dest_dir> <src_file>
  local dest_dir="$1" src="$2"
  [[ -e "$src" ]] || return 0
  local base srcabs tgt cur
  base="$(basename "$src")"
  srcabs="$(readlink -f "$src")"
  tgt="$dest_dir/$base"
  if [[ -e "$tgt" || -L "$tgt" ]]; then
    cur="$(readlink -f "$tgt" 2>/dev/null || echo "")"
    if [[ "$cur" != "$srcabs" ]]; then
      echo "[COLLISION] $base"
      echo "    have: $cur"
      echo "    new : $srcabs"
      COLLISIONS=$((COLLISIONS + 1))
    fi
    return 0
  fi
  if [[ "$MERGE_MODE" == "copy" ]]; then
    cp -pv "$srcabs" "$tgt"
  else
    ln -s "$srcabs" "$tgt"
  fi
}

# ---- mosdepth --------------------------------------------------------------
echo ">>> mosdepth ..."
mos_n=0
for d in "${MOSDEPTH_DIRS[@]}"; do
  [[ -d "$d" ]] || { echo "[WARN] missing mosdepth dir: $d"; continue; }
  for f in "$d"/*.mosdepth.summary.txt \
           "$d"/*.regions.bed.gz "$d"/*.regions.bed.gz.csi \
           "$d"/*.thresholds.bed.gz "$d"/*.thresholds.bed.gz.csi \
           "$d"/*.mosdepth.global.dist.txt "$d"/*.mosdepth.region.dist.txt; do
    [[ -e "$f" ]] || continue
    place_one "$COMBINED_DIR/mosdepth" "$f"
    mos_n=$((mos_n + 1))
  done
done
echo "    placed $mos_n files; samples (by summary.txt): $(ls -1 "$COMBINED_DIR"/mosdepth/*.mosdepth.summary.txt 2>/dev/null | wc -l)"

# ---- map QC ----------------------------------------------------------------
echo ">>> mapping QC ..."
for d in "${MAPQC_DIRS[@]}"; do
  [[ -d "$d" ]] || { echo "[WARN] missing qc_metrics_map dir: $d"; continue; }
  for sub in flagstat stats idxstats; do
    [[ -d "$d/$sub" ]] || continue
    for f in "$d/$sub"/*; do
      [[ -e "$f" ]] || continue
      place_one "$COMBINED_DIR/qc_metrics_map/$sub" "$f"
    done
  done
done
echo "    flagstat=$(ls -1 "$COMBINED_DIR"/qc_metrics_map/flagstat 2>/dev/null | wc -l) stats=$(ls -1 "$COMBINED_DIR"/qc_metrics_map/stats 2>/dev/null | wc -l) idxstats=$(ls -1 "$COMBINED_DIR"/qc_metrics_map/idxstats 2>/dev/null | wc -l)"

# ---- het table (concatenate, single header) --------------------------------
echo ">>> het table ..."
out_het="$COMBINED_DIR/het_extracted/all_het_snps.tsv"
: > "$out_het"
ref_header=""
first=1
for t in "${HET_TABLES[@]}"; do
  [[ -f "$t" ]] || { echo "[WARN] missing het table: $t"; continue; }
  h="$(head -1 "$t")"
  if (( first )); then
    ref_header="$h"
    cat "$t" >> "$out_het"
    first=0
  else
    if [[ "$h" != "$ref_header" ]]; then
      echo "[WARN] het header differs from first table -- check schema:"
      echo "    first: $ref_header"
      echo "    this : $h    ($t)"
    fi
    tail -n +2 "$t" >> "$out_het"
  fi
done
if [[ -s "$out_het" ]]; then
  n_rows=$(( $(wc -l < "$out_het") - 1 ))
  n_samp=$(tail -n +2 "$out_het" | cut -f1 | sort -u | wc -l)
  echo "    het rows: $n_rows across $n_samp samples -> $out_het"
else
  echo "[ERROR] no het rows written; check HET_TABLES paths." >&2
  exit 1
fi

# ---- summary ---------------------------------------------------------------
echo
if (( COLLISIONS > 0 )); then
  echo "=============================================================="
  echo "[WARN] $COLLISIONS filename collision(s): same basename from different sources."
  echo "       Usually means a sample name is duplicated across batches."
  echo "       The first occurrence is kept; resolve before trusting cohort stats."
  echo "=============================================================="
fi
echo ">>> Merge complete."
echo ">>> Set analysis_config.R:  BATCH_TAG <- \"${COMBINED_TAG}\"   (COMBINED_DIR = $COMBINED_DIR)"
echo ">>> Then run hetero_qc.R, then fig6.R."
