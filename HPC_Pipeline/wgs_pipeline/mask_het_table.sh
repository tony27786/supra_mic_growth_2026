#!/bin/bash
# =============================================================================
# mask_het_table.sh  --  apply the mappability mask to an EXISTING het table.
#
# Works on ANY all_het_snps.tsv (old-pipeline cohort, new replenish batch, or
# the MERGED table) because the mask is reference-coordinate-based and pipeline-
# agnostic. This is the right place to mask the CURRENT analysis: your fig6.r
# input is a merge of two pipelines, so masking the merged table applies one
# uniform ruler to old + new at once.
#
# It does two things:
#   1) writes masked_all_het_snps.tsv  (same schema, only het in unique regions)
#   2) writes per_sample_mask_summary.tsv  (n_het_all vs n_het_masked, % kept,
#      mean AF before/after) -> tells you immediately whether the ~5000 cohort
#      baseline collapses toward 0 (artifact) and how Cn160 behaves.
#
# Run inside the bwasam_snake env (needs bedtools on PATH):
#   bash mask_het_table.sh <all_het_snps.tsv> <mask.bed> [out_dir]
# Defaults: mask = $MASK_BED if 00_config.sh is sourceable; out_dir = ./mask_out
#
# Input schema (header required, as written by 04_extract_het.slurm):
#   sample chrom pos ref alt qual filter_st gt dp ad
#   ad = "ref_count,alt_count"  (used to compute alt AF = alt/(ref+alt))
# =============================================================================
set -eo pipefail

IN="${1:-}"
MASK="${2:-${MASK_BED:-}}"
OUTDIR="${3:-./mask_out}"

usage() { echo "usage: bash mask_het_table.sh <all_het_snps.tsv> <mask.bed> [out_dir]" >&2; exit 1; }
[[ -n "$IN"   && -f "$IN"   ]] || { echo "[ERROR] het table not found: '${IN:-<unset>}'" >&2; usage; }
[[ -n "$MASK" && -f "$MASK" ]] || { echo "[ERROR] mask BED not found: '${MASK:-<unset>}'" >&2; usage; }
command -v bedtools >/dev/null 2>&1 || { echo "[ERROR] bedtools not on PATH (activate bwasam_snake)" >&2; exit 1; }

mkdir -p "$OUTDIR"
OUT="${OUTDIR}/masked_all_het_snps.tsv"
SUMMARY="${OUTDIR}/per_sample_mask_summary.tsv"
HEADER="sample	chrom	pos	ref	alt	qual	filter_st	gt	dp	ad"

echo "[$(date '+%F %T')] >>> input : $IN"
echo "[$(date '+%F %T')] >>> mask  : $MASK"
echo "[$(date '+%F %T')] >>> out   : $OUTDIR"

# ---- 1) intersect het positions with the unique-region mask -----------------
# het table is 1-based POS -> BED (pos-1, pos); carry all non-coordinate fields
# through the intersect, then reformat back to the original column order.
# bedtools loads the (small) mask into memory, so -a need not be pre-sorted.
{
  echo -e "$HEADER"
  tail -n +2 "$IN" \
  | awk 'BEGIN{OFS="\t"} NF>=10 {print $2, $3-1, $3, $1, $4, $5, $6, $7, $8, $9, $10}' \
  | bedtools intersect -a - -b "$MASK" -u \
  | awk 'BEGIN{OFS="\t"} {print $4, $1, $3, $5, $6, $7, $8, $9, $10, $11}'
} > "$OUT"

# ---- 2) per-sample before/after summary ------------------------------------
# Reads IN then OUT (arg order); buckets rows by FILENAME. AF = alt/(ref+alt).
awk -F'\t' -v INF="$IN" -v OUTF="$OUT" '
  function afc(s,   a){ split(s, a, ","); if (a[1]=="" || a[2]=="" || (a[1]+a[2])<=0) return -1; return a[2]/(a[1]+a[2]) }
  $1=="sample" { next }                       # skip header rows in both files
  {
    s=$1; v=afc($10)
    if (FILENAME==INF) { na[s]++; if (v>=0){ sa[s]+=v; ca[s]++ } }
    else               { nm[s]++; if (v>=0){ sm[s]+=v; cm[s]++ } }
  }
  END{
    for (s in na){
      m  = (s in nm) ? nm[s] : 0
      pr = (na[s]>0) ? 100*m/na[s] : 0
      aa = (ca[s]>0) ? sa[s]/ca[s] : -1
      am = (cm[s]>0) ? sm[s]/cm[s] : -1
      printf "%s\t%d\t%d\t%.1f\t%s\t%s\n", s, na[s], m, pr,
             (aa>=0 ? sprintf("%.3f",aa) : "NA"),
             (am>=0 ? sprintf("%.3f",am) : "NA")
    }
  }
' "$IN" "$OUT" \
| sort -t$'\t' -k1,1V \
| awk 'BEGIN{print "sample\tn_het_all\tn_het_masked\tpct_retained\tmeanAF_all\tmeanAF_masked"} {print}' \
> "$SUMMARY"

NALL=$(($(wc -l < "$IN") - 1))
NMSK=$(($(wc -l < "$OUT") - 1))
echo "=============================================================="
echo "[$(date '+%F %T')] >>> het SNPs total : $NALL  -> masked : $NMSK"
echo "[$(date '+%F %T')] >>> masked table   : $OUT"
echo "[$(date '+%F %T')] >>> summary        : $SUMMARY"
echo "--------------------------------------------------------------"
echo ">>> per-sample (n_het_all -> n_het_masked, % kept):"
column -t -s$'\t' "$SUMMARY" 2>/dev/null || cat "$SUMMARY"
echo "=============================================================="
echo ">>> Read: cohort rows should drop hard (artifact baseline gone);"
echo ">>>       Cn160/Cn160m may stay high (genome-wide divergence, not repeats)."
exit 0
