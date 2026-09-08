#!/bin/bash
#SBATCH -J cn04_het
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -o 04_extract_het.out
# =============================================================================
# Stage 04: extract biallelic heterozygous SNPs from filtered per-sample VCFs.
# Reads : $VCF_OUT/vcf/filtered/{sample}.filtered.vcf.gz
# Writes: $HET_OUT/all_het_snps.tsv   (one combined table; the fig6.r input)
#
# NEW: if MASK_BED (00_config.sh) is set and exists, het is counted ONLY inside
# the uniquely-mappable regions (bcftools query -T $MASK_BED), stripping the
# paralogous / mis-mapping het baseline. Set MASK_BED="" to disable and exactly
# reproduce the original genome-wide behavior. Re-running this stage alone is
# cheap (filtered VCFs already exist; nothing is re-called).
#
# NOTE: this masks ONLY the het produced by THIS pipeline. The current fig6.r
# input is a MERGE of this batch with the old-pipeline cohort, so to put one
# uniform ruler on old + new, run mask_het_table.sh on the MERGED table instead
# (or in addition -- intersecting twice with the same BED is idempotent).
# =============================================================================

PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${PIPELINE_DIR}/00_config.sh"
# NOTE: '-u' intentionally NOT enabled here (mirrors original) but trap/pipefail are.
set -eo pipefail
trap 'code=$?; echo -e "\n[ERROR] Exit $code at line $LINENO" >&2' ERR
IFS=$'\n\t'
shopt -s nullglob

VCF_DIR="$VCF_FILTERED_DIR"
HET_LOG_DIR="$HET_OUT/logs"
FAIL_DIR="$HET_OUT/.failed"
mkdir -p "$HET_OUT" "$HET_LOG_DIR" "$FAIL_DIR"
rm -f "$FAIL_DIR"/* 2>/dev/null || true

cd "$HET_OUT"
echo "Working Directory = $(pwd)"
activate_conda "$CONDA_ENV_BWASAM"
print_config_banner
echo "[$(date '+%F %T')] >>> VCF_DIR = $VCF_DIR"

# ---- Mask gate: restrict het to uniquely-mappable regions if configured -----
QUERY_TARGETS=()
if [[ -n "${MASK_BED:-}" && -f "${MASK_BED}" ]]; then
  QUERY_TARGETS=(-T "${MASK_BED}")
  echo "[$(date '+%F %T')] >>> MASK ON  : counting het only within $MASK_BED"
else
  echo "[$(date '+%F %T')] >>> MASK OFF : MASK_BED unset/missing -> genome-wide het (original behavior)"
fi

# ---- Het quality gate (MixInfect2-style, DIVISION-FREE) ---------------------
# A credible het needs the minor allele backed by enough reads AND a non-trivial
# fraction, else it is sequencing/contamination noise. Mirrors MixInfect2's QC
# (min minor-allele depth 5, min minor-allele proportion 0.05). Fraction bounds
# are written as cross-multiplications (alt >= AF*(ref+alt)) so there is no
# division in the bcftools expression -> no parser ambiguity / div-by-zero.
# NOTE: this cleans SEQUENCING-ERROR het only; it does NOT remove divergence-vs-
# H99 het (those have real read support) -> that is handled downstream by popFreq
# (shared) + k-mer ploidy / single-colony (private). Re-run 04 on ALL batches
# (all + allmono) so old and new share one ruler; then re-merge + mask_het_table.sh.
# Set APPLY_HET_QC=0 to reproduce the original (genotype-only) behavior.
APPLY_HET_QC="${APPLY_HET_QC:-1}"
MIN_DP="${MIN_DP:-10}"               # total depth floor
MIN_MINOR_DP="${MIN_MINOR_DP:-5}"    # each allele >= this many reads (=> minor allele >= 5)
MIN_MINOR_AF="${MIN_MINOR_AF:-0.05}" # minor-allele fraction floor (alt frac kept in [AF, 1-AF])
HET_FILTER='GT="0/1" && TYPE="snp" && N_ALT=1'
if [[ "$APPLY_HET_QC" == "1" ]]; then
  MAX_MINOR_AF="$(awk "BEGIN{printf \"%.4f\", 1 - $MIN_MINOR_AF}")"
  HET_FILTER="${HET_FILTER} && FMT/DP>=${MIN_DP}"
  HET_FILTER="${HET_FILTER} && FMT/AD[0:0]>=${MIN_MINOR_DP} && FMT/AD[0:1]>=${MIN_MINOR_DP}"
  HET_FILTER="${HET_FILTER} && FMT/AD[0:1]>=${MIN_MINOR_AF}*(FMT/AD[0:0]+FMT/AD[0:1])"
  HET_FILTER="${HET_FILTER} && FMT/AD[0:1]<=${MAX_MINOR_AF}*(FMT/AD[0:0]+FMT/AD[0:1])"
  echo "[$(date '+%F %T')] >>> Het QC ON : DP>=${MIN_DP}, each-allele>=${MIN_MINOR_DP}, alt-frac in [${MIN_MINOR_AF},${MAX_MINOR_AF}]"
else
  echo "[$(date '+%F %T')] >>> Het QC OFF: APPLY_HET_QC=0 -> genotype-only (original behavior)"
fi
echo "[$(date '+%F %T')] >>> het -i filter: ${HET_FILTER}"

[[ -d "$VCF_DIR" ]] || { echo "[ERROR] filtered VCF dir not found: $VCF_DIR" >&2; exit 1; }
VCFS=( "${VCF_DIR}"/*.filtered.vcf.gz )
(( ${#VCFS[@]} == 0 )) && { echo "[ERROR] No *.filtered.vcf.gz under $VCF_DIR" >&2; exit 1; }
echo "[$(date '+%F %T')] >>> Detected ${#VCFS[@]} filtered VCFs."

# ---- Parallel (FIFO) ----
THREADS="${SLURM_CPUS_PER_TASK:-16}"
MAX_JOBS="$HET_MAX_JOBS"; (( MAX_JOBS > THREADS )) && MAX_JOBS="$THREADS"
echo "[$(date '+%F %T')] >>> parallel jobs: $MAX_JOBS"
SEM_FIFO="$(mktemp -u)"; mkfifo "$SEM_FIFO"; exec 3<>"$SEM_FIFO"; rm -f "$SEM_FIFO"
for ((i=0;i<MAX_JOBS;i++)); do echo >&3; done
run_parallel() { read -u 3; ( trap 'echo >&3' EXIT; "$@" ) & }

# ---- Task ----
extract_het_one() {
  local vcf="$1" sample out_file
  sample="$(basename "$vcf" .filtered.vcf.gz)"
  out_file="${HET_OUT}/${sample}.het.tsv"

  # "${QUERY_TARGETS[@]}" is empty when no mask -> query stays genome-wide.
  # MASK_BED is named *.bed so bcftools -T treats it as 0-based (matches POS).
  if ! bcftools query \
        "${QUERY_TARGETS[@]}" \
        -i "$HET_FILTER" \
        -f "${sample}\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT]\t[%DP]\t[%AD]\n" \
        "$vcf" > "$out_file" 2>"${HET_LOG_DIR}/${sample}.err"; then
    echo "[ERROR] bcftools query failed: $sample"
    touch "$FAIL_DIR/$sample"; rm -f "$out_file"; return 1
  fi
  local n; n=$(wc -l < "$out_file")
  echo "[$(date '+%F %T')] Done: $sample ($n het SNPs)"
}

for vcf in "${VCFS[@]}"; do
  run_parallel extract_het_one "$vcf"
done
wait || true
exec 3>&- 3<&- || true

# ---- Failure report ----
mapfile -t FAILED < <(find "$FAIL_DIR" -maxdepth 1 -type f -printf '%f\n' 2>/dev/null | sort)
if (( ${#FAILED[@]} > 0 )); then
  printf '%s\n' "${FAILED[@]}" > "$HET_OUT/FAILED_SAMPLES.txt"
  echo "=============================================================="
  echo "[WARN] ${#FAILED[@]} sample(s) FAILED het extraction:"
  printf '   %s\n' "${FAILED[@]}"
  echo "       full list: $HET_OUT/FAILED_SAMPLES.txt"
  echo "=============================================================="
fi

# ---- Merge (guarded: empty glob would otherwise make 'cat' read stdin) ----
echo "[$(date '+%F %T')] >>> Merging per-sample tables..."
PART_FILES=( "${HET_OUT}"/*.het.tsv )
if (( ${#PART_FILES[@]} == 0 )); then
  echo "[ERROR] No per-sample het files produced; nothing to merge." >&2
  exit 1
fi

OUT_TSV="${HET_OUT}/all_het_snps.tsv"
echo -e "sample\tchrom\tpos\tref\talt\tqual\tfilter_st\tgt\tdp\tad" > "$OUT_TSV"
cat "${PART_FILES[@]}" >> "$OUT_TSV"
rm -f "${PART_FILES[@]}"

TOTAL=$(( $(wc -l < "$OUT_TSV") - 1 ))
echo "[$(date '+%F %T')] >>> All done! Total het SNPs: $TOTAL"
echo "[$(date '+%F %T')] >>> Output: $OUT_TSV"
echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
exit 0
