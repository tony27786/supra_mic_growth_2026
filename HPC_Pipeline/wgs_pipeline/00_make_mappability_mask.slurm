#!/bin/bash
#SBATCH -J cn00_mask
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -o 00_make_mappability_mask.out
# =============================================================================
# 00_make_mappability_mask.slurm  --  ONE-TIME reference-level prep.
#
# Builds a BED of uniquely-mappable (single-copy) regions of REF_FA using GenMap.
# This is NOT part of the per-sample 00->04 chain: run it ONCE per reference
# (like building the bwa/faidx index). The resulting $MASK_BED is reused by:
#   - stage 04 (04_extract_het.slurm), via `bcftools query -T $MASK_BED`
#   - the standalone mask_het_table.sh, via `bedtools intersect -b $MASK_BED`
#
# Purpose: het counted against H99 is dominated by reads that mis-map at
# paralogous / repetitive / divergent loci (fake biallelic sites). Restricting
# het counting to the unique fraction removes that artifact baseline, so a clean
# clone trends toward ~0 het and only genuine within-sample heterozygosity /
# mixture survives.
#
# Coordinates: GenMap bedGraph and bedtools merge are 0-based half-open; the
# output is named *.bed so bcftools `-T` also treats it as 0-based. Consistent
# with the het table's 1-based POS once mask_het_table.sh does pos-1 -> pos.
#
# Run:
#   cd ~/fig6_wgs_pipeline
#   sbatch 00_make_mappability_mask.slurm
#   # tune k / errors / strictness if needed:
#   #   export MASK_K=150 MASK_E=1 MASK_MIN_MAPPABILITY=1
# =============================================================================

PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${PIPELINE_DIR}/00_config.sh"
apply_strict_mode

activate_conda "$CONDA_ENV_BWASAM"   # genmap + bedtools live here
print_config_banner

echo "[$(date '+%F %T')] >>> MASK_K               = $MASK_K"
echo "[$(date '+%F %T')] >>> MASK_E               = $MASK_E"
echo "[$(date '+%F %T')] >>> MASK_MIN_MAPPABILITY = $MASK_MIN_MAPPABILITY"
echo "[$(date '+%F %T')] >>> GENMAP_DIR           = $GENMAP_DIR"
echo "[$(date '+%F %T')] >>> GENMAP_INDEX         = $GENMAP_INDEX"
echo "[$(date '+%F %T')] >>> GENMAP_PREFIX        = $GENMAP_PREFIX"
echo "[$(date '+%F %T')] >>> MASK_BED (output)    = $MASK_BED"

command -v genmap   >/dev/null 2>&1 || { echo "[ERROR] genmap not on PATH (env $CONDA_ENV_BWASAM)"   >&2; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "[ERROR] bedtools not on PATH (env $CONDA_ENV_BWASAM)" >&2; exit 1; }
[[ -f "$REF_FA" ]] || { echo "[ERROR] REF_FA not found: $REF_FA" >&2; exit 1; }

mkdir -p "$GENMAP_DIR"

# samtools .fai is needed for the genome-fraction report (built by 01_map on
# first run; build it here too if missing so this script is self-contained).
if [[ ! -f "${REF_FA}.fai" ]]; then
  echo "[$(date '+%F %T')] >>> building ${REF_FA}.fai"
  samtools faidx "$REF_FA"
fi

# ---- 1) GenMap index -------------------------------------------------------
# genmap index refuses to write into an existing directory, so clear it first.
# Guard the rm so a misconfigured path can never nuke something unrelated.
case "$GENMAP_INDEX" in
  *genmap*) : ;;
  *) echo "[ERROR] GENMAP_INDEX path looks unsafe, refusing to rm: $GENMAP_INDEX" >&2; exit 1 ;;
esac
rm -rf "$GENMAP_INDEX"

echo "[$(date '+%F %T')] >>> genmap index ..."
genmap index -F "$REF_FA" -I "$GENMAP_INDEX"

# ---- 2) GenMap mappability (bedGraph) --------------------------------------
# -K = k-mer length (match read length); -E = mismatches allowed when deciding
# whether a k-mer is unique. Output: ${GENMAP_PREFIX}.bedgraph (chrom start end mappability).
echo "[$(date '+%F %T')] >>> genmap map -K $MASK_K -E $MASK_E ..."
genmap map -K "$MASK_K" -E "$MASK_E" -I "$GENMAP_INDEX" -O "$GENMAP_PREFIX" -bg

BG="${GENMAP_PREFIX}.bedgraph"
[[ -s "$BG" ]] || { echo "[ERROR] GenMap produced no bedgraph: $BG" >&2; exit 1; }

# ---- 3) Threshold to unique regions + merge -> MASK_BED --------------------
# mappability == 1 means the k-mer occurs exactly once (strictly unique). Keep
# positions with mappability >= MASK_MIN_MAPPABILITY (default 1), then merge
# adjacent intervals.
echo "[$(date '+%F %T')] >>> threshold (>= $MASK_MIN_MAPPABILITY) + merge -> MASK_BED"
awk -v m="$MASK_MIN_MAPPABILITY" 'BEGIN{OFS="\t"} ($4+0)>=m {print $1, $2, $3}' "$BG" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - > "$MASK_BED"

[[ -s "$MASK_BED" ]] || { echo "[ERROR] empty MASK_BED after thresholding: $MASK_BED" >&2; exit 1; }

# ---- 4) Report fraction of genome retained ---------------------------------
MASK_BP=$(awk '{s += $3-$2} END{print s+0}' "$MASK_BED")
GENOME_BP=$(awk '{s += $2} END{print s+0}' "${REF_FA}.fai")
echo "=============================================================="
echo "[$(date '+%F %T')] >>> mappable bp        : $MASK_BP"
echo "[$(date '+%F %T')] >>> genome bp (.fai)   : $GENOME_BP"
awk -v a="$MASK_BP" -v b="$GENOME_BP" 'BEGIN{ if(b>0) printf ">>> uniquely-mappable : %.2f%% of genome\n", 100*a/b }'
echo ">>> MASK_BED = $MASK_BED"
echo "=============================================================="
echo ">>> Next: run  bash mask_het_table.sh <merged all_het_snps.tsv> $MASK_BED <outdir>"
echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
exit 0
