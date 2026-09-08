#!/usr/bin/env bash
# =============================================================================
# submit_all.sh -- one-shot submission of the whole upstream pipeline.
#
# Freezes a single RUN_TAG, writes it to run_tag.txt, exports it to every job,
# and chains the stages with SLURM afterok dependencies:
#
#       00_fastp
#         '--- 01_map
#                |--- 02_mosdepth          (depends on 01)
#                '--- 03_vcf -> 04_het     (03 depends on 01, 04 depends on 03)
#
# Usage:
#   cd ~/fig6_wgs_pipeline
#   bash submit_all.sh                # RUN_TAG = today (yymmdd)
#   bash submit_all.sh 260702         # pin a specific RUN_TAG
#
# To override roots/params for the whole run, export before calling, e.g.:
#   export VC_MAX_DEPTH=2000
#   export RAW_DIR=/path/to/fastq
#   bash submit_all.sh
# =============================================================================
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PIPELINE_DIR"

RUN_TAG="${1:-$(date +%y%m%d)}"
echo "$RUN_TAG" > "${PIPELINE_DIR}/run_tag.txt"
echo ">>> RUN_TAG = $RUN_TAG  (written to run_tag.txt)"

# Pass current env (ALL) plus the two anchors explicitly.
EXPORTS="ALL,PIPELINE_DIR=${PIPELINE_DIR},RUN_TAG=${RUN_TAG}"

jid_fastp=$(sbatch --parsable --export="${EXPORTS}" 00_fastp.slurm)
echo ">>> 00_fastp      : $jid_fastp"

jid_map=$(sbatch --parsable --dependency=afterok:${jid_fastp} --export="${EXPORTS}" 01_map.slurm)
echo ">>> 01_map        : $jid_map   (afterok:$jid_fastp)"

jid_mos=$(sbatch --parsable --dependency=afterok:${jid_map} --export="${EXPORTS}" 02_mosdepth.slurm)
echo ">>> 02_mosdepth   : $jid_mos   (afterok:$jid_map)"

jid_vcf=$(sbatch --parsable --dependency=afterok:${jid_map} --export="${EXPORTS}" 03_vcf.slurm)
echo ">>> 03_vcf        : $jid_vcf   (afterok:$jid_map)"

jid_het=$(sbatch --parsable --dependency=afterok:${jid_vcf} --export="${EXPORTS}" 04_extract_het.slurm)
echo ">>> 04_extract_het: $jid_het   (afterok:$jid_vcf)"

WR="${WORK_ROOT:-${WORK_ROOT}}"
cat <<EOF

>>> Dependency graph
    00_fastp ($jid_fastp)
      '-- 01_map ($jid_map)
            |-- 02_mosdepth ($jid_mos)
            '-- 03_vcf ($jid_vcf) -> 04_extract_het ($jid_het)

>>> Monitor:  squeue -u \$USER
>>> VCF progress:
    tail -f ${PIPELINE_DIR}/03_cn_psvcf_snake.out \\
      | grep --line-buffered -oE '[0-9]+ of [0-9]+ steps \([0-9]+%\) done'

>>> Final results (when 02 and 04 finish):
    CNV/aneuploidy : ${WR}/cn_replenish_mosdepth_${RUN_TAG}/mosdepth/
    het SNP table  : ${WR}/cn_re_per_sample_vcf_${RUN_TAG}/het_extracted/all_het_snps.tsv
    -> move/copy these into your fig6.r working area when ready.

>>> Check each stage for partial failures:
    cat ${WR}/cn_replenish_fastpout_${RUN_TAG}/FAILED_SAMPLES.txt        2>/dev/null
    cat ${WR}/cn_replenish_map_${RUN_TAG}/FAILED_SAMPLES.txt             2>/dev/null
    cat ${WR}/cn_replenish_mosdepth_${RUN_TAG}/FAILED_SAMPLES.txt        2>/dev/null
    cat ${WR}/cn_re_per_sample_vcf_${RUN_TAG}/het_extracted/FAILED_SAMPLES.txt 2>/dev/null
EOF
