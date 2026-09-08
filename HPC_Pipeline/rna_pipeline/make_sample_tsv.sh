#!/usr/bin/env bash
# =============================================================================
# make_sample_tsv.sh -- build samples.tsv for cn_rnaseq.slurm.
#
#   RAW_DIR=/path/to/rna_fastq bash make_sample_tsv.sh
#
# Expects one pair per library in a flat directory:
#   {sample}_R1.fastq.gz   and   {sample}_R2.fastq.gz
#
# Writes samples.tsv next to this script, with the header row that
# cn_rnaseq.slurm expects (it reads the file with `tail -n +2`):
#   sample <tab> fastq_1 <tab> fastq_2
#
# Sample names must match the `sample` column of ../manifest/rna_samples.tsv
# (60 libraries: 20 isolates x Monoclonal/Selected/Passage).
# =============================================================================
set -euo pipefail

RAW_DIR="${RAW_DIR:?set RAW_DIR to the directory holding the RNA-seq FASTQ pairs}"
[[ -d "$RAW_DIR" ]] || { echo "[ERROR] RAW_DIR not found: $RAW_DIR" >&2; exit 1; }

OUT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/samples.tsv"

{
  printf "sample\tfastq_1\tfastq_2\n"
  find "$RAW_DIR" -maxdepth 1 -type f -name "*_R1.fastq.gz" -printf "%f\n" \
    | sort \
    | while read -r r1; do
        sample="${r1%_R1.fastq.gz}"
        r2="${sample}_R2.fastq.gz"
        if [[ -f "$RAW_DIR/$r2" ]]; then
          printf "%s\t%s\t%s\n" "$sample" "$r1" "$r2"
        else
          echo "WARNING: missing R2 for $r1" >&2
        fi
      done
} > "$OUT"

echo "written: $OUT"
echo "libraries: $(( $(wc -l < "$OUT") - 1 ))  (expected 60)"
head -3 "$OUT"
