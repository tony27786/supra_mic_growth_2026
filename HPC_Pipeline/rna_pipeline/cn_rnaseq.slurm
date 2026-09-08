#!/bin/bash
#SBATCH -J cn_rnaseq                              # [CN-CHANGE] job name
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -N 1
#SBATCH --cpus-per-task=48
#SBATCH --mem=400G                                # [CN-NOTE] huge headroom; this run peaked at ~10 GB. Safe to drop to ~96G if scheduling is slow.
#SBATCH --output=%x_%j.log

# =============================================================================
# fastp -> STAR -> featureCounts   (paired-end)
# [CN-CHANGE] Target: Cryptococcus neoformans var. grubii H99
#             genome/annotation = NCBI RefSeq GCF_000149245.1 (CNA3)
#             SAME reference as the Aim 3 WGS pipeline (NC_ seqids), so RNA-seq
#             and WGS share one coordinate system.
# Reads samples.tsv: sample <tab> fastq_1 <tab> fastq_2   (flat RAW_DIR)
# Final gene count matrix -> read_counts/gene_counts.txt   (gene_id = CNAG_xxxxx)
# =============================================================================

# ---- Strict mode ------------------------------------------------------------
set -eo pipefail
trap 'code=$?; echo -e "\n[ERROR] Exit $code at line $LINENO" >&2' ERR
IFS=$'\n\t'
shopt -s nullglob

# ---- Initialize Conda environment -------------------------------------------
CONDA_BASE="$(${CONDA_LAUNCHER} info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
set +u
conda activate prernaseq
set -u
echo "CONDA_PREFIX=$CONDA_PREFIX"
echo "==> Active conda environment: ${CONDA_DEFAULT_ENV:-unknown}"
echo "==> conda list (full package inventory of the activated env):"
conda list
echo "==> Key tool versions:"
for t in fastp STAR featureCounts samtools; do
  if command -v "$t" >/dev/null 2>&1; then
    printf '    %-14s %s\n' "$t" "$(command -v "$t")"
  else
    echo "    [MISSING] $t not found in PATH (env may be wrong)"; exit 1
  fi
done
echo ""

# ---- SITE SETTINGS -- export these before submitting -------------------------
#   export WORK_ROOT=/scratch/$USER/hetero
#   export RAW_DIR=/data/hetero/rna_fastq      # {sample}_R1.fastq.gz / _R2
#   export REFDIR=/refs                        # holds the two GCF_000149245.1 files
#   export SLURM_PARTITION=<your partition>
# Sample names must match ../manifest/rna_samples.tsv (60 libraries); build
# samples.tsv with make_sample_tsv.sh before submitting.
WORK_ROOT="${WORK_ROOT:?set WORK_ROOT to a writable scratch directory}"
PROJ="${PROJ:-$WORK_ROOT/cn_rna}"                  # this run's project dir
RAW_DIR="${RAW_DIR:?set RAW_DIR to the directory holding the RNA-seq FASTQ pairs}"
SAMPLES="$PROJ/samples.tsv"
CLEAN="$PROJ/clean_fastq"
ALIGN="$PROJ/star_align"
COUNTS="$PROJ/read_counts"                        # <-- raw counts
# Reference: NCBI RefSeq GCF_000149245.1 (C. neoformans H99, assembly CNA3).
# Same reference and coordinate system as the WGS pipeline.
REFDIR="${REFDIR:-$WORK_ROOT/h99ref}"
GENOME_FASTA="$REFDIR/GCF_000149245.1_CNA3_genomic.fna"
GTF_FILE="$REFDIR/GCF_000149245.1_CNA3_genomic.gtf"
STAR_INDEX="$PROJ/star_index"
STRAND=0                                          # [CN-NOTE] 0=unstranded 1=fwd 2=rev -- VERIFY after first run (see footer); many dUTP kits are 2
BAMSORT_RAM=15000000000                            # [CN-FIX] explicit BAM sort RAM per STAR job (~8 GB); the small H99 genome makes STAR's default (=genome index size) too small -> "not enough memory for BAM sorting"

THREADS_TOTAL="${SLURM_CPUS_PER_TASK:-48}"
MAX_JOBS_FASTP=8                                  # fastp concurrency
MAX_JOBS_STAR=8                                   # [CN-NOTE] per-job peak this run: ~1.4 GB mapping + up to BAMSORT_RAM sorting; 8 * 8 GB = 64 GB << --mem. Safe to raise.
FASTP_STAGGER=3                                   # seconds between fastp launches (avoid concurrent-finalize deadlock)
STAR_STAGGER=2                                    # seconds between STAR launches (smooth the genome-load I/O burst)
# -----------------------------------------------------------------------------

mkdir -p "$CLEAN" "$ALIGN" "$COUNTS" "$REFDIR"

echo "==> Resolved configuration:"
echo "    PROJ        = $PROJ"
echo "    RAW_DIR     = $RAW_DIR"
echo "    SAMPLES     = $SAMPLES"
echo "    GENOME_FASTA= $GENOME_FASTA"
echo "    GTF_FILE    = $GTF_FILE"
echo "    STAR_INDEX  = $STAR_INDEX"
echo "    COUNTS      = $COUNTS"
echo "    STRAND      = $STRAND   (0=unstranded 1=fwd 2=rev)"
echo "    BAMSORT_RAM = $BAMSORT_RAM bytes per STAR job"
echo "    THREADS_TOTAL=$THREADS_TOTAL  MAX_JOBS_FASTP=$MAX_JOBS_FASTP  MAX_JOBS_STAR=$MAX_JOBS_STAR"
[[ -s "$SAMPLES" ]] || { echo "[FATAL] samples.tsv not found: $SAMPLES"; exit 1; }
N_SAMPLES=$(($(wc -l < "$SAMPLES") - 1))
echo "    samples in sheet: $N_SAMPLES"
echo ""
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ts(){ date '+%F %T'; }

# ---- FIFO semaphore ----------------------------------------------------------
SEM_FIFO=""; exec 3>&- 2>/dev/null || true
sem_init(){
  local n="$1"
  [[ -n "$SEM_FIFO" ]] && exec 3>&- 3<&- 2>/dev/null || true
  SEM_FIFO="$(mktemp -u)"; mkfifo "$SEM_FIFO"; exec 3<>"$SEM_FIFO"; rm -f "$SEM_FIFO"
  local i; for ((i=0;i<n;i++)); do echo >&3; done
}
run_parallel(){ read -u 3; ( trap 'echo >&3' EXIT; "$@" ) & }

# =============================================================================
# STEP 0  Build STAR index (once; auto-skip if present, auto-detect overhang)
# [CN-CHANGE] --genomeSAindexNbases 11 for the small (~19 Mb) H99 genome:
#            min(14, log2(GenomeLength)/2 - 1) = min(14, ~11.1) = 11.
# =============================================================================
if [[ ! -s "$STAR_INDEX/SAindex" ]]; then
  echo "[$(ts)] >>> Building STAR index"
  echo "[$(ts)]     STAR version: $(STAR --version)"
  [[ -s "$GENOME_FASTA" ]] || { echo "[FATAL] Missing genome fasta: $GENOME_FASTA (gunzip it first)"; exit 1; }
  [[ -s "$GTF_FILE" ]]     || { echo "[FATAL] Missing GTF: $GTF_FILE (gunzip it first)"; exit 1; }
  echo "[$(ts)]     genome fasta: $GENOME_FASTA ($(du -h "$GENOME_FASTA" | cut -f1))"
  echo "[$(ts)]     gtf file    : $GTF_FILE ($(du -h "$GTF_FILE" | cut -f1))"
  echo "[$(ts)]     memory available:"; free -h | sed 's/^/        /'

  # read length from first sample -> sjdbOverhang = readlen-1
  # NOTE: awk reads via process substitution and exits itself, so the upstream
  #       zcat SIGPIPE never propagates through a pipeline -> safe under pipefail.
  first_r1=$(tail -n +2 "$SAMPLES" | head -n 1 | cut -f2)
  echo "[$(ts)]     probing read length from: $RAW_DIR/$first_r1"
  [[ -s "$RAW_DIR/$first_r1" ]] || { echo "[FATAL] first R1 not found: $RAW_DIR/$first_r1"; exit 1; }
  RLEN=$(awk 'NR%4==2{if(length>m)m=length} NR>=4000{exit} END{print m+0}' <(zcat -- "$RAW_DIR/$first_r1" 2>/dev/null))
  [[ "${RLEN:-0}" -gt 0 ]] || RLEN=150
  OVERHANG=$(( RLEN - 1 ))
  echo "[$(ts)]     detected read length=${RLEN}, sjdbOverhang=${OVERHANG}"

  mkdir -p "$STAR_INDEX"
  echo "[$(ts)]     launching genomeGenerate with $THREADS_TOTAL threads ..."
  STAR --runThreadN "$THREADS_TOTAL" --runMode genomeGenerate \
       --genomeDir "$STAR_INDEX" --genomeFastaFiles "$GENOME_FASTA" \
       --sjdbGTFfile "$GTF_FILE" --sjdbOverhang "$OVERHANG" \
       --genomeSAindexNbases 11 \
       --outFileNamePrefix "$STAR_INDEX/_genomeGenerate."
  [[ -s "$STAR_INDEX/SAindex" ]] || { echo "[FATAL] STAR index build failed (no SAindex). See $STAR_INDEX/_genomeGenerate.Log.out"; exit 1; }
  echo "[$(ts)] >>> STAR index built OK -> $STAR_INDEX"
else
  echo "[$(ts)] >>> STAR index exists, skip build"
fi

# =============================================================================
# STEP 1  fastp (parallel)
# =============================================================================
fastp_job(){
  local s="$1" r1="$2" r2="$3"
  local od="$CLEAN/$s"; mkdir -p "$od"
  local o1="$od/${s}_R1.clean.fq.gz" o2="$od/${s}_R2.clean.fq.gz"
  if [[ -s "$o1" && -s "$o2" ]]; then echo "[$(ts)] [skip fastp] $s"; return 0; fi
  echo "[$(ts)] [fastp] $s"
  fastp -i "$RAW_DIR/$r1" -I "$RAW_DIR/$r2" -o "$o1" -O "$o2" \
        -w "$THREADS_PER_JOB" -h "$od/${s}.fastp.html" -j "$od/${s}.fastp.json"
}
echo "[$(ts)] === STAGE 1: fastp ==="
THREADS_PER_JOB=$(( (THREADS_TOTAL + MAX_JOBS_FASTP - 1) / MAX_JOBS_FASTP )); ((THREADS_PER_JOB<1))&&THREADS_PER_JOB=1
export THREADS_PER_JOB
sem_init "$MAX_JOBS_FASTP"
fastp_pids=()
while IFS=$'\t' read -r s r1 r2; do
  [[ -z "$s" ]] && continue
  run_parallel fastp_job "$s" "$r1" "$r2"
  fastp_pids+=("$!")
  sleep "$FASTP_STAGGER"     # stagger launches; semaphore still caps concurrency at MAX_JOBS_FASTP
done < <(tail -n +2 "$SAMPLES")
fp_fail=0
for pid in "${fastp_pids[@]}"; do wait "$pid" || fp_fail=1; done
(( fp_fail == 0 )) || { echo "[FATAL] one or more fastp jobs failed"; exit 1; }
echo "[$(ts)] STAGE 1 complete: ${#fastp_pids[@]} fastp jobs finished."

# =============================================================================
# STEP 2  STAR alignment (parallel)
# [CN-CHANGE] --alignIntronMax 3000 : H99 introns are tiny (median 56 bp, longest
#            real intron 2124 bp across 45,622 introns; 0 exceed 3000). Covers
#            every real intron while blocking cross-gene spurious spliced aligns.
#            --alignIntronMin left at STAR default 21 (= observed minimum intron).
# [CN-CHANGE] --quantMode GeneCounts : emits ReadsPerGene.out.tab so the footer
#            strandedness self-check works + free cross-check vs featureCounts.
# [CN-FIX]   --limitBAMsortRAM $BAMSORT_RAM : without it, on this 19 Mb genome
#            STAR sets the sort buffer to the (tiny) genome index size, and
#            larger samples die with "not enough memory for BAM sorting" (exit 1,
#            no signal). This was the cause of 43/57 failures in run 2021907.
# [CN-FIX]   per-sample isolation: a failed sample is recorded in FAILED_SAMPLES.txt
#            and does NOT abort the batch; the run reports the full failure list.
# =============================================================================
FAILED_LIST="$ALIGN/FAILED_SAMPLES.txt"
: > "$FAILED_LIST"                                  # reset for this run

star_job(){
  local s="$1"
  local id="$CLEAN/$s"
  local r1=( "$id"/*_R1.clean.fq.gz ) r2=( "$id"/*_R2.clean.fq.gz )
  local od="$ALIGN/$s"; mkdir -p "$od"
  local bam="$od/${s}.Aligned.sortedByCoord.out.bam"
  if [[ -s "$bam" ]]; then echo "[$(ts)] [skip STAR] $s"; return 0; fi
  if (( ${#r1[@]} == 0 || ${#r2[@]} == 0 )); then
    echo "[$(ts)] [STAR] $s : missing clean fastq in $id" >&2; return 1
  fi
  echo "[$(ts)] [STAR] $s"
  STAR --runThreadN "$THREADS_PER_JOB" --genomeDir "$STAR_INDEX" \
       --readFilesIn "${r1[0]}" "${r2[0]}" --readFilesCommand zcat \
       --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
       --alignIntronMax 3000 \
       --limitBAMsortRAM "$BAMSORT_RAM" \
       --quantMode GeneCounts \
       --outFileNamePrefix "$od/${s}."
}

# isolate each sample: record failures, never abort the batch (Track A pattern)
star_one(){
  local s="$1"
  if star_job "$s"; then
    return 0
  fi
  echo "$s" >> "$FAILED_LIST"
  echo "[$(ts)] [STAR FAILED] $s  (see $ALIGN/$s/$s.Log.out)" >&2
  return 0
}

echo "[$(ts)] === STAGE 2: STAR ==="
THREADS_PER_JOB=$(( (THREADS_TOTAL + MAX_JOBS_STAR - 1) / MAX_JOBS_STAR )); ((THREADS_PER_JOB<1))&&THREADS_PER_JOB=1
export THREADS_PER_JOB
sem_init "$MAX_JOBS_STAR"
star_pids=()
while IFS=$'\t' read -r s r1 r2; do
  [[ -z "$s" ]] && continue
  run_parallel star_one "$s"
  star_pids+=("$!")
  sleep "$STAR_STAGGER"      # spread out the genome-load bursts
done < <(tail -n +2 "$SAMPLES")
for pid in "${star_pids[@]}"; do wait "$pid"; done   # star_one always exits 0; failures are in FAILED_LIST
exec 3>&- 3<&- 2>/dev/null || true

if [[ -s "$FAILED_LIST" ]]; then
  n_fail=$(wc -l < "$FAILED_LIST")
  echo "[$(ts)] [WARN] ${n_fail} STAR job(s) FAILED:"
  sed 's/^/        - /' "$FAILED_LIST"
  echo "[FATAL] Refusing to run featureCounts with ${n_fail} failed sample(s)."
  echo "        Inspect each <sample>.Log.out, fix, and re-submit -- completed samples are skipped automatically."
  exit 1
fi
echo "[$(ts)] STAGE 2 complete: all ${#star_pids[@]} STAR jobs finished, 0 failures."

# =============================================================================
# STEP 3  featureCounts (single multithreaded run over all BAMs)
# [CN-NOTE] -t exon -g gene_id is UNCHANGED: standard GTF2.2 with exon features
#           and gene_id "CNAG_xxxxx" on every line -> counts per CNAG locus tag.
# =============================================================================
echo "[$(ts)] === STAGE 3: featureCounts ==="
mapfile -t BAMS < <(find "$ALIGN" -type f -name "*.Aligned.sortedByCoord.out.bam" | sort)
echo "[$(ts)] Found ${#BAMS[@]} BAMs (expected $N_SAMPLES)"
if (( ${#BAMS[@]} != N_SAMPLES )); then
  echo "[FATAL] BAM count (${#BAMS[@]}) != expected sample count ($N_SAMPLES)."
  echo "        Refusing to run featureCounts on an incomplete BAM set."
  echo "        Inspect $ALIGN for failed/incomplete STAR jobs and re-run."
  exit 1
fi

featureCounts -T "$THREADS_TOTAL" -s "$STRAND" \
  -p --countReadPairs \
  -t exon -g gene_id \
  -a "$GTF_FILE" \
  -o "$COUNTS/gene_counts.txt" \
  "${BAMS[@]}"

echo "[$(ts)] >>> DONE. Count matrix: $COUNTS/gene_counts.txt"
echo "[$(ts)] >>> Strandedness check: for each sample compare the column sums of"
echo "          ${ALIGN}/<sample>/<sample>.ReadsPerGene.out.tab"
echo "          col2=unstranded, col3=forward, col4=reverse; the largest column is the true strandedness."
echo "          If col4 dominates -> set STRAND=2 and re-run STAGE 3 (featureCounts) only."
echo ">>> SLURM BATCH TERMINATED AT [$(ts)]"
