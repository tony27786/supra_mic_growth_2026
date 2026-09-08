#!/bin/bash
#SBATCH -J cn01_map
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -N 1
#SBATCH --mem=300G
#SBATCH --cpus-per-task=48
#SBATCH -o 01_cn_map_bwa.out
# =============================================================================
# Stage 01: bwa mem -> samtools sort -> index -> flagstat/stats/idxstats
# Outputs: $MAP_OUT/{sample}.bam(.bai)  +  $MAP_OUT/qc_metrics_map/*
# =============================================================================

# ---- Load shared config (PIPELINE_DIR = where sbatch was invoked) ----
PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${PIPELINE_DIR}/00_config.sh"
apply_strict_mode

THREADS="${SLURM_CPUS_PER_TASK:-4}"

# ---- Output layout ----
METRICS_OUT="$MAP_OUT/qc_metrics_map"
FLAGSTAT_DIR="$METRICS_OUT/flagstat"
STATS_DIR="$METRICS_OUT/stats"
IDXSTATS_DIR="$METRICS_OUT/idxstats"
LOG_DIR="$METRICS_OUT/logs"
FAIL_DIR="$MAP_OUT/.failed"
mkdir -p "$MAP_OUT" "$METRICS_OUT" "$FLAGSTAT_DIR" "$STATS_DIR" \
         "$IDXSTATS_DIR" "$LOG_DIR" "$FAIL_DIR"
rm -f "$FAIL_DIR"/* 2>/dev/null || true

cd "$MAP_OUT"
echo "Working Directory = $(pwd)"
activate_conda "$CONDA_ENV_BWASAM"
print_config_banner
echo "[$(date '+%F %T')] >>> SLURM_CPUS  = $THREADS"

# ---- Dependencies & reference ----
command -v bwa       >/dev/null 2>&1 || { echo "[ERROR] bwa not in PATH" >&2; exit 1; }
command -v samtools  >/dev/null 2>&1 || { echo "[ERROR] samtools not in PATH" >&2; exit 1; }
[[ -f "$REF_FA" ]] || { echo "[ERROR] Reference FASTA not found: $REF_FA" >&2; exit 1; }

echo "[$(date '+%F %T')] >>> Checking reference index"
if [[ ! -f "${REF_FA}.bwt" && ! -f "${REF_FA}.0123" ]]; then
  echo "[$(date '+%F %T')] >>> bwa index"
  bwa index "$REF_FA"
else
  echo "[$(date '+%F %T')] >>> bwa index exists, skip"
fi
if [[ ! -f "${REF_FA}.fai" ]]; then
  echo "[$(date '+%F %T')] >>> samtools faidx"
  samtools faidx "$REF_FA"
else
  echo "[$(date '+%F %T')] >>> fai exists, skip"
fi

# ---- Per-job thread budget ----
MAX_JOBS="$MAP_MAX_JOBS"
JOB_CPUS=$(( THREADS / MAX_JOBS )); (( JOB_CPUS < 1 )) && JOB_CPUS=1
SORT_THREADS="$MAP_SORT_THREADS"; (( JOB_CPUS <= SORT_THREADS )) && SORT_THREADS=1
BWA_THREADS=$(( JOB_CPUS - SORT_THREADS )); (( BWA_THREADS < 1 )) && BWA_THREADS=1
INDEX_THREADS="$MAP_INDEX_THREADS"; (( INDEX_THREADS > JOB_CPUS )) && INDEX_THREADS="$JOB_CPUS"
THREADS_PER_JOB="$JOB_CPUS"
echo "[$(date '+%F %T')] >>> Total CPUs: $THREADS; parallel jobs: $MAX_JOBS"
echo "[$(date '+%F %T')] >>> Per-job: bwa=$BWA_THREADS sort=$SORT_THREADS index=$INDEX_THREADS"

# ---- FIFO semaphore ----
SEM_FIFO="$(mktemp -u)"; mkfifo "$SEM_FIFO"; exec 3<>"$SEM_FIFO"; rm -f "$SEM_FIFO"
for ((i=0;i<MAX_JOBS;i++)); do echo >&3; done
run_parallel() {
  read -u 3
  (
    trap 'echo >&3' EXIT
    export OMP_NUM_THREADS="$THREADS_PER_JOB" THREADS="$THREADS_PER_JOB"
    "$@"
  ) &
}

# ---- Discover paired FASTQs ----
cd "$CLEAN_DIR" || { echo "[ERROR] CLEAN_DIR not found: $CLEAN_DIR" >&2; exit 1; }
SAMPLES=(); R1_LIST=(); R2_LIST=()
for d in */; do
  sample="${d%/}"
  r1="${CLEAN_DIR}/${sample}/${sample}_clean_R1.fq.gz"
  r2="${CLEAN_DIR}/${sample}/${sample}_clean_R2.fq.gz"
  if [[ -f "$r1" && -f "$r2" ]]; then
    SAMPLES+=("$sample"); R1_LIST+=("$r1"); R2_LIST+=("$r2")
  else
    echo "[WARN] Missing R1/R2 for sample: $sample"
  fi
done
cd - >/dev/null

TOTAL_JOBS="${#SAMPLES[@]}"
echo "[$(date '+%F %T')] >>> Detected ${TOTAL_JOBS} samples for mapping."
(( TOTAL_JOBS == 0 )) && { echo "[ERROR] No samples detected." >&2; exit 1; }

# ---- Progress bar ----
PROGRESS_WIDTH=45
COUNTER_FILE="$(mktemp)"; LOCK_FILE="$(mktemp)"; echo 0 > "$COUNTER_FILE"
update_progress() {
  local sample="$1"
  {
    flock -x 9
    local done_jobs; done_jobs="$(<"$COUNTER_FILE")"; done_jobs=$((done_jobs + 1))
    echo "$done_jobs" > "$COUNTER_FILE"
    local percent=$(( done_jobs * 100 / TOTAL_JOBS ))
    local bar_len=$(( done_jobs * PROGRESS_WIDTH / TOTAL_JOBS )) bar="" i
    for ((i=0; i<bar_len; i++)); do bar+="#"; done
    for ((i=bar_len; i<PROGRESS_WIDTH; i++)); do bar+="."; done
    printf "[PROGRESS] [%s] %3d%% (%d/%d) latest: %s\n" \
      "$bar" "$percent" "$done_jobs" "$TOTAL_JOBS" "$sample"
  } 9>>"$LOCK_FILE"
}

# ---- Mapping task ----
map_job() {
  local sample="$1" r1="$2" r2="$3"
  local bam_out="$MAP_OUT/${sample}.bam"
  local flag_out="$FLAGSTAT_DIR/${sample}.flagstat.txt"
  local stats_out="$STATS_DIR/${sample}.stats.txt"
  local idxstats_out="$IDXSTATS_DIR/${sample}.idxstats.txt"
  local log_out="$LOG_DIR/${sample}.map.log"; : > "$log_out"
  local rg="@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA"
  local tmp_dir="$MAP_OUT/tmp/${sample}"; mkdir -p "$tmp_dir"

  echo "[$(date '+%F %T')] >>> START map: $sample" | tee -a "$log_out"

  # --- To resume a partial run, uncomment to skip completed samples: ---
  # if [[ -s "$bam_out" && -s "${bam_out}.bai" && -s "$flag_out" ]]; then
  #   echo "[$(date '+%F %T')] >>> skip (done): $sample"; update_progress "$sample"; return 0
  # fi

  if ! bwa mem -R "$rg" -t "$BWA_THREADS" "$REF_FA" "$r1" "$r2" 2>>"$log_out" \
       | samtools sort -@ "$SORT_THREADS" -T "${tmp_dir}/${sample}" -o "$bam_out" - 2>>"$log_out"; then
    echo "[ERROR] bwa/sort failed: $sample" | tee -a "$log_out"
    touch "$FAIL_DIR/$sample"; rm -rf "$tmp_dir"; update_progress "$sample"; return 1
  fi
  if ! samtools index -@ "$INDEX_THREADS" "$bam_out" 2>>"$log_out"; then
    echo "[ERROR] index failed: $sample" | tee -a "$log_out"
    touch "$FAIL_DIR/$sample"; rm -rf "$tmp_dir"; update_progress "$sample"; return 1
  fi
  if ! samtools quickcheck "$bam_out"; then
    echo "[ERROR] quickcheck failed: $sample" | tee -a "$log_out"
    touch "$FAIL_DIR/$sample"; rm -rf "$tmp_dir"; update_progress "$sample"; return 1
  fi

  samtools flagstat -@ "$INDEX_THREADS" "$bam_out" > "$flag_out"     2>>"$log_out" || true
  samtools stats    -@ "$INDEX_THREADS" "$bam_out" > "$stats_out"    2>>"$log_out" || true
  samtools idxstats                     "$bam_out" > "$idxstats_out" 2>>"$log_out" || true

  rm -rf "$tmp_dir"
  echo "[$(date '+%F %T')] >>> DONE  map: $sample" | tee -a "$log_out"
  update_progress "$sample"
}

# ---- Submit ----
for i in "${!SAMPLES[@]}"; do
  run_parallel map_job "${SAMPLES[$i]}" "${R1_LIST[$i]}" "${R2_LIST[$i]}"
done
wait || true
exec 3>&- 3<&- || true

# ---- Failure report ----
mapfile -t FAILED < <(find "$FAIL_DIR" -maxdepth 1 -type f -printf '%f\n' 2>/dev/null | sort)
n_ok=$(find "$MAP_OUT" -maxdepth 1 -name '*.bam' 2>/dev/null | wc -l)
if (( ${#FAILED[@]} > 0 )); then
  printf '%s\n' "${FAILED[@]}" > "$MAP_OUT/FAILED_SAMPLES.txt"
  echo "=============================================================="
  echo "[WARN] ${#FAILED[@]} sample(s) FAILED mapping:"
  printf '   %s\n' "${FAILED[@]}"
  echo "       full list: $MAP_OUT/FAILED_SAMPLES.txt"
  echo "=============================================================="
fi
echo "[$(date '+%F %T')] >>> Mapping complete. BAMs produced: ${n_ok} / ${TOTAL_JOBS}"
(( n_ok == 0 )) && { echo "[ERROR] No BAM produced; halting pipeline." >&2; exit 1; }

echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
exit 0
