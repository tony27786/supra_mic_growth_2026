#!/bin/bash
#SBATCH -J cn00_fastp
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -N 1
#SBATCH --mem=80G
#SBATCH --cpus-per-task=40
#SBATCH -o 00_fastp.out
# =============================================================================
# Stage 00: fastp adapter/quality trimming.
# Reads : $RAW_DIR/{sample}_R1.fq.gz + {sample}_R2.fq.gz   (flat layout)
# Writes: $CLEAN_DIR/{sample}/{sample}_clean_R1.fq.gz (+ _R2, .json, .html)
# CLEAN_DIR is exactly what stage 01 (mapping) discovers, so the chain links up.
# =============================================================================

PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${PIPELINE_DIR}/00_config.sh"
apply_strict_mode

THREADS_TOTAL="${SLURM_CPUS_PER_TASK:-4}"
FASTP_OUT="$CLEAN_DIR"
FAIL_DIR="$FASTP_OUT/.failed"
mkdir -p "$FASTP_OUT" "$FAIL_DIR"
rm -f "$FAIL_DIR"/* 2>/dev/null || true

cd "$FASTP_OUT"
echo "Working Directory = $(pwd)"
activate_conda "$CONDA_ENV_FASTP"
print_config_banner

command -v fastp >/dev/null 2>&1 || { echo "[ERROR] fastp not in PATH" >&2; exit 1; }
[[ -d "$RAW_DIR" ]] || { echo "[ERROR] RAW_DIR not found: $RAW_DIR" >&2; exit 1; }

# ---- Per-job thread budget ----
MAX_JOBS="$FASTP_MAX_JOBS"
THREADS_PER_JOB=$(( (THREADS_TOTAL + MAX_JOBS - 1) / MAX_JOBS ))   # ceil
(( THREADS_PER_JOB < 1 )) && THREADS_PER_JOB=1
# fastp accepts at most 16 worker threads
(( THREADS_PER_JOB > 16 )) && THREADS_PER_JOB=16
echo "[$(date '+%F %T')] >>> CPUs: $THREADS_TOTAL  jobs: $MAX_JOBS  threads/job: $THREADS_PER_JOB"

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

# ---- Discover raw paired FASTQs ----
cd "$RAW_DIR" || { echo "[ERROR] cannot cd RAW_DIR: $RAW_DIR" >&2; exit 1; }
SAMPLES=(); R1_LIST=(); R2_LIST=()
for r1 in *_R1.fq.gz; do
  sample="${r1%_R1.fq.gz}"
  r2="${sample}_R2.fq.gz"
  if [[ -f "$r2" ]]; then
    SAMPLES+=("$sample"); R1_LIST+=("$RAW_DIR/$r1"); R2_LIST+=("$RAW_DIR/$r2")
  else
    echo "[WARN] R2 not found for sample: $sample"
  fi
done
cd - >/dev/null

TOTAL_JOBS="${#SAMPLES[@]}"
echo "[$(date '+%F %T')] >>> Detected ${TOTAL_JOBS} samples for fastp."
(( TOTAL_JOBS == 0 )) && { echo "[ERROR] No raw *_R1.fq.gz / *_R2.fq.gz pairs in $RAW_DIR" >&2; exit 1; }

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

# ---- fastp task ----
fastp_job() {
  local sample="$1" r1="$2" r2="$3"
  local odir="$FASTP_OUT/$sample"; mkdir -p "$odir"
  local log="$odir/${sample}.fastp.log"
  echo "[$(date '+%F %T')] >>> START fastp: $sample"
  if ! fastp \
        -i "$r1" -I "$r2" \
        -o "$odir/${sample}_clean_R1.fq.gz" \
        -O "$odir/${sample}_clean_R2.fq.gz" \
        -j "$odir/${sample}.json" \
        -h "$odir/${sample}.html" \
        -w "$THREADS_PER_JOB" >"$log" 2>&1; then
    echo "[ERROR] fastp failed: $sample"
    touch "$FAIL_DIR/$sample"; update_progress "$sample"; return 1
  fi
  echo "[$(date '+%F %T')] >>> DONE  fastp: $sample"
  update_progress "$sample"
}

# ---- Submit ----
for i in "${!SAMPLES[@]}"; do
  run_parallel fastp_job "${SAMPLES[$i]}" "${R1_LIST[$i]}" "${R2_LIST[$i]}"
done
wait || true
exec 3>&- 3<&- || true

# ---- Failure report ----
mapfile -t FAILED < <(find "$FAIL_DIR" -maxdepth 1 -type f -printf '%f\n' 2>/dev/null | sort)
n_ok=$(find "$FASTP_OUT" -mindepth 2 -maxdepth 2 -name '*_clean_R1.fq.gz' 2>/dev/null | wc -l)
if (( ${#FAILED[@]} > 0 )); then
  printf '%s\n' "${FAILED[@]}" > "$FASTP_OUT/FAILED_SAMPLES.txt"
  echo "=============================================================="
  echo "[WARN] ${#FAILED[@]} sample(s) FAILED fastp:"
  printf '   %s\n' "${FAILED[@]}"
  echo "       full list: $FASTP_OUT/FAILED_SAMPLES.txt"
  echo "=============================================================="
fi
echo "[$(date '+%F %T')] >>> fastp complete. Cleaned pairs: ${n_ok} / ${TOTAL_JOBS}"
(( n_ok == 0 )) && { echo "[ERROR] No cleaned reads produced; halting pipeline." >&2; exit 1; }

echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
exit 0
