#!/bin/bash
#SBATCH -J cn02_mosdepth
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -N 1
#SBATCH --mem=200G
#SBATCH --cpus-per-task=48
#SBATCH -o 02_cn_mosdepth.out
# =============================================================================
# Stage 02: mosdepth per-sample depth for CNV / aneuploidy.
# Consumes BAMs from stage 01. Key outputs used downstream by fig6.r:
#   $MOS_OUT/mosdepth/{sample}.mosdepth.summary.txt   (per-chrom mean depth)
#   $MOS_OUT/mosdepth/{sample}.regions.bed.gz         ({MOS_WINDOW_SIZE}-bp bins)
# =============================================================================

PIPELINE_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
source "${PIPELINE_DIR}/00_config.sh"
apply_strict_mode

THREADS_TOTAL="${SLURM_CPUS_PER_TASK:-4}"
MAP_DIR="$MAP_OUT"
FAIL_DIR="$MOS_OUT/.failed"
mkdir -p "$MOS_OUT"/{mosdepth,logs,tmp} "$FAIL_DIR"
rm -f "$FAIL_DIR"/* 2>/dev/null || true

cd "$MOS_OUT"
echo "Working Directory = $(pwd)"
activate_conda "$CONDA_ENV_MOSDEPTH"
print_config_banner
echo "[$(date '+%F %T')] >>> MAP_DIR     = $MAP_DIR"
echo "[$(date '+%F %T')] >>> WINDOW_SIZE = $MOS_WINDOW_SIZE"
echo "[$(date '+%F %T')] >>> THRESHOLDS  = ${MOS_THRESHOLDS:-<none>}"
echo "[$(date '+%F %T')] >>> FAST_MODE   = $MOS_FAST_MODE  NO_PER_BASE = $MOS_NO_PER_BASE"

command -v mosdepth >/dev/null 2>&1 || { echo "[ERROR] mosdepth not in PATH" >&2; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "[ERROR] samtools not in PATH" >&2; exit 1; }

# ---- Discover BAMs ----
BAMS=( "$MAP_DIR"/*.bam )
(( ${#BAMS[@]} == 0 )) && { echo "[ERROR] No BAMs under $MAP_DIR" >&2; exit 1; }
SAMPLES=(); BAM_LIST=()
for bam in "${BAMS[@]}"; do
  sample="$(basename "$bam" .bam)"
  [[ -s "$bam" ]]        || { echo "[WARN] empty BAM, skip: $bam"; continue; }
  [[ -s "${bam}.bai" ]]  || { echo "[WARN] missing .bai, skip: $bam"; continue; }
  SAMPLES+=( "$sample" ); BAM_LIST+=( "$bam" )
done
TOTAL_JOBS="${#SAMPLES[@]}"
echo "[$(date '+%F %T')] >>> Detected ${TOTAL_JOBS} BAMs for mosdepth."
(( TOTAL_JOBS == 0 )) && { echo "[ERROR] No usable BAM/BAI pairs." >&2; exit 1; }

# ---- Thread budget ----
MAX_JOBS="$MOS_MAX_JOBS"
(( MAX_JOBS > THREADS_TOTAL )) && MAX_JOBS="$THREADS_TOTAL"
(( MAX_JOBS > TOTAL_JOBS ))    && MAX_JOBS="$TOTAL_JOBS"
MOS_THREADS=$(( THREADS_TOTAL / MAX_JOBS )); (( MOS_THREADS < 1 )) && MOS_THREADS=1
echo "[$(date '+%F %T')] >>> parallel jobs: $MAX_JOBS  (mosdepth -t $MOS_THREADS)"

# ---- FIFO semaphore ----
SEM_FIFO="$(mktemp -u)"; mkfifo "$SEM_FIFO"; exec 3<>"$SEM_FIFO"; rm -f "$SEM_FIFO"
for ((i=0;i<MAX_JOBS;i++)); do echo >&3; done
run_parallel() { read -u 3; ( trap 'echo >&3' EXIT; "$@" ) & }

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

# ---- mosdepth task ----
mos_job() {
  local sample="$1" bam="$2"
  local prefix="$MOS_OUT/mosdepth/${sample}"
  local log="$MOS_OUT/logs/${sample}.log"; : > "$log"
  local tmp_dir="$MOS_OUT/tmp/$sample"; mkdir -p "$tmp_dir"
  export TMPDIR="$tmp_dir"

  echo "[$(date '+%F %T')] >>> START mosdepth: $sample" | tee -a "$log"

  local mos_args=( -t "$MOS_THREADS" --by "$MOS_WINDOW_SIZE" )
  [[ -n "$MOS_THRESHOLDS" ]]    && mos_args+=( --thresholds "$MOS_THRESHOLDS" )
  (( MOS_FAST_MODE   == 1 ))    && mos_args+=( --fast-mode )
  (( MOS_NO_PER_BASE == 1 ))    && mos_args+=( --no-per-base )

  if ! mosdepth "${mos_args[@]}" "$prefix" "$bam" >>"$log" 2>&1; then
    echo "[ERROR] mosdepth failed: $sample" | tee -a "$log"
    touch "$FAIL_DIR/$sample"; rm -rf "$tmp_dir"; update_progress "$sample"; return 1
  fi
  if [[ ! -s "${prefix}.mosdepth.summary.txt" || ! -s "${prefix}.regions.bed.gz" ]]; then
    echo "[ERROR] missing summary/regions output: $sample" | tee -a "$log"
    touch "$FAIL_DIR/$sample"; rm -rf "$tmp_dir"; update_progress "$sample"; return 1
  fi

  rm -rf "$tmp_dir"
  echo "[$(date '+%F %T')] >>> DONE  mosdepth: $sample" | tee -a "$log"
  update_progress "$sample"
}

# ---- Submit ----
for i in "${!SAMPLES[@]}"; do
  run_parallel mos_job "${SAMPLES[$i]}" "${BAM_LIST[$i]}"
done
wait || true
exec 3>&- 3<&- || true

# ---- Failure report ----
mapfile -t FAILED < <(find "$FAIL_DIR" -maxdepth 1 -type f -printf '%f\n' 2>/dev/null | sort)
n_ok=$(find "$MOS_OUT/mosdepth" -maxdepth 1 -name '*.mosdepth.summary.txt' 2>/dev/null | wc -l)
if (( ${#FAILED[@]} > 0 )); then
  printf '%s\n' "${FAILED[@]}" > "$MOS_OUT/FAILED_SAMPLES.txt"
  echo "=============================================================="
  echo "[WARN] ${#FAILED[@]} sample(s) FAILED mosdepth:"
  printf '   %s\n' "${FAILED[@]}"
  echo "       full list: $MOS_OUT/FAILED_SAMPLES.txt"
  echo "=============================================================="
fi
echo "[$(date '+%F %T')] >>> mosdepth complete. Summaries: ${n_ok} / ${TOTAL_JOBS}"
echo "Outputs under: $MOS_OUT/mosdepth/"
echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
exit 0
