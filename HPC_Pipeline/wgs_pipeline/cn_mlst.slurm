#!/bin/bash
#SBATCH -J cn_all_mlst
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -N 1
#SBATCH --mem=50G
#SBATCH --cpus-per-task=10
#SBATCH -o cn_mlst.out
# =============================================================================
# In-silico MLST (stringMLST, k-mer based -> runs directly from fastp reads,
# NO assembly needed). Optimised to the project's established pattern:
#   FIFO semaphore | per-sample .ok success marker + failure isolation
#   | merge ONLY .ok samples | n_ok/N_total report | STRICT toggle | resume.
# Run:  sbatch cn_mlst.slurm            (lenient: finish, report failures)
#       STRICT=1 sbatch cn_mlst.slurm   (exit 1 if any sample failed)
# =============================================================================

cd "$HOME/working/hetero_fig6_wf" || exit 1
echo "Working Directory = $(pwd)"

CONDA_BASE="$(${CONDA_LAUNCHER} info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate stringmlst
echo "CONDA_PREFIX=$CONDA_PREFIX"

set -uo pipefail                       # NOTE: no -e; failures are handled per-sample, not fatally
trap 'code=$?; echo -e "\n[WARN] shell exit $code at line $LINENO" >&2' ERR
IFS=$'\n\t'
shopt -s nullglob
export PYTHONUNBUFFERED=1

# ---- config -----------------------------------------------------------------
INPUT_ROOT="$LUSTRE/working/cn_80_fastpout_260702"
DB_PREFIX="$HOME/mlstdbs/stringmlst/CN"
OUT_DIR="$HOME/working/cn_80_stringmlst_$(date +%y%m%d)"
FINAL_REPORT="$OUT_DIR/Cn_all_merged_results.txt"
FAILED_LIST="$OUT_DIR/FAILED_SAMPLES.txt"
MAX_JOBS=10
STRICT="${STRICT:-0}"
KMER=35

mkdir -p "$OUT_DIR"/{results,logs,tmp}
: > "$FAILED_LIST"

echo "[$(date '+%F %T')] INPUT_ROOT=$INPUT_ROOT | DB_PREFIX=$DB_PREFIX | OUT_DIR=$OUT_DIR | MAX_JOBS=$MAX_JOBS | STRICT=$STRICT"
command -v stringMLST.py >/dev/null 2>&1 || { echo "[ERROR] stringMLST.py not found" >&2; exit 1; }
for ext in "$DB_PREFIX"_*; do :; done
[[ -e "${DB_PREFIX}_profile.txt" || -n "$(echo "$DB_PREFIX"_* 2>/dev/null)" ]] || echo "[WARN] DB files for prefix $DB_PREFIX not obviously present"

# ---- FIFO semaphore (token on FD 3) ----------------------------------------
SEM_FIFO="$(mktemp -u)"; mkfifo "$SEM_FIFO"; exec 3<>"$SEM_FIFO"; rm -f "$SEM_FIFO"
for ((i=0;i<MAX_JOBS;i++)); do echo >&3; done
run_parallel() { read -u 3; ( trap 'echo >&3' EXIT; "$@" ) & }

# ---- discover samples -------------------------------------------------------
SAMPLES=(); DIR_LIST=()
for d in "$INPUT_ROOT"/Cn*; do
  [[ -d "$d" ]] || continue
  sn="$(basename "$d")"
  if compgen -G "$d/*.fq.gz" >/dev/null; then SAMPLES+=( "$sn" ); DIR_LIST+=( "$d" )
  else echo "[WARN] no .fq.gz in $d, skipping"; fi
done
TOTAL_JOBS="${#SAMPLES[@]}"
echo "[$(date '+%F %T')] detected ${TOTAL_JOBS} samples"
(( TOTAL_JOBS == 0 )) && { echo "[ERROR] no samples found" >&2; exit 1; }

# ---- progress bar (flock-guarded counter) ----------------------------------
PROGRESS_WIDTH=45
COUNTER_FILE="$(mktemp)"; LOCK_FILE="$(mktemp)"; echo 0 > "$COUNTER_FILE"
update_progress() {
  local sample="$1"
  { flock -x 9
    local done_jobs; done_jobs="$(<"$COUNTER_FILE")"; done_jobs=$((done_jobs + 1)); echo "$done_jobs" > "$COUNTER_FILE"
    local percent=$(( done_jobs * 100 / TOTAL_JOBS )) bar_len=$(( done_jobs * PROGRESS_WIDTH / TOTAL_JOBS )) bar="" i
    for ((i=0;i<bar_len;i++)); do bar+="#"; done; for ((i=bar_len;i<PROGRESS_WIDTH;i++)); do bar+="."; done
    printf "[PROGRESS] [%s] %3d%% (%d/%d) latest: %s\n" "$bar" "$percent" "$done_jobs" "$TOTAL_JOBS" "$sample"
  } 9>>"$LOCK_FILE"
}

# ---- per-sample job: explicit success (.ok) / failure isolation ------------
mlst_job() {
  local sample="$1" dir="$2"
  local out="$OUT_DIR/results/${sample}_mlst.tsv" ok="$OUT_DIR/results/${sample}.ok" log="$OUT_DIR/logs/${sample}.log"
  local r1="${dir}/${sample}_clean_R1.fq.gz" r2="${dir}/${sample}_clean_R2.fq.gz"

  if [[ -s "$out" && -f "$ok" ]]; then update_progress "$sample (SKIP)"; return 0; fi   # resume
  if [[ ! -f "$r1" || ! -f "$r2" ]]; then
    echo "[ERROR] missing FASTQ for $sample ($r1 / $r2)" | tee -a "$log"
    { flock -x 8; echo "$sample	missing_fastq" >> "$FAILED_LIST"; } 8>>"$LOCK_FILE"
    update_progress "$sample (FAILED)"; return 0
  fi

  rm -f "$ok"
  echo "[$(date '+%F %T')] start: $sample" | tee -a "$log"
  if stringMLST.py --predict -1 "$r1" -2 "$r2" -P "$DB_PREFIX" -k "$KMER" -o "$out" >>"$log" 2>&1 && [[ -s "$out" ]]; then
    touch "$ok"; echo "[$(date '+%F %T')] done:  $sample" | tee -a "$log"; update_progress "$sample"
  else
    echo "[ERROR] stringMLST failed: $sample" | tee -a "$log"
    { flock -x 8; echo "$sample	stringmlst_error" >> "$FAILED_LIST"; } 8>>"$LOCK_FILE"
    update_progress "$sample (FAILED)"
  fi
  return 0                              # always 0 so the subshell releases its token cleanly
}

# ---- submit -----------------------------------------------------------------
for i in "${!SAMPLES[@]}"; do run_parallel mlst_job "${SAMPLES[$i]}" "${DIR_LIST[$i]}"; done
wait
exec 3>&- 3<&-
echo "[$(date '+%F %T')] all jobs returned"

# ---- merge ONLY successful (.ok) samples -----------------------------------
FILES_TO_MERGE=()
for f in "$OUT_DIR"/results/*_mlst.tsv; do
  s="$(basename "$f" _mlst.tsv)"
  [[ -f "$OUT_DIR/results/${s}.ok" && -s "$f" ]] && FILES_TO_MERGE+=( "$f" )
done
if ((${#FILES_TO_MERGE[@]})); then
  mapfile -t FILES_TO_MERGE < <(printf '%s\n' "${FILES_TO_MERGE[@]}" | sort -V)
  head -n 1 "${FILES_TO_MERGE[0]}" > "$FINAL_REPORT"
  for f in "${FILES_TO_MERGE[@]}"; do tail -n +2 "$f" >> "$FINAL_REPORT"; done
fi

n_ok=${#FILES_TO_MERGE[@]}
n_fail=$( [[ -s "$FAILED_LIST" ]] && wc -l < "$FAILED_LIST" || echo 0 )
echo "----------------------------------------------------------------"
echo "[SUMMARY] merged n_ok=$n_ok / total=$TOTAL_JOBS ; failed=$n_fail -> $FINAL_REPORT"
(( n_fail > 0 )) && { echo "[SUMMARY] failed samples:"; cat "$FAILED_LIST"; }
[[ -s "$FINAL_REPORT" ]] && { echo "--- preview ---"; head -n 6 "$FINAL_REPORT"; }
rm -f "$COUNTER_FILE" "$LOCK_FILE"

if [[ "$STRICT" == "1" && "$n_fail" -gt 0 ]]; then
  echo "[STRICT] $n_fail sample(s) failed -> exit 1"; exit 1
fi
echo ">>> DONE  [$(date '+%F %T')]"
