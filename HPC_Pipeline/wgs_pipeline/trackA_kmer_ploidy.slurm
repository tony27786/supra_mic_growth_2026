#!/bin/bash
#SBATCH -J trackA_kmer
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -N 1
#SBATCH --cpus-per-task=48
#SBATCH --mem=240G
#SBATCH -o trackA_kmer_ploidy.out
# =============================================================================
# Track A -- reference-free ploidy / heterozygosity ("what is Cn160").
#   Haploid (het ~0%)  -> ~18000 het vs H99 = reference-divergence artifact (not
#                         real het, not a mixture).
#   Diploid + high het  -> real heterozygosity (AD hybrid / diploid) -> then
#                         H99-vs-JEC21 to name the A/D haplotypes.
#
#   ALL-FastK (one DB per sample, k=21):
#     FastK -> Histex -G -> GenomeScope2   DECISIVE (heterozygosity %)
#     FastK -> smudgeplot                  2D ploidy visual (best-effort)
#
#   Robustness (v2, after first run):
#     * NO `set -u` (empty optional arg must not abort).
#     * Per-sample isolation: a sample that fails is logged + recorded in
#       FAILED_SAMPLES.txt; the run CONTINUES (one bad sample never aborts all).
#     * GenomeScope called WITHOUT -n, so it writes summary.txt.
#     * FastK uses a per-sample temp dir (-P) + cleanup.
#     * On any step failure, the tail of its log is printed to this main log.
#
#   Parallel (v3, this version):
#     * Bash FIFO semaphore (token on FD 3), same template as extract_het.
#     * MAX_JOBS samples run concurrently; THREADS_PER_JOB = ceil(CPUS / MAX_JOBS).
#       48 cores / 8 jobs = 6 threads per sample (FastK -T / smudgeplot -t).
#     * run_parallel exports THREADS per job, so process_one needs no change.
#     * `sleep 1` between launches: stagger FastK temp-dir creation / IO bursts
#       instead of firing all 8 at the same instant.
#     * SUMMARY / FAILED_SAMPLES appends are single short lines (atomic under
#       PIPE_BUF); summary is re-sorted by sample after `wait`.
#
#   ! MEMORY: FastK -M is PER JOB. Total FastK target ~= MAX_JOBS * MEM_GB.
#     Default 8 * 16 = 128G under the 160G #SBATCH --mem cap. If you raise
#     MAX_JOBS or MEM_GB, raise #SBATCH --mem too (and confirm the rna
#     partition node actually has 48 cores).
#
#   INPUT IS A FOLDER. Auto-detects either layout:
#     flat:           <dir>/<sample>_R1.fq.gz   <dir>/<sample>_R2.fq.gz
#     per-sample dir: <dir>/<sample>/*_R1*.fq.gz  *_R2*.fq.gz
#
#   Run:
#     sbatch trackA_kmer_ploidy.slurm <input_dir> [sample1,sample2,...]
#   optional env: OUT_DIR=... MEM_GB=16 K=21 MAX_JOBS=8 RUN_SMUDGE=0 HET_VERDICT_PCT=1.0
# =============================================================================

# ---- setup: fail fast on config errors (note: NO -u, per run-1 feedback) ----
set -eo pipefail
trap 'echo "[ERROR] setup failed at line $LINENO" >&2' ERR
export PYTHONUNBUFFERED=1
shopt -s nullglob

INDIR="${1:?usage: sbatch trackA_kmer_ploidy.slurm <input_dir> [sample1,sample2,...]}"
FILTER="${2:-}"
OUT_ROOT="${OUT_DIR:-$PWD/trackA_kmer_$(date +%y%m%d)}"
K="${K:-21}"
FASTK_T="${FASTK_T:-4}"           # FastK table threshold: keep k-mers with count >= this.
                                  # .hist (-> GenomeScope) is ALWAYS complete regardless of -t.
RUN_SMUDGE="${RUN_SMUDGE:-0}"     # 0 (default) = SKIP the whole smudgeplot/hetmers leg.
                                  #   hetmers (PloidyPlot) is the slow, NON-decisive leg;
                                  #   GenomeScope het% above is decisive. smudge_ploidy stays NA
                                  #   (downstream already tolerates NA, same as a smudge failure).
                                  #   The FastK DB is kept on disk, so you can build smudge plots
                                  #   later (per-strain, many cores) without recomputing FastK.
                                  # 1 = run hetmers inline (then SMUDGE_PLOT decides plotting).
SMUDGE_PLOT="${SMUDGE_PLOT:-0}"   # only used when RUN_SMUDGE=1. 1 = also run `smudgeplot all` (needs pandas<3).
HET_VERDICT_PCT="${HET_VERDICT_PCT:-1.0}"
[[ -d "$INDIR" ]] || { echo "[ERROR] not a directory: $INDIR" >&2; exit 1; }

# ------------------------------
# Parallel Settings (FIFO) -- same template as extract_het
# ------------------------------
THREADS_TOTAL="${SLURM_CPUS_PER_TASK:-48}"
MAX_JOBS="${MAX_JOBS:-16}"
THREADS_PER_JOB=$(( (THREADS_TOTAL + MAX_JOBS - 1) / MAX_JOBS ))   # ceil(48/16)=3
(( THREADS_PER_JOB < 1 )) && THREADS_PER_JOB=1
MEM_GB="${MEM_GB:-16}"            # PER JOB FastK memory cap (-M). Total ~= MAX_JOBS*MEM_GB.

CONDA_LAUNCHER="${CONDA_LAUNCHER:-${CONDA_LAUNCHER}}"
# shellcheck source=/dev/null
source "$("$CONDA_LAUNCHER" info --base)/etc/profile.d/conda.sh"
conda activate hetero_kmer
echo "[$(date '+%F %T')] env=$CONDA_DEFAULT_ENV cpus=$THREADS_TOTAL max_jobs=$MAX_JOBS threads/job=$THREADS_PER_JOB mem/job=${MEM_GB}G k=$K run_smudge=$RUN_SMUDGE"

for t in FastK Histex genomescope.R smudgeplot; do
  command -v "$t" >/dev/null 2>&1 || { echo "[ERROR] '$t' not on PATH in env hetero_kmer." >&2;
    echo "        FastK/Histex come from: mamba install -c bioconda fastk" >&2; exit 1; }
done

# ---- discover samples (both layouts) ----
declare -A R1MAP R2MAP
for d in "$INDIR"/*/; do
  s="$(basename "${d%/}")"
  c1=( "$d"*_R1*.f*q.gz "$d"*_1.f*q.gz ); c2=( "$d"*_R2*.f*q.gz "$d"*_2.f*q.gz )
  if [[ ${#c1[@]} -ge 1 && ${#c2[@]} -ge 1 ]]; then R1MAP["$s"]="${c1[0]}"; R2MAP["$s"]="${c2[0]}"; fi
done
for r1 in "$INDIR"/*_R1.f*q.gz; do
  [[ -f "$r1" ]] || continue
  bn="$(basename "$r1")"; s="${bn%%_R1*}"; r2="${r1/_R1/_R2}"
  [[ -f "$r2" ]] && { R1MAP["$s"]="$r1"; R2MAP["$s"]="$r2"; }
done
(( ${#R1MAP[@]} > 0 )) || { echo "[ERROR] no paired reads under $INDIR" >&2; exit 1; }

declare -A KEEP; USE_FILTER=0
if [[ -n "$FILTER" ]]; then
  USE_FILTER=1; IFS=',' read -ra arr <<< "$FILTER"
  for x in "${arr[@]}"; do KEEP["$(echo "$x" | xargs)"]=1; done
fi

mkdir -p "$OUT_ROOT"
SUMMARY="$OUT_ROOT/trackA_ploidy_summary.tsv"
FAILED_LIST="$OUT_ROOT/FAILED_SAMPLES.txt"
printf "sample\tgs_heterozygous_ab\tgs_model_fit\tgs_haploid_len\tsmudge_ploidy\tverdict_hint\n" > "$SUMMARY"
: > "$FAILED_LIST"

mapfile -t SAMPLES < <(printf "%s\n" "${!R1MAP[@]}" | sort -V)

# resolve the effective run list (apply filter now) so the progress bar denominator is correct
RUN_SAMPLES=()
for S in "${SAMPLES[@]}"; do
  (( USE_FILTER )) && [[ -z "${KEEP[$S]:-}" ]] && continue
  RUN_SAMPLES+=("$S")
done
TOTAL_JOBS="${#RUN_SAMPLES[@]}"
echo "[$(date '+%F %T')] discovered ${#SAMPLES[@]} sample(s); will run ${TOTAL_JOBS}; filter=${FILTER:-<none>}"
(( TOTAL_JOBS == 0 )) && { echo "[ERROR] no samples left after filter '${FILTER}'" >&2; exit 1; }

# ---- progress bar (flock-locked shared counter, same as 00_fastp) ----
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

# print last N lines of a log, indented, if it exists
tail_log() { [[ -f "$1" ]] && tail -n "${2:-6}" "$1" | sed 's/^/      | /'; }

# ---- per-sample work (returns 1 only if NO usable GenomeScope result) -------
process_one() {
  local S="$1" R1="$2" R2="$3"
  local SDIR="$OUT_ROOT/$S" TMP="$OUT_ROOT/$S/tmp"
  mkdir -p "$SDIR" "$TMP"
  local FK="$SDIR/${S}_fastk"
  rm -f "${FK}".* "$SDIR/.${S}_fastk".* 2>/dev/null
  rm -rf "${TMP:?}/"* 2>/dev/null

  # ---- FastK (feeds both legs); -P keeps temp local; -T/-M are per-job ----
  echo "[$(date '+%F %T')] $S: FastK k=$K (T=$THREADS M=${MEM_GB}G) ..."
  if ! FastK -v -t"$FASTK_T" -k"$K" -M"$MEM_GB" -T"$THREADS" -P"$TMP" "$R1" "$R2" -N"$FK" \
        > "$SDIR/${S}_fastk.log" 2>&1; then
    echo "[ERROR] $S: FastK failed -- log tail:"; tail_log "$SDIR/${S}_fastk.log"
    return 1
  fi

  # ---- Histex -> GenomeScope-format histogram ----
  local HIST="$SDIR/${S}_k${K}.hist"
  if ! Histex -G "$FK" > "$HIST" 2> "$SDIR/${S}_histex.log" || [[ ! -s "$HIST" ]]; then
    echo "[ERROR] $S: Histex failed/empty -- log tail:"; tail_log "$SDIR/${S}_histex.log"
    return 1
  fi

  # ---- leg 1: GenomeScope2 (decisive). NO -n -> writes summary.txt ----------
  echo "[$(date '+%F %T')] $S: GenomeScope2 (p2 + p1) ..."
  genomescope.R -i "$HIST" -o "$SDIR/genomescope_p2" -k "$K" -p 2 > "$SDIR/${S}_gs_p2.log" 2>&1 \
    || echo "[WARN] $S: GenomeScope p2 returned nonzero"
  genomescope.R -i "$HIST" -o "$SDIR/genomescope_p1" -k "$K" -p 1 > "$SDIR/${S}_gs_p1.log" 2>&1 \
    || echo "[WARN] $S: GenomeScope p1 returned nonzero"

  local GS het="NA" fit="NA" hlen="NA"
  GS="$(ls "$SDIR/genomescope_p2/"*summary.txt 2>/dev/null | head -1)"
  if [[ -z "$GS" ]]; then
    echo "[ERROR] $S: GenomeScope produced no summary -- p2 log tail:"; tail_log "$SDIR/${S}_gs_p2.log"
    return 1
  fi
  het="$(grep -iE "Heterozygous|Heterozygosity" "$GS" | head -1 | grep -oE "[0-9.]+ *%" | head -1 | tr -d ' ')"
  fit="$(grep -i  "Model Fit"      "$GS" | head -1 | grep -oE "[0-9.]+ *%"  | head -1 | tr -d ' ')"
  hlen="$(grep -iE "Haploid Length" "$GS" | head -1 | grep -oE "[0-9,]+ *bp" | head -1 | tr -d ' ')"
  : "${het:=NA}"; : "${fit:=NA}"; : "${hlen:=NA}"

  # ---- leg 2: smudgeplot (same FastK DB; best-effort, never fatal) ----------
  # RUN_SMUDGE=0 (default) skips this entirely -- it is the slow, NON-decisive
  # leg. smud stays "NA"; FastK DB at ${FK}.* is kept for a later per-strain run.
  local smud="NA"
  if [[ "$RUN_SMUDGE" == "1" ]]; then
    local L SMU
    L="$(smudgeplot cutoff "$HIST" L 2>/dev/null | tr -dc '0-9')"; [[ -z "$L" ]] && L=10
    echo "[$(date '+%F %T')] $S: smudgeplot hetmers -L $L ..."
    if smudgeplot hetmers -L "$L" -t "$THREADS" -o "$SDIR/${S}_pairs" --verbose "$FK" \
          > "$SDIR/${S}_smudge.log" 2>&1; then
      SMU="$(ls "$SDIR/${S}_pairs"*.smu 2>/dev/null | head -1)"
      if [[ "$SMUDGE_PLOT" != "1" ]]; then
        echo "[$(date '+%F %T')] $S: .smu built ($SMU); SMUDGE_PLOT=0 -> skipping plot (make it later in a pandas<3 env: smudgeplot all -o ${S}_smudge $SMU)"
      elif [[ -n "$SMU" ]] && smudgeplot all -o "$SDIR/${S}_smudge" "$SMU" >> "$SDIR/${S}_smudge.log" 2>&1; then
        sp="$(grep -rhiE "ploidy" "$SDIR/${S}_smudge"*summary* 2>/dev/null | head -1 | grep -oE "[0-9]+" | head -1)"
        [[ -n "${sp:-}" ]] && smud="$sp"
      else
        echo "[WARN] $S: smudgeplot plot step failed (likely pandas>=3; .smu is fine, re-plot in a pandas<3 env) -- log tail:"; tail_log "$SDIR/${S}_smudge.log"
      fi
    else
      echo "[WARN] $S: smudgeplot hetmers failed -- log tail:"; tail_log "$SDIR/${S}_smudge.log"
    fi
  else
    echo "[$(date '+%F %T')] $S: RUN_SMUDGE=0 -> skip smudgeplot/hetmers (non-decisive; FastK DB kept at ${FK}.*)"
  fi

  # ---- verdict hint (heuristic; summary.txt + plots are authoritative) ------
  local verdict="check_plot" hnum
  hnum="$(echo "${het}" | grep -oE "[0-9.]+" | head -1)"
  [[ -n "$hnum" ]] && verdict="$(awk -v h="$hnum" -v thr="$HET_VERDICT_PCT" 'BEGIN{print (h+0>=thr)?"diploid/het candidate":"haploid-like (het~0)"}')"

  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$S" "$het" "$fit" "$hlen" "$smud" "$verdict" >> "$SUMMARY"
  echo "[$(date '+%F %T')] $S done: het(ab)=$het fit=$fit -> $verdict"
  return 0
}

# ---- wrapper run inside each FIFO worker: isolates failure per sample -------
run_one_sample() {
  local S="$1"
  echo "=============================================================="
  echo "[$(date '+%F %T')] >>> SAMPLE = $S  (T=$THREADS)"
  echo "    R1=${R1MAP[$S]}"; echo "    R2=${R2MAP[$S]}"
  if ! process_one "$S" "${R1MAP[$S]}" "${R2MAP[$S]}"; then
    echo "[FAIL] $S recorded in FAILED_SAMPLES.txt; continuing."
    echo "$S" >> "$FAILED_LIST"
    printf "%s\tFAILED\tFAILED\tFAILED\tNA\tFAILED\n" "$S" >> "$SUMMARY"
  fi
  update_progress "$S"      # count once per sample, whether it succeeded or failed
}

# ---- FIFO semaphore (tokens on FD 3) ---------------------------------------
SEM_FIFO="$(mktemp -u)"
mkfifo "$SEM_FIFO"
exec 3<>"$SEM_FIFO"
rm -f "$SEM_FIFO"
for ((i=0;i<MAX_JOBS;i++)); do echo >&3; done

run_parallel() {
  read -u 3                                   # acquire token (blocks if none free)
  (
    trap 'echo >&3' EXIT                       # always return the token
    export OMP_NUM_THREADS="$THREADS_PER_JOB"
    export THREADS="$THREADS_PER_JOB"          # process_one reads $THREADS for FastK -T / smudgeplot -t
    "$@"
  ) &
}

# ---- submit loop: stagger by 1s, never abort the whole run -----------------
trap - ERR
set +e
for S in "${RUN_SAMPLES[@]}"; do
  run_parallel run_one_sample "$S"
  sleep 1                                      # wait 1s before submitting the next
done

wait
exec 3>&- 3<&-
rm -f "$COUNTER_FILE" "$LOCK_FILE" 2>/dev/null

# ---- re-sort summary body by sample (completion order is nondeterministic) --
if [[ -s "$SUMMARY" ]]; then
  { head -1 "$SUMMARY"; tail -n +2 "$SUMMARY" | sort -V; } > "$SUMMARY.tmp" && mv "$SUMMARY.tmp" "$SUMMARY"
fi

echo "=============================================================="
echo ">>> combined verdict table: $SUMMARY"
echo "--------------------------------------------------------------"
column -t -s$'\t' "$SUMMARY" 2>/dev/null || cat "$SUMMARY"
nfail=$(grep -c . "$FAILED_LIST" 2>/dev/null || echo 0)
(( nfail > 0 )) && { echo "--------------------------------------------------------------"; echo ">>> $nfail sample(s) FAILED (see $FAILED_LIST):"; sort -V "$FAILED_LIST" | sed 's/^/    /'; }
echo "=============================================================="
echo ">>> Authoritative per-sample: <out>/<sample>/genomescope_p2/{summary.txt,linear_plot.png,log_plot.png}"
echo ">>> het(ab) ~0%  -> HAPLOID (het vs H99 = divergence artifact, not real, not a mixture)."
echo ">>> het(ab) several % + clear het peak -> DIPLOID/HETEROZYGOUS -> then H99-vs-JEC21."
echo ">>> Cn160 vs Cn160m should give the SAME verdict."
echo ">>> SLURM BATCH TERMINATED AT [$(date '+%F %T')]"
exit 0
