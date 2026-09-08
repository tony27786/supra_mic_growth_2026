#!/bin/bash
#SBATCH -J cn_ksnp4
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -N 1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=48
#SBATCH -o cn_ksnp4.out
# =============================================================================
# Reference-free SNP tree of the 20 taxa (19 working parents + H99) with kSNP4.
#   Rationale: reference-free (no H99 bias) -> consistent with the paper's thesis;
#   SNP-level with support -> stronger than a k-mer distance (Mash) tree; handles
#   the divergent (Cn160/150) and het-diploid (Cn147) strains without a reference.
#   At 200x coverage raw reads are noisy for kSNP -> ASSEMBLE first (shovill), then
#   kSNP4 on clean contigs.  k chosen by Kchooser4 (fallback 21, matches FastK).
#
#   Stage 1  shovill assemble  (FIFO semaphore | .ok isolation)
#   Stage 2  build kSNP4 infile from .ok assemblies
#   Stage 3  Kchooser4 -> optimum k
#   Stage 4  kSNP4 -core -ML  -> Newick trees
#
# Run:  sbatch cn_ksnp4.slurm      /      STRICT=1 sbatch cn_ksnp4.slurm
# VERIFY: conda env names (ASM_ENV/KSNP_ENV); kSNP4/Kchooser4/shovill CLI flags vs your build.
# =============================================================================
ASM_ENV="${ASM_ENV:-ksnp}"
KSNP4_HOME="${KSNP4_HOME:-$HOME/kSNP4.1pkg}"   # official kSNP4.1 package
K_OVERRIDE="${K_OVERRIDE:-}"           # set e.g. 21 to skip Kchooser4

mkdir -p "$HOME/working/cn_ksnp4" && cd "$HOME/working/cn_ksnp4" || exit 1
CONDA_BASE="$(${CONDA_LAUNCHER} info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"

set -uo pipefail
trap 'code=$?; echo -e "\n[WARN] shell exit $code at line $LINENO" >&2' ERR
IFS=$'\n\t'; shopt -s nullglob; export PYTHONUNBUFFERED=1

# ---- config -----------------------------------------------------------------
INPUT_ROOT="$LUSTRE/working/cn_80_fastpout_260702"
OUT_DIR="$LUSTRE/working/cn_ksnp4_$(date +%y%m%d)"
ASM_DIR="$OUT_DIR/assemblies"; FASTA_DIR="$OUT_DIR/fasta"; KOUT="$OUT_DIR/ksnp4_out"
INFILE="$OUT_DIR/kSNP4_infile"; FAILED_LIST="$OUT_DIR/FAILED_SAMPLES.txt"; LOCK_FILE="$(mktemp)"
mkdir -p "$ASM_DIR" "$FASTA_DIR" "$OUT_DIR/logs"; : > "$FAILED_LIST"
K_DEFAULT=21; MAX_JOBS=8; ASM_CPUS=6; STRICT="${STRICT:-0}"
TOTAL_CPUS="${SLURM_CPUS_PER_TASK:-16}"

TAXA=(Cn101 Cn106 Cn126 Cn128 Cn160 Cn162 Cn208 Cn213 Cn215 Cn217 \
      Cn147 Cn150 Cn169 Cn174 Cn225 Cn226 Cn196 Cn197 Cn198 CnH99)
TOTAL="${#TAXA[@]}"

# ---- FIFO semaphore ---------------------------------------------------------
SEM_FIFO="$(mktemp -u)"; mkfifo "$SEM_FIFO"; exec 3<>"$SEM_FIFO"; rm -f "$SEM_FIFO"
for ((i=0;i<MAX_JOBS;i++)); do echo >&3; done
run_parallel() { read -u 3; ( trap 'echo >&3' EXIT; "$@" ) & }

# ============================ Stage 1: assemble ==============================
conda activate "$ASM_ENV"; echo "[ASM] env=$CONDA_PREFIX"
command -v shovill >/dev/null 2>&1 || { echo "[ERROR] shovill not found in '$ASM_ENV'"; exit 1; }

assemble_job() {
  local s="$1" dir="$INPUT_ROOT/$1"
  local r1="$dir/${s}_clean_R1.fq.gz" r2="$dir/${s}_clean_R2.fq.gz"
  local fa="$FASTA_DIR/${s}.fasta" ok="$FASTA_DIR/${s}.ok" log="$OUT_DIR/logs/${s}.asm.log"
  if [[ -s "$fa" && -f "$ok" ]]; then echo "[SKIP] $s"; return 0; fi
  if [[ ! -f "$r1" || ! -f "$r2" ]]; then
    echo "[ERROR] missing FASTQ: $s" | tee -a "$log"
    { flock -x 8; echo -e "$s\tmissing_fastq" >> "$FAILED_LIST"; } 8>>"$LOCK_FILE"; return 0
  fi
  rm -f "$ok"
  if shovill --R1 "$r1" --R2 "$r2" --outdir "$ASM_DIR/$s" --cpus "$ASM_CPUS" --force >>"$log" 2>&1 \
       && [[ -s "$ASM_DIR/$s/contigs.fa" ]]; then
    cp -f "$ASM_DIR/$s/contigs.fa" "$fa"; touch "$ok"; echo "[OK] $s"
  else
    echo "[ERROR] shovill failed: $s" | tee -a "$log"
    { flock -x 8; echo -e "$s\tassembly_error" >> "$FAILED_LIST"; } 8>>"$LOCK_FILE"
  fi
  return 0
}
for s in "${TAXA[@]}"; do run_parallel assemble_job "$s"; done
wait; exec 3>&- 3<&-
echo "[$(date '+%F %T')] assemblies returned"

# ============================ Stage 2: infile ================================
: > "$INFILE"
n_ok=0
for s in "${TAXA[@]}"; do
  fa="$FASTA_DIR/${s}.fasta"
  if [[ -f "$FASTA_DIR/${s}.ok" && -s "$fa" ]]; then
    printf '%s\t%s\n' "$(readlink -f "$fa")" "$s" >> "$INFILE"; n_ok=$((n_ok+1))
  fi
done
n_fail=$( [[ -s "$FAILED_LIST" ]] && wc -l < "$FAILED_LIST" || echo 0 )
echo "[SUMMARY] assembled n_ok=$n_ok / $TOTAL ; failed=$n_fail"
(( n_ok < 4 )) && { echo "[ERROR] <4 assemblies -- cannot build a meaningful tree"; exit 1; }

# ============================ Stage 3: choose k ==============================
[[ -d "$KSNP4_HOME" ]] || { echo "[ERROR] KSNP4_HOME not found: $KSNP4_HOME"; exit 1; }
export PATH="$KSNP4_HOME:$PATH"
hash -r 2>/dev/null || true
# chmod +x
if [[ ! -x "$KSNP4_HOME/kSNP4" || ! -x "$KSNP4_HOME/Kchooser4" ]]; then
  chmod +x "$KSNP4_HOME"/* 2>/dev/null || true
fi
[[ -x "$KSNP4_HOME/kSNP4" ]] || { echo "[ERROR] kSNP4 not executable: $KSNP4_HOME/kSNP4"; exit 1; }
[[ -x "$KSNP4_HOME/Kchooser4" ]] || { echo "[ERROR] Kchooser4 not executable: $KSNP4_HOME/Kchooser4"; exit 1; }
command -v jellyfish >/dev/null 2>&1 || { echo "[ERROR] jellyfish not found (kSNP4 needs it)"; exit 1; }
echo "[KSNP] using official package: $KSNP4_HOME"
echo "[KSNP] kSNP4=$(command -v kSNP4)"
echo "[KSNP] Kchooser4=$(command -v Kchooser4)"
echo "[KSNP] FastTreeMP=$(command -v FastTreeMP || true)"
echo "[KSNP] jellyfish=$(command -v jellyfish) | tcsh=$(command -v tcsh || true)"
echo "[KSNP] mummer=$(command -v mummer || true)"
if [[ -n "$K_OVERRIDE" ]]; then
  K="$K_OVERRIDE"; echo "[Kchooser4] skipped; using K_OVERRIDE=$K"
else
  ( cd "$OUT_DIR" && "$KSNP4_HOME/Kchooser4" -in "$(basename "$INFILE")" >"$OUT_DIR/logs/kchooser4.log" 2>&1 ) || true
  K="$(grep -oiP 'optimum value of k is\D*\K[0-9]+' "$OUT_DIR"/Kchooser4_*.report 2>/dev/null | head -1)"
  if [[ -z "$K" ]]; then K="$K_DEFAULT"; echo "[Kchooser4] could not parse optimum k -> fallback K=$K"
  else echo "[Kchooser4] optimum k = $K"; fi
fi

# ============================ Stage 4: kSNP4 =================================
rm -rf "$KOUT"; mkdir -p "$KOUT"
echo "[$(date '+%F %T')] running kSNP4 (k=$K, core, ML, CPU=$TOTAL_CPUS)"
kSNP4 -in "$INFILE" -outdir "$KOUT" -k "$K" -core -ML -CPU "$TOTAL_CPUS" \
      >"$OUT_DIR/logs/ksnp4.log" 2>&1 || { echo "[ERROR] kSNP4 failed -- see logs/ksnp4.log"; exit 1; }

echo "----------------------------------------------------------------"
echo "[SUMMARY] kSNP4 done. Trees in $KOUT :"
ls -1 "$KOUT"/*.tre 2>/dev/null || echo "  (no .tre found -- check logs/ksnp4.log)"
echo "[SUMMARY] taxa=$n_ok/$TOTAL ; k=$K ; failed=$n_fail"
(( n_fail > 0 )) && { echo "[SUMMARY] failed:"; cat "$FAILED_LIST"; }
echo "[NEXT] copy $KOUT/tree.core.ML.tre (or tree.core.parsimony.tre) locally -> fig6_ksnp4_tree.R"
rm -f "$LOCK_FILE"
if [[ "$STRICT" == "1" && "$n_fail" -gt 0 ]]; then echo "[STRICT] failures -> exit 1"; exit 1; fi
echo ">>> DONE  [$(date '+%F %T')]"
