#!/usr/bin/env bash
#SBATCH -J fig6_smallvar_v2
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -N 1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH -o 05_smallvar_v2.out
# =============================================================================
# 05_smallvar_4state.reviewed_v2.slurm  --  STEP 1, second-review hardened.
#
# OUTPUT IS A CANDIDATE VARIANT-EVIDENCE TABLE, NOT high-confidence calls.
# It provides, per strain, 4-state DP/AD/GQ/VAF at jointly discovered sites plus
# a per-gene JOINT callability track so a target-negative is supportable.
#
# Fixes vs reviewed (v1):
#   [blocker] process_strain is NO LONGER called in an `if` condition. Bash
#     suppresses `set -e`/ERR *through* a subshell when the caller is a tested
#     command, so a failed strain could still write `.ok`. Now: `set +e;
#     process_strain; rc=$?; set -e` and require rc==0 AND a fresh `.ok`.
#   [blocker] callability no longer paste()s mosdepth regions+thresholds (the
#     thresholds file has a header -> off-by-one). Now uses `samtools depth`
#     over the 4 BAMs at once and reports the 4-state COMMON callable bp
#     (min(DP_parent,DP01,DP02,DP03) >= threshold) — the same bases readable in
#     ALL states, which is what a "covered but no variant" claim requires.
#   [fix] missing AD (e.g. ".,8") -> VAF=NA, not 1.
#   [fix] soft FILTER + per-sample GQ emitted; output labelled candidate.
#   [fix] mask is now FAIL-CLOSED (set ALLOW_NO_MASK=1 to bypass explicitly).
#   [fix] if GFF3 given and csq fails, the STRAIN fails (no silent unannotated
#         .ok -> consistent annotation state across strains).
#
# Still NOT excluded (report as limits): low-frequency parent subpopulations,
# regulatory-region variants, structural variants; VAF must be read against
# per-chromosome copy number.
#
# Inputs: ALL settings come from the SV_* block in 00_config.sh --
#   SV_PARENT_BAM_DIR / SV_CLONE_BAM_DIR (parents and clones are in SEPARATE
#   run dirs on this cluster), SV_TARGET_BED (fig6_panel25_cds.bed, col4=CNAG),
#   SV_TAG (fresh), optional SV_GFF3 (validated Ensembl-format GFF3 only).
# Submit with no arguments:  sbatch 05_smallvar_4state.reviewed_v2.slurm
# =============================================================================

set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
source ./00_config.sh
apply_strict_mode
activate_conda "$CONDA_ENV_BWASAM"      # bcftools / samtools / bedtools

# All stage-05 settings come from 00_config.sh (SV_* block); override at submit
# time with --export=ALL,SV_TAG=...,SV_MAX_JOBS=... if needed.
# BAMs are in TWO dirs on this cluster: parents vs clones (see 00_config.sh).
PARENT_BAM_DIR="${SV_PARENT_BAM_DIR:?set SV_PARENT_BAM_DIR in 00_config.sh}"
CLONE_BAM_DIR="${SV_CLONE_BAM_DIR:?set SV_CLONE_BAM_DIR in 00_config.sh}"
SMALLVAR_TAG="${SV_TAG:?set SV_TAG (fresh; not the pipeline RUN_TAG)}"
[[ -d "$PARENT_BAM_DIR" ]] || { echo "[ERROR] SV_PARENT_BAM_DIR not a dir: $PARENT_BAM_DIR" >&2; exit 1; }
[[ -d "$CLONE_BAM_DIR"  ]] || { echo "[ERROR] SV_CLONE_BAM_DIR not a dir: $CLONE_BAM_DIR" >&2; exit 1; }

SV_MAX_DEPTH="${SV_MAX_DEPTH:-4000}"
SV_MIN_MAPQ="${SV_MIN_MAPQ:-$VC_MIN_MAPQ}"
SV_MIN_BQ="${SV_MIN_BQ:-$VC_MIN_BQ}"
SV_DP_CALLABLE="${SV_DP_CALLABLE:-10}"
SV_QUAL_SOFT="${SV_QUAL_SOFT:-20}"         # QUAL below this -> FILTER=LowQual (soft)
TARGET_BED="${SV_TARGET_BED:-}"
GFF3="${SV_GFF3:-}"
INCLUDE_H99="${SV_INCLUDE_H99:-0}"
ALLOW_NO_MASK="${SV_ALLOW_NO_MASK:-0}"
FORCE="${FORCE:-1}"

# ---- parallelism (FIFO semaphore, same pattern as 01_map / 02_mosdepth) ------
# Each strain is ~2 concurrent procs (mpileup | call), then 1 (samtools depth),
# and is I/O-bound rather than CPU-bound, so 16 strains on 24 CPUs is a sensible
# oversubscription. Memory is the thing to watch: mpileup -d $SV_MAX_DEPTH over
# 4 BAMs x MAX_JOBS strains shares --mem (64G/16 ~ 4G per strain here). If jobs
# get OOM-killed, lower SV_MAX_JOBS rather than raising --mem.
THREADS="${SLURM_CPUS_PER_TASK:-24}"
MAX_JOBS="${SV_MAX_JOBS:-20}"

# ---- mask: FAIL-CLOSED ------------------------------------------------------
if [[ -z "${MASK_BED:-}" || ! -s "$MASK_BED" ]]; then
  if [[ "$ALLOW_NO_MASK" == "1" ]]; then echo "[WARN] no mask; ALLOW_NO_MASK=1 -> genome-wide."; MASK_BED=""
  else echo "[ERROR] mappability mask missing (${MASK_BED:-unset}). Set ALLOW_NO_MASK=1 to bypass." >&2; exit 1; fi
fi
# ---- target BED (required for the callability the negative claim needs) ------
[[ -n "$TARGET_BED" && -s "$TARGET_BED" ]] || { echo "[ERROR] TARGET_BED (fig6_panel25_cds.bed) required." >&2; exit 1; }

OUT="$WORK_ROOT/fig6_smallvar_${SMALLVAR_TAG}"
[[ -e "$OUT" && "$FORCE" != "1" ]] && { echo "[ERROR] $OUT exists; new SMALLVAR_TAG or FORCE=1." >&2; exit 1; }
VCF_DIR="$OUT/vcf"; TAB_DIR="$OUT/tables"; LOG_DIR="$OUT/logs"; PS_DIR="$OUT/perstrain"
mkdir -p "$VCF_DIR" "$TAB_DIR" "$LOG_DIR" "$PS_DIR"
print_config_banner
print_smallvar_banner

TARGETS_OPT=(); [[ -n "$MASK_BED" ]] && TARGETS_OPT=(-T "$MASK_BED")
TGT_MAP_BED="$OUT/target_mappable.bed"
# NOTE: build the intersection EXPLICITLY as BED4. The downstream callability
# awk keys on column 4 (CNAG); if col4 were lost the job would still "succeed"
# but silently report mappable_bp=0 / callability 0-NA for every gene.
if [[ -n "$MASK_BED" ]]; then
  bedtools intersect -a "$TARGET_BED" -b "$MASK_BED" -wa -wb \
    | awk 'BEGIN{OFS="\t"}{
        s=($2>$6?$2:$6);
        e=($3<$7?$3:$7);
        if(s<e) print $1,s,e,$4
      }' > "$TGT_MAP_BED"
else
  cp "$TARGET_BED" "$TGT_MAP_BED"
fi

[[ -s "$TGT_MAP_BED" ]] || { echo "[ERROR] empty target x mask intersection" >&2; exit 1; }
awk 'NF!=4 || $4==""{exit 1}' "$TGT_MAP_BED" || { echo "[ERROR] target_mappable.bed is not valid BED4" >&2; exit 1; }
echo "[$(date '+%F %T')] target x mask -> $TGT_MAP_BED ($(wc -l < "$TGT_MAP_BED") BED4 intervals, $(awk '{s+=$3-$2} END{print s+0}' "$TGT_MAP_BED") bp)"

DO_CSQ=0
if [[ -n "$GFF3" ]]; then
  [[ -s "$GFF3" ]] && DO_CSQ=1 || { echo "[ERROR] GFF3 not found: $GFF3" >&2; exit 1; }
fi

STRAINS=(101 106 126 128 147 150 160 162 169 174 196 197 198 208 213 215 217 225 226)
(( INCLUDE_H99 )) && STRAINS+=(H99)

sample_id() { [[ "$2" == "parent" ]] && echo "Cn$1" || echo "Cn$1$2"; }
# parent BAMs and clone (01/02/03) BAMs live in different run dirs.
bam_path()  {
  if [[ "$2" == "parent" ]]; then echo "$PARENT_BAM_DIR/$(sample_id "$1" "$2").bam"
  else                            echo "$CLONE_BAM_DIR/$(sample_id "$1" "$2").bam"; fi
}
bai_ok()    { [[ -s "$1.bai" || -s "${1%.bam}.bai" ]]; }
rg_sm()     { samtools view -H "$1" 2>/dev/null | awk -F'\t' '/^@RG/{for(i=1;i<=NF;i++) if($i ~ /^SM:/){sub(/^SM:/,"",$i);print $i;exit}}'; }

# ---- STRICT 4-state preflight ----------------------------------------------
MANIFEST="$OUT/input_manifest.tsv"
echo -e "strain\tstate\tsample\tbam\texists\tquickcheck\tindexed\trg_sm_ok" > "$MANIFEST"
declare -a RUN_STRAINS=() INCOMPLETE=()
for s in "${STRAINS[@]}"; do
  ok_all=1
  for st in parent 01 02 03; do
    b="$(bam_path "$s" "$st")"; sid="$(sample_id "$s" "$st")"; ex=0 qc=0 ix=0 rg=0
    [[ -s "$b" ]] && ex=1
    if (( ex )); then samtools quickcheck "$b" 2>/dev/null && qc=1; bai_ok "$b" && ix=1; [[ "$(rg_sm "$b")" == "$sid" ]] && rg=1; fi
    echo -e "${s}\t${st}\t${sid}\t${b}\t${ex}\t${qc}\t${ix}\t${rg}" >> "$MANIFEST"
    (( ex && qc && ix && rg )) || ok_all=0
  done
  if (( ok_all )); then RUN_STRAINS+=("$s"); else INCOMPLETE+=("$s"); echo "[WARN] $s incomplete -> excluded."; fi
done
echo "[$(date '+%F %T')] runnable (${#RUN_STRAINS[@]}): ${RUN_STRAINS[*]:-none}"
(( ${#RUN_STRAINS[@]} )) || { echo "[ERROR] no strain passed strict 4-state preflight." >&2; exit 1; }

# ---- per-strain (subshell set -e HOLDS because we never call this in an if) --
process_strain() {
  local s="$1" sdir="$PS_DIR/Cn${s}"
  mkdir -p "$sdir"; rm -f "$sdir"/*.part "$sdir/.ok"
  (
    set -e
    trap 'echo "[ERR] Cn'"$s"' failed near line $LINENO" >&2' ERR
    local bp bo b1 b2 b3
    bp="$(bam_path "$s" parent)"; b1="$(bam_path "$s" 01)"; b2="$(bam_path "$s" 02)"; b3="$(bam_path "$s" 03)"
    local raw="$sdir/raw.vcf.gz" norm="$sdir/norm.vcf.gz" flt="$sdir/flt.vcf.gz" final CSQ_TOK='.'

    bcftools mpileup -f "$REF_FA" ${TARGETS_OPT[@]+"${TARGETS_OPT[@]}"} \
        -q "$SV_MIN_MAPQ" -Q "$SV_MIN_BQ" -d "$SV_MAX_DEPTH" \
        -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP -Ou "$bp" "$b1" "$b2" "$b3" \
      | bcftools call -m -v -G - -f GQ --ploidy "${VC_PLOIDY:-2}" -Oz -o "$raw"
    bcftools index -t "$raw"
    bcftools norm -f "$REF_FA" -m -any --keep-sum AD -Oz -o "$norm" "$raw"
    bcftools filter -m + -s LowQual -e "QUAL<${SV_QUAL_SOFT}" -Oz -o "$flt" "$norm"
    bcftools index -t "$flt"
    final="$flt"

    if (( DO_CSQ )); then
      local ann="$sdir/csq.vcf.gz"
      bcftools csq --local-csq -f "$REF_FA" -g "$GFF3" -Oz -o "$ann" "$flt"   # failure -> set -e aborts strain
      bcftools index -t "$ann"; final="$ann"; CSQ_TOK='%INFO/BCSQ'
    fi

    local fmt='[%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%QUAL\t%FILTER\t'"$CSQ_TOK"'\t%SAMPLE\t%GT\t%GQ\t%DP\t%AD\t%ADF\t%ADR\t%SP\n]'
    bcftools query -f "$fmt" "$final" \
      | awk -v S="$s" 'BEGIN{OFS="\t"}{
          samp=$9; st="NA";
          if(samp=="Cn"S)st="parent"; else if(samp=="Cn"S"01")st="01";
          else if(samp=="Cn"S"02")st="02"; else if(samp=="Cn"S"03")st="03";
          n=split($13,ad,","); r=ad[1]; a=ad[2];
          if(r!~/^[0-9]+$/||a!~/^[0-9]+$/){ rr="NA"; aa="NA"; vaf="NA" }
          else{ rr=r; aa=a; vaf=(a+r>0)?sprintf("%.4f",a/(a+r)):"NA" }
          split($14,f,","); fa=(f[2]~/^[0-9]+$/)?f[2]:"NA";
          split($15,v,","); ra=(v[2]~/^[0-9]+$/)?v[2]:"NA";
          print S,$1,$2,$3,$4,$5,$6,$7,$8,samp,st,$10,$11,$12,rr,aa,fa,ra,$16,vaf
        }' > "$sdir/long.part"

    # ---- 4-state JOINT callability over target x mask (samtools depth) -------
    # depth cols: chrom pos DPparent DP01 DP02 DP03 (bam order below)
    samtools depth -a -q "$SV_MIN_BQ" -Q "$SV_MIN_MAPQ" -b "$TGT_MAP_BED" "$bp" "$b1" "$b2" "$b3" \
      | awk -v S="$s" -v THR="$SV_DP_CALLABLE" -v CDSF="$TARGET_BED" -v MAPF="$TGT_MAP_BED" '
          BEGIN{ OFS="\t";   # REQUIRED: without it print uses space and the
                             # callability table becomes space-delimited while
                             # its header is tab-delimited -> unparsable in R.
            while((getline l < CDSF)>0){ split(l,c,"\t"); CDS[c[4]]+=c[3]-c[2] }
            while((getline l < MAPF)>0){ split(l,m,"\t"); MAP[m[4]]+=m[3]-m[2];
              for(p=m[2]+1;p<=m[3];p++) G[m[1] SUBSEP p]=m[4] }
          }
          { g=G[$1 SUBSEP $2]; if(g=="") next;
            if($3>=THR)cp[g]++; if($4>=THR)c1[g]++; if($5>=THR)c2[g]++; if($6>=THR)c3[g]++;
            mn=$3; if($4<mn)mn=$4; if($5<mn)mn=$5; if($6<mn)mn=$6; if(mn>=THR)cc[g]++ }
          END{ for(g in CDS){ cds=CDS[g]; mp=MAP[g]+0;
            print S,g,cds,mp,(cds>0)?sprintf("%.3f",mp/cds):"NA",
                  cp[g]+0,c1[g]+0,c2[g]+0,c3[g]+0,cc[g]+0,
                  (mp>0)?sprintf("%.3f",(cc[g]+0)/mp):"NA",
                  (cds>0)?sprintf("%.3f",(cc[g]+0)/cds):"NA" } }' \
      > "$sdir/callable.part"

    mv "$sdir/long.part"     "$sdir/Cn${s}.long.tsv"
    mv "$sdir/callable.part" "$sdir/Cn${s}.callable.tsv"
    cp "$final" "$VCF_DIR/Cn${s}.vcf.gz"; cp "$final.tbi" "$VCF_DIR/Cn${s}.vcf.gz.tbi"
    : > "$sdir/.ok"
  ) >>"$LOG_DIR/Cn${s}.log" 2>&1
}

# ---- parallel execution: FIFO semaphore --------------------------------------
# FAILURE SEMANTICS IN PARALLEL: a background job's exit code cannot be captured
# by the parent, so each worker writes its own .rc, and process_strain still only
# publishes .ok after every step succeeded. A strain counts as OK iff rc==0 AND
# .ok exists (checked after wait). process_strain is STILL never called from an
# `if`/`&&` condition -- that would suppress `set -e` inside its subshell and let
# a failed strain publish .ok (the v1 blocker).
SEM_FIFO="$(mktemp -u)"; mkfifo "$SEM_FIFO"; exec 3<>"$SEM_FIFO"; rm -f "$SEM_FIFO"
for ((i=0;i<MAX_JOBS;i++)); do echo >&3; done

run_one() {
  local s="$1" sdir="$PS_DIR/Cn${s}"
  mkdir -p "$sdir"; rm -f "$sdir/.rc"
  set +e
  process_strain "$s"            # plain call: keeps the inner subshell's set -e live
  local rc=$?
  set -e
  printf '%s\n' "$rc" > "$sdir/.rc"
  if (( rc == 0 )); then echo "[$(date '+%F %T')]     Cn${s} done (rc=0)"
  else                   echo "[$(date '+%F %T')]     Cn${s} FAILED (rc=$rc)"; fi
}

echo "[$(date '+%F %T')] launching ${#RUN_STRAINS[@]} strains, ${MAX_JOBS} in parallel on ${THREADS} CPUs"
for s in "${RUN_STRAINS[@]}"; do
  read -u 3
  (
    trap 'echo >&3' EXIT
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    run_one "$s"
  ) &
done
wait
exec 3>&-
echo "[$(date '+%F %T')] all workers finished; verifying rc + .ok"

# ---- post-hoc verification (rc AND .ok) -------------------------------------
declare -a FAILED=()
for s in "${RUN_STRAINS[@]}"; do
  rc="$(cat "$PS_DIR/Cn${s}/.rc" 2>/dev/null || echo 99)"
  if [[ "$rc" == "0" && -f "$PS_DIR/Cn${s}/.ok" ]]; then
    echo "[$(date '+%F %T')] Cn${s} OK"
  else
    echo "[ERROR] Cn${s} FAILED rc=${rc} (see $LOG_DIR/Cn${s}.log)"
    FAILED+=("$s"); rm -f "$PS_DIR/Cn${s}"/*.part "$PS_DIR/Cn${s}/.ok"
  fi
done

# ---- merge ONLY this run's .ok ----------------------------------------------
LONG_HEADER="strain	chrom	pos	ref	alt	type	qual	filter	bcsq	sample	state	gt	gq	dp	ad_ref	ad_alt	adf_alt	adr_alt	sp	vaf"
CALL_HEADER="strain	cnag	cds_bp	mappable_bp	mappable_frac	call_parent	call_01	call_02	call_03	common_callable_bp	common_frac_mappable	common_frac_cds"
: > "$TAB_DIR/long.body"; : > "$TAB_DIR/call.body"
for s in "${RUN_STRAINS[@]}"; do
  [[ -f "$PS_DIR/Cn${s}/.ok" ]] || continue
  cat "$PS_DIR/Cn${s}/Cn${s}.long.tsv"     >> "$TAB_DIR/long.body"
  cat "$PS_DIR/Cn${s}/Cn${s}.callable.tsv" >> "$TAB_DIR/call.body"
done
{ echo -e "$LONG_HEADER"; cat "$TAB_DIR/long.body"; } > "$TAB_DIR/smallvar_4state_long.tsv"
{ echo -e "$CALL_HEADER"; cat "$TAB_DIR/call.body"; } > "$TAB_DIR/target_panel_callability.tsv"
rm -f "$TAB_DIR/long.body" "$TAB_DIR/call.body"

FAIL_TSV="$OUT/FAILED_SAMPLES.tsv"; echo -e "strain\treason" > "$FAIL_TSV"
# NB bare "${arr[@]}" on an EMPTY array is an "unbound variable" error under
# set -u in bash < 4.4 (this cluster). The +"" idiom expands to nothing instead.
# This only ever fired when EVERY strain succeeded, which is why it survived testing.
for s in ${INCOMPLETE[@]+"${INCOMPLETE[@]}"}; do echo -e "${s}\tincomplete_preflight" >> "$FAIL_TSV"; done
for s in ${FAILED[@]+"${FAILED[@]}"};         do echo -e "${s}\tprocessing_failed"   >> "$FAIL_TSV"; done

echo "[$(date '+%F %T')] long: $(( $(wc -l < "$TAB_DIR/smallvar_4state_long.tsv") - 1 )) rows"
echo "[$(date '+%F %T')] callability: $(( $(wc -l < "$TAB_DIR/target_panel_callability.tsv") - 1 )) gene-rows"
echo "[$(date '+%F %T')] manifest=$MANIFEST failures=$FAIL_TSV"
if (( ${#INCOMPLETE[@]} || ${#FAILED[@]} )); then
  echo "[$(date '+%F %T')] NON-ZERO EXIT: incomplete=${#INCOMPLETE[@]} failed=${#FAILED[@]}" >&2; exit 1
fi
echo "[$(date '+%F %T')] all target strains OK."
