#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# ============================================================================
# STAGE 0 — Calibrate simulation to the REAL data (both BAMs: normal + abnormal)
# Measures mean depth, insert size (mean/sd) and read length so every simulation
# parameter traces back to the real NovaSeq data (not invented). Light compute.
# Reproducible + resumable: checkpoints, logs, fixed sampling regions.
# ============================================================================
set -euo pipefail

# --- paths (edit here only) -------------------------------------------------
BASE="${SIM_DIR}"
REF="${REF_FASTA}"
BAM_NORMAL="${BAM_DIR}/CH-NORMS1.dedup.bam"
BAM_ABNORMAL="${BAM_DIR}/Ab-NormS2.dedup.bam"
SAMTOOLS="samtools"

PARAMS="$BASE/params.txt"
LOG="$BASE/run.log"
N_CONTIGS=5          # sample the N largest contigs
WIN=10000000         # 10 Mb window per contig (fast, representative)

# --- checkpoint: skip if already done --------------------------------------
if grep -q "^# STAGE0_DONE" "$PARAMS" 2>/dev/null; then
  echo "[s00] already complete -> $PARAMS (delete the file to re-run)"; exit 0
fi

log(){ echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }
log "===== STAGE 0 calibrate START ====="
log "samtools: $("$SAMTOOLS" --version | head -1)"

# --- pick the N largest contigs from the .fai; build a sampling BED ---------
mapfile -t TOPC < <(sort -k2,2nr "$REF.fai" | head -n "$N_CONTIGS" | cut -f1)
log "sampling contigs: ${TOPC[*]} (first ${WIN} bp each)"

SBED="$BASE/outputs/calibration_regions.bed"
FIRST_REGION=""
: > "$SBED"
for c in "${TOPC[@]}"; do
  len=$(awk -v c="$c" '$1==c{print $2}' "$REF.fai")
  end=$(( len < WIN ? len : WIN ))
  printf "%s\t0\t%s\n" "$c" "$end" >> "$SBED"     # BED is 0-based half-open
  [ -z "$FIRST_REGION" ] && FIRST_REGION="$c:1-$end"
done

# --- measure one BAM -------------------------------------------------------
measure(){
  local label="$1" bam="$2"
  log "--- $label: $bam ---"
  # mean depth over the sampled BED regions via bedcov (sum depth / sum bp)
  local depth
  depth=$("$SAMTOOLS" bedcov "$SBED" "$bam" 2>>"$LOG" \
          | awk '{d+=$4; b+=($3-$2)} END{if(b>0) printf "%.1f", d/b; else print "NA"}')
  # insert size + read length from samtools stats over first sampled region
  local stats ins insd rlen
  stats=$("$SAMTOOLS" stats "$bam" "$FIRST_REGION" 2>>"$LOG")
  ins=$(awk -F'\t' '/^SN\tinsert size average:/{print $3}'            <<<"$stats")
  insd=$(awk -F'\t' '/^SN\tinsert size standard deviation:/{print $3}' <<<"$stats")
  rlen=$(awk -F'\t' '/^SN\taverage length:/{print $3}'                 <<<"$stats")
  log "$label -> depth=${depth}x  insert=${ins}±${insd}  readlen=${rlen}"
  printf "%s\tdepth=%s\tinsert_mean=%s\tinsert_sd=%s\treadlen=%s\n" \
         "$label" "$depth" "$ins" "$insd" "$rlen" >> "$PARAMS.tmp"
}

# --- run (samtools coverage needs -r before EACH region; expand array) ------
# NOTE: "${REGIONS[@]/#/-r }" turns each "reg" into "-r reg"
: > "$PARAMS.tmp"
{
  echo "# STAGE 0 calibration — measured from real BAMs $(date '+%F %T')"
  echo "# sampled contigs: ${TOPC[*]}  window: ${WIN} bp"
} >> "$PARAMS.tmp"

measure "NORMAL"   "$BAM_NORMAL"
measure "ABNORMAL" "$BAM_ABNORMAL"

# --- decide: one simulation or two? ----------------------------------------
dN=$(awk -F'depth=' '/NORMAL/{split($2,a,"\t"); print a[1]}'   "$PARAMS.tmp")
dA=$(awk -F'depth=' '/ABNORMAL/{split($2,a,"\t"); print a[1]}' "$PARAMS.tmp")
decision=$(awk -v n="$dN" -v a="$dA" 'BEGIN{
  d=(n>a)?n-a:a-n; rel=(n>0)?d/n:1;
  if(rel<=0.20) print "SIMILAR -> ONE simulation at mean(depth)";
  else          print "DIFFERENT -> TWO simulations, one per sample depth" }')
echo "# DEPTH_DECISION: normal=${dN}x abnormal=${dA}x -> $decision" >> "$PARAMS.tmp"
log "DECISION: normal=${dN}x abnormal=${dA}x -> $decision"

echo "# STAGE0_DONE $(date '+%F %T')" >> "$PARAMS.tmp"
mv "$PARAMS.tmp" "$PARAMS"
log "===== STAGE 0 DONE -> $PARAMS ====="
echo
echo ">>> Review $PARAMS, then we proceed to s01_subset.sh"
