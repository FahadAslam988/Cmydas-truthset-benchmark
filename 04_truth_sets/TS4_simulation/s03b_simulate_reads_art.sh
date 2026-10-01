#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# ============================================================================
# STAGE 3b — Simulate PE reads with ART (custom C. mydas error profile).
# hap1 + hap2 each at HALF the target depth (35x) -> pooled ~70x diploid.
# Parallelized per-contig across 80 threads. Reproducible (per-chunk seeds).
# Output: outputs/02_simulated_reads/sim_R1.fq.gz , sim_R2.fq.gz
# Resumable + logged.
# ============================================================================
set -uo pipefail

BASE="${SIM_DIR}"
PARAMS="$BASE/params.txt"
GDIR="$BASE/outputs/01_mutated_genome"
RDIR="$BASE/outputs/02_simulated_reads"
WORK="$RDIR/chunks"
LOG="$BASE/run.log"
BGZIP=bgzip
SEEDBASE=20260714
THREADS=80
log(){ echo "[$(date '+%F %T')] [s03b] $*" | tee -a "$LOG"; }

getp(){ grep -E "^$1=" "$PARAMS" | head -1 | cut -d= -f2- | awk '{print $1}'; }
PROFR1=$(getp ART_PROFILE_R1); PROFR2=$(getp ART_PROFILE_R2)
INSMEAN=$(getp SIM_INSERT_MEAN); INSSD=$(getp SIM_INSERT_SD)
COV=$(getp SIM_COVERAGE); DEPTH_HALF=$(( COV / 2 ))
R1="$RDIR/sim_R1.fq.gz"; R2="$RDIR/sim_R2.fq.gz"

if [ -s "$R1" ] && [ -s "$R2" ]; then log "reads already simulated -> $R1 (delete to re-run)"; exit 0; fi
[ -s "$PROFR1" ] && [ -s "$PROFR2" ] || { log "ERROR: custom profiles missing (run s03a)"; exit 2; }
mkdir -p "$WORK"; : > "$WORK/.failures"
log "===== STAGE 3b ART reads START (per-hap ${DEPTH_HALF}x, insert ${INSMEAN}±${INSSD}, ${THREADS} threads) ====="

# --- build job list: one line per (hap, contig); TAG keeps names unique -----
JOBS="$WORK/jobs.tsv"; TAGS="$WORK/tags.order"; : > "$JOBS"; : > "$TAGS"
mkjobs(){ # hap-label  genome.fa
  local lab="$1" gen="$2" i=0
  while read -r c _; do
    i=$((i+1)); local tag="${lab}_$(printf '%05d' $i)"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$gen" "$c" "$tag" "$WORK" "$PROFR1" "$PROFR2" "$DEPTH_HALF" "$INSMEAN" "$INSSD" "$((SEEDBASE+i))" >> "$JOBS"
    echo "$tag" >> "$TAGS"
  done < "$gen.fai"
}
mkjobs h1 "$GDIR/hap1.simseq.genome.fa"
mkjobs h2 "$GDIR/hap2.simseq.genome.fa"
NJOBS=$(wc -l < "$JOBS")
log "launching $NJOBS ART jobs (72 contigs x 2 haplotypes) across $THREADS threads..."

# --- run all chunks in parallel --------------------------------------------
# each line has 10 tab-separated fields; convert BOTH tabs and newlines to NUL so
# xargs -0 -n10 groups exactly 10 args per worker call (paths here contain no spaces)
tr '\t\n' '\0\0' < "$JOBS" | xargs -0 -n 10 -P "$THREADS" bash "$BASE/scripts/art_chunk_worker.sh"

if [ -s "$WORK/.failures" ]; then
  log "ERROR: some ART chunks failed:"; cat "$WORK/.failures" | tee -a "$LOG"; exit 3
fi
log "all ART chunks done. Concatenating (order-preserving) + compressing..."

# --- concatenate in tag order (R1 and R2 identical order -> pairing preserved)
: > "$RDIR/sim_R1.fq"; : > "$RDIR/sim_R2.fq"
while read -r t; do
  cat "$WORK/${t}_R1.fq" >> "$RDIR/sim_R1.fq"
  cat "$WORK/${t}_R2.fq" >> "$RDIR/sim_R2.fq"
done < "$TAGS"
"$BGZIP" -@ 16 -f "$RDIR/sim_R1.fq"
"$BGZIP" -@ 16 -f "$RDIR/sim_R2.fq"

# --- verify + record --------------------------------------------------------
NR1=$(( $(zcat "$R1" | wc -l) / 4 )); NR2=$(( $(zcat "$R2" | wc -l) / 4 ))
GBP=$(awk '{s+=$2}END{print s}' "$GDIR/subset.fa.fai")
EFFCOV=$(awk -v n=$NR1 -v g=$GBP 'BEGIN{printf "%.1f", (n*2*151)/g}')   # R1+R2 x151 / genome
log "reads: R1=$NR1 R2=$NR2 pairs; effective pooled coverage ~= ${EFFCOV}x (target ${COV}x)"
{ echo "# ---- STAGE 3b ART reads ----"
  echo "SIM_R1=$R1"; echo "SIM_R2=$R2"; echo "SIM_READ_PAIRS=$NR1"; echo "SIM_EFF_COVERAGE=$EFFCOV"
  echo "# STAGE3b_DONE $(date '+%F %T')"; } >> "$PARAMS"
rm -rf "$WORK"
log "===== STAGE 3b DONE -> $R1 , $R2 ====="
