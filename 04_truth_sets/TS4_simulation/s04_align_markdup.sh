#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# ============================================================================
# STAGE 4 — Align simulated reads to the FULL reference, matching the REAL pipeline:
#   bwa mem  ->  samtools sort  ->  (read group)  ->  GATK MarkDuplicates  ->  index
# Real versions: BWA 0.7.17, samtools 1.22.1, GATK 4.6.2.0. Here: bwa 0.7.18 / samtools
# 1.19 (alignment+sort are functionally identical across these patches; documented) and
# GATK 4.6.2.0 (exact, from gatk-env). RG added inline via bwa -R (equivalent to the
# real pipeline's separate AddReadGroups step). Uses 80 threads. Resumable + logged.
# ============================================================================
set -uo pipefail
BASE="${SIM_DIR}"
REF="${REF_FASTA}"
RDIR="$BASE/outputs/02_simulated_reads"
ADIR="$BASE/outputs/03_aligned"
PARAMS="$BASE/params.txt"; LOG="$BASE/run.log"
BWA=bwa
SAMTOOLS=samtools
THREADS=80
SORTED="$ADIR/sim_bwa_sorted.bam"
DEDUP="$ADIR/sim.dedup.bam"
mkdir -p "$ADIR"
log(){ echo "[$(date '+%F %T')] [s05] $*" | tee -a "$LOG"; }

if [ -s "$DEDUP" ]; then log "already aligned+dedup -> $DEDUP (delete to re-run)"; exit 0; fi
log "===== STAGE 4 align START (bwa mem -t$THREADS -> sort -> MarkDuplicates) ====="
log "bwa: $($BWA 2>&1 | grep -i version | head -1); samtools: $($SAMTOOLS --version | head -1)"

# --- bwa mem (RG inline) | sort ---------------------------------------------
if [ ! -s "$SORTED" ]; then
  log "aligning ~48.9M pairs with bwa mem + sorting..."
  "$BWA" mem -t "$THREADS" -R '@RG\tID:sim\tSM:SIMSAMPLE\tPL:ILLUMINA\tLB:sim\tPU:sim' \
      "$REF" "$RDIR/sim_R1.fq.gz" "$RDIR/sim_R2.fq.gz" 2>>"$LOG" \
    | "$SAMTOOLS" view -bS - 2>>"$LOG" \
    | "$SAMTOOLS" sort -@ 48 -m 6G -o "$SORTED" - 2>>"$LOG"   # 48x6G=288G in-RAM, avoids disk spill (640G box)
  log "sorted BAM written: $SORTED"
else
  log "sorted BAM exists, skipping bwa/sort"
fi

# --- GATK MarkDuplicates (exact version match, gatk-env) --------------------
log "running GATK MarkDuplicates (gatk-env, v4.6.2.0)..."
source ${CONDA_SH}; set +u; conda activate gatk-env; set -u
gatk --java-options "-Xmx128g" MarkDuplicates \
    -I "$SORTED" -O "$DEDUP" -M "$ADIR/sim.dedup.metrics.txt" \
    --CREATE_INDEX true --TMP_DIR "$ADIR/tmp" 2>&1 | tee -a "$LOG"
set +u; conda deactivate; set -u
[ -s "$DEDUP" ] || { log "ERROR: MarkDuplicates produced no BAM"; exit 3; }

# --- verify -----------------------------------------------------------------
log "verifying (flagstat)..."
"$SAMTOOLS" flagstat "$DEDUP" 2>>"$LOG" | tee -a "$LOG"
MAPPED=$("$SAMTOOLS" flagstat "$DEDUP" 2>/dev/null | awk '/mapped \(/{print $1; exit}')
DUPS=$("$SAMTOOLS" flagstat "$DEDUP" 2>/dev/null | awk '/duplicates/{print $1; exit}')
{ echo "# ---- STAGE 4 alignment ----"
  echo "SIM_BAM=$DEDUP"; echo "SIM_MAPPED=$MAPPED"; echo "SIM_DUPLICATES=$DUPS"
  echo "# STAGE4_DONE $(date '+%F %T')"; } >> "$PARAMS"
rm -rf "$ADIR/tmp"
log "===== STAGE 4 DONE -> $DEDUP (mapped=$MAPPED dups=$DUPS) ====="
