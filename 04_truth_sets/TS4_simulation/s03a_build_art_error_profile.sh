#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# ============================================================================
# STAGE 3a — Build a CUSTOM ART error profile from the REAL C. mydas NovaSeq reads,
# so simulated sequencing errors are empirically derived from the actual data
# (gold-standard realism). Steps:
#   1. extract a sample of real read pairs from CH-NORMS1.dedup.bam (samtools fastq
#      RESTORES original read orientation -> correct per-cycle quality);
#   2. keep only full-length 151 bp reads (99% of data; profiler needs uniform length);
#   3. art_profiler_illumina -> R1/R2 quality profiles (uses 80 threads).
# Fallback: if the profiler fails, Stage 3b uses built-in HSXn (documented limitation).
# Resumable + logged.
# ============================================================================
set -uo pipefail

BASE="${SIM_DIR}"
BAM="${BAM_DIR}/CH-NORMS1.dedup.bam"
SAMTOOLS=samtools
PROFDIR="$BASE/outputs/02_simulated_reads/error_profile"
FQDIR="$PROFDIR/real_fastq"
PROFNAME="$PROFDIR/cmydas_novaseq"      # -> cmydas_novaseq_R1.txt / _R2.txt
LOG="$BASE/run.log"
PARAMS="$BASE/params.txt"
REGION="NC_057855.1:1-15000000"          # ~15 Mb sample region
MAXREADS=4000000                          # cap per mate (plenty for a stable profile)
THREADS=80
log(){ echo "[$(date '+%F %T')] [s03a] $*" | tee -a "$LOG"; }

if [ -s "${PROFNAME}R1.txt" ] && [ -s "${PROFNAME}R2.txt" ]; then
  log "custom profile already built (${PROFNAME}R1/R2.txt) — skipping"; exit 0
fi
mkdir -p "$FQDIR"
log "===== STAGE 3a build custom error profile START ====="

# --- 1+2. extract real read pairs (oriented) and keep only 151 bp ----------
log "extracting real reads from $REGION (samtools collate|fastq), filtering to 151 bp..."
"$SAMTOOLS" view -b -f 2 -F 0x900 "$BAM" "$REGION" 2>>"$LOG" \
  | "$SAMTOOLS" collate -u -O -@ 16 - 2>>"$LOG" \
  | "$SAMTOOLS" fastq -n -1 "$FQDIR/tmp_1.fq" -2 "$FQDIR/tmp_2.fq" -0 /dev/null -s /dev/null -@ 16 - 2>>"$LOG"

keep151(){ # in.fq out.fq  (keep 4-line records whose SEQ length==151, cap MAXREADS)
  awk -v max="$MAXREADS" 'NR%4==1{h=$0} NR%4==2{s=$0} NR%4==3{p=$0}
       NR%4==0{q=$0; if(length(s)==151){c++; print h"\n"s"\n"p"\n"q; if(c>=max) exit}}' "$1" > "$2"
}
log "filtering R1/R2 to 151 bp (cap ${MAXREADS})..."
keep151 "$FQDIR/tmp_1.fq" "$FQDIR/real_1.fq"
keep151 "$FQDIR/tmp_2.fq" "$FQDIR/real_2.fq"
rm -f "$FQDIR/tmp_1.fq" "$FQDIR/tmp_2.fq"
N1=$(( $(wc -l < "$FQDIR/real_1.fq") / 4 )); N2=$(( $(wc -l < "$FQDIR/real_2.fq") / 4 ))
log "kept R1=$N1  R2=$N2 reads (151 bp)"
if [ "$N1" -lt 100000 ] || [ "$N2" -lt 100000 ]; then
  log "ERROR: too few 151bp reads for a profile (R1=$N1 R2=$N2). Use HSXn fallback in 3b."; exit 3
fi

# --- 3. build the profile --------------------------------------------------
source ${CONDA_SH}; set +u; conda activate art-env; set -u
log "running art_profiler_illumina ($THREADS threads)..."
art_profiler_illumina "$PROFNAME" "$FQDIR" fq "$THREADS" 2>&1 | tee -a "$LOG"

if [ -s "${PROFNAME}R1.txt" ] && [ -s "${PROFNAME}R2.txt" ]; then
  log "custom profile built: ${PROFNAME}R1.txt / ${PROFNAME}R2.txt"
  { echo "# ---- STAGE 3a custom ART error profile ----"
    echo "ART_PROFILE_R1=${PROFNAME}R1.txt"; echo "ART_PROFILE_R2=${PROFNAME}R2.txt"
    echo "ART_PROFILE_READS_R1=$N1"; echo "ART_PROFILE_READS_R2=$N2"
    echo "# STAGE3a_DONE $(date '+%F %T')"; } >> "$PARAMS"
  # free the raw fastq sample (keep only the profiles)
  rm -f "$FQDIR/real_1.fq" "$FQDIR/real_2.fq"; rmdir "$FQDIR" 2>/dev/null || true
  log "===== STAGE 3a DONE ====="
else
  log "ERROR: art_profiler_illumina did not produce R1/R2 profiles -> Stage 3b will fall back to HSXn."; exit 4
fi
