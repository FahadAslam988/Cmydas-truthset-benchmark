#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# Per-contig ART worker (called in parallel by s03b). Simulates PE reads from one
# contig of a haplotype genome using the CUSTOM error profile.
# args: GENOME CONTIG TAG OUTDIR PROFR1 PROFR2 DEPTH INSMEAN INSSD SEED
set -uo pipefail
GENOME="$1"; CONTIG="$2"; TAG="$3"; OUTDIR="$4"; PROFR1="$5"; PROFR2="$6"
DEPTH="$7"; INSMEAN="$8"; INSSD="$9"; SEED="${10}"
SAMTOOLS=samtools
ART=art_illumina

CHUNK="$OUTDIR/$TAG.fa"
# extract the contig, rename header to the unique TAG (keeps pooled read names unique)
"$SAMTOOLS" faidx "$GENOME" "$CONTIG" | sed "1s/.*/>$TAG/" > "$CHUNK"
# simulate: custom R1/R2 profiles, 2x151 PE, given fold-coverage + insert size
"$ART" -1 "$PROFR1" -2 "$PROFR2" -na -p -l 151 -f "$DEPTH" \
       -m "$INSMEAN" -s "$INSSD" -rs "$SEED" -i "$CHUNK" -o "$OUTDIR/${TAG}_R" \
       > "$OUTDIR/${TAG}.artlog" 2>&1
rc=$?
rm -f "$CHUNK"
if [ $rc -ne 0 ] || [ ! -s "$OUTDIR/${TAG}_R1.fq" ]; then
  echo "WORKER_FAIL $TAG rc=$rc" >> "$OUTDIR/.failures"
fi
exit 0
