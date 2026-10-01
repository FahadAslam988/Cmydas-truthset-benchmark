#!/usr/bin/env bash
# Adapter/quality trimming with fastp 0.23.4 (command as run; default quality settings).
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
set -euo pipefail; mkdir -p "$TRIM_DIR"
for pair in "CH-NORMS1:W1-A" "Ab-NormS2:S1"; do s=${pair%%:*}; f=${pair##*:}
  fastp -i "$FASTQ_DIR/${f}_1.fastq.gz" -I "$FASTQ_DIR/${f}_2.fastq.gz" \
        -o "$TRIM_DIR/${f}_1_trimmed.fastq.gz" -O "$TRIM_DIR/${f}_2_trimmed.fastq.gz" \
        --detect_adapter_for_pe --thread "$THREADS" --html "$TRIM_DIR/${s}_fastp.html" --json "$TRIM_DIR/${s}_fastp.json"
done
