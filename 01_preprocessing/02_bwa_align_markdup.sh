#!/usr/bin/env bash
# Alignment and duplicate marking, as recorded in the @PG lines of the deduplicated BAM headers:
#   bwa mem 0.7.17-r1188 (-t 15, default settings) | samtools 1.10 view -bS | samtools sort
#   read group (ID/SM = sample, LB lib1, PL illumina, PU unit1) -> GATK 4.6.2.0 MarkDuplicates
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
set -euo pipefail; mkdir -p "$BAM_DIR"
for pair in "CH-NORMS1:W1-A" "Ab-NormS2:S1"; do s=${pair%%:*}; f=${pair##*:}
  bwa mem -t "$THREADS" "$REF_FASTA" "$TRIM_DIR/${f}_1_trimmed.fastq.gz" "$TRIM_DIR/${f}_2_trimmed.fastq.gz" \
    | samtools view -bS - | samtools sort -@ "$THREADS" -o "$BAM_DIR/${s}_sorted.bam"
  gatk AddOrReplaceReadGroups -I "$BAM_DIR/${s}_sorted.bam" -O "$BAM_DIR/${s}_sorted_rg.bam" \
    --RGID "$s" --RGSM "$s" --RGLB lib1 --RGPL illumina --RGPU unit1
  gatk MarkDuplicates -I "$BAM_DIR/${s}_sorted_rg.bam" -O "$BAM_DIR/${s}.dedup.bam" \
    -M "$BAM_DIR/${s}.metrics.txt" --CREATE_INDEX true
done
