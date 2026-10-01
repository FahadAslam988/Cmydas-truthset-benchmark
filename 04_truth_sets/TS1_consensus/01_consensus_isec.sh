#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# Stage 4 (research standard): CONSENSUS via allele-aware bcftools isec on NORMALIZED SNP files.
# Levels: >=3 (majority, primary), >=4, =5. (Exact =3/=4 dropped: not valid truth sets.)
# Refs: Zook 2014/2019 (GIAB consensus); Danecek 2021 (bcftools isec).
set -uo pipefail
THREADS=16

NSNP="${WORK_DIR}/04_norm_snps"
OUT="${WORK_DIR}/05_consensus"
mkdir -p "$OUT"
LOG="$OUT/consensus_log.txt"; : > "$LOG"
SUM="$OUT/consensus_counts.tsv"
printf "Group\tLevel\tConsensus_SNP_sites\n" > "$SUM"

do_group() {
  local group="$1"; shift
  local files=( "$@" )
  # level label pairs
  for pair in "+3:3plus" "+4:4plus" "=5:5tools"; do
    local n="${pair%%:*}"; local label="${pair##*:}"
    local d="$OUT/isec_${group}_${label}"
    rm -rf "$d"
    echo "[$(date +%H:%M:%S)] isec $group $n -> $label" >> "$LOG"
    bcftools isec --threads $THREADS -n "$n" -p "$d" "${files[@]}" 2>>"$LOG"
    local c; c=$(wc -l < "$d/sites.txt" 2>/dev/null)
    printf "%s\t%s\t%s\n" "$group" "$label" "${c:-ERR}" >> "$SUM"
    echo "[$(date +%H:%M:%S)] DONE $group $label sites=${c:-ERR}" >> "$LOG"
  done
}

do_group normal \
  "$NSNP/BCF_normal_SNPs.vcf.gz" "$NSNP/Deepvariant_normal_SNPs.vcf.gz" \
  "$NSNP/Freebayes_normal_SNPs.vcf.gz" "$NSNP/Gatk_normal_SNPs.vcf.gz" \
  "$NSNP/Varscan_normal_SNPs.vcf.gz" &

do_group abnormal \
  "$NSNP/BCF_abnormal_SNPs.vcf.gz" "$NSNP/Deepvariant_abnormal_SNPs.vcf.gz" \
  "$NSNP/Freebayes_abnormal_SNPs.vcf.gz" "$NSNP/Gatk_abnormal_SNPs.vcf.gz" \
  "$NSNP/Varscan_abnormal_SNPs.vcf.gz" &
wait
echo "[$(date +%H:%M:%S)] ALL DONE" >> "$LOG"
sort "$SUM" -o "$SUM"
