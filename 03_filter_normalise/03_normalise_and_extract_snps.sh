#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# Stage 3 (research standard): NORMALIZE filtered VCFs, then re-extract SNPs.
# Standard: bcftools norm -f REF -m -any  (split multiallelics + left-align) -> required before cross-caller consensus.
# Refs: Tan 2015; Bayat 2017; Danecek 2021.
set -uo pipefail
THREADS=8
REF="${REF_FASTA}"

FILT="${WORK_DIR}/01_filtered"
NORMDIR="${WORK_DIR}/03_normalized"
SNPDIR="${WORK_DIR}/04_norm_snps"
mkdir -p "$NORMDIR" "$SNPDIR"
LOG="$NORMDIR/normalize_log.txt"; : > "$LOG"
rm -f "$SNPDIR"/.count_*.tsv

run() {
  local name="$1"
  local infile="$FILT/${name}.filtered.vcf.gz"
  local norm="$NORMDIR/${name}.norm.vcf.gz"
  local snp="$SNPDIR/${name}_SNPs.vcf.gz"
  echo "[$(date +%H:%M:%S)] START $name" >> "$LOG"
  if [ ! -f "$infile" ]; then echo "[ERR] missing $infile" >> "$LOG"; return; fi
  # 1) normalize: split multiallelics (-m -any), left-align, warn (not fail) on ref mismatch
  bcftools norm -f "$REF" -m -any -c w --threads $THREADS "$infile" -Oz -o "$norm" 2>>"$LOG"
  bcftools index -f --threads $THREADS "$norm" 2>>"$LOG"
  # 2) re-extract SNPs from normalized
  bcftools view --threads $THREADS -v snps "$norm" -Oz -o "$snp" 2>>"$LOG"
  bcftools index -f --threads $THREADS "$snp" 2>>"$LOG"
  local n; n=$(bcftools index -n "$snp" 2>>"$LOG")
  local multi; multi=$(bcftools view -H -m3 "$snp" 2>/dev/null | wc -l)
  printf "%s\t%s\t%s\n" "$name" "${n:-ERR}" "${multi}" > "$SNPDIR/.count_${name}.tsv"
  echo "[$(date +%H:%M:%S)] DONE  $name normSNPs=${n:-ERR} remaining_multiallelic=${multi}" >> "$LOG"
}

for s in BCF_normal BCF_abnormal Deepvariant_normal Deepvariant_abnormal \
         Freebayes_normal Freebayes_abnormal Gatk_normal Gatk_abnormal \
         Varscan_normal Varscan_abnormal; do
  run "$s" &
done
wait

SUM="$SNPDIR/norm_snp_counts.tsv"
printf "Callset\tNormalized_SNPs\tRemaining_Multiallelic\n" > "$SUM"
cat "$SNPDIR"/.count_*.tsv 2>/dev/null | sort >> "$SUM"
rm -f "$SNPDIR"/.count_*.tsv
echo "[$(date +%H:%M:%S)] ALL DONE" >> "$LOG"
