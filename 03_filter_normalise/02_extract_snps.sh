#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# Stage 2 re-run: extract SNPs from filtered VCFs, compare to documented SNP counts.
# Command: bcftools view -v snps (same as original make_snps.sh)
set -uo pipefail
THREADS=8

IN="${WORK_DIR}/01_filtered"
OUT="${WORK_DIR}/02_snps"
mkdir -p "$OUT"
LOG="$OUT/snps_log.txt"; : > "$LOG"
rm -f "$OUT"/.count_*.tsv

run() {
  local name="$1" target="$2"
  local infile="$IN/${name}.filtered.vcf.gz"
  local out="$OUT/${name}_SNPs.vcf.gz"
  echo "[$(date +%H:%M:%S)] START $name" >> "$LOG"
  if [ ! -f "$infile" ]; then echo "[ERR] missing $infile" >> "$LOG"; return; fi
  bcftools view --threads $THREADS -v snps "$infile" -Oz -o "$out" 2>>"$LOG"
  bcftools index --threads $THREADS "$out" 2>>"$LOG"
  local n; n=$(bcftools index -n "$out" 2>>"$LOG")
  printf "%s\t%s\t%s\n" "$name" "${n:-ERROR}" "$target" > "$OUT/.count_${name}.tsv"
  echo "[$(date +%H:%M:%S)] DONE  $name count=${n:-ERROR} target=$target" >> "$LOG"
}

run BCF_normal          12678209 &
run BCF_abnormal        12726145 &
run Deepvariant_normal  12058893 &
run Deepvariant_abnormal 12103451 &
run Freebayes_normal    11247911 &
run Freebayes_abnormal  11300566 &
run Gatk_normal         12629773 &
run Gatk_abnormal       12663916 &
run Varscan_normal      12228987 &
run Varscan_abnormal    12288474 &
wait

SUMMARY="$OUT/snp_counts.tsv"
printf "Callset\tNewCount\tTarget\tMatch\n" > "$SUMMARY"
cat "$OUT"/.count_*.tsv 2>/dev/null | while IFS=$'\t' read -r name n target; do
  if [ "$n" = "$target" ]; then m="YES"; else m="NO(diff=$(( ${n:-0} - target )))"; fi
  printf "%s\t%s\t%s\t%s\n" "$name" "$n" "$target" "$m" >> "$SUMMARY"
done
rm -f "$OUT"/.count_*.tsv
echo "[$(date +%H:%M:%S)] ALL DONE" >> "$LOG"
