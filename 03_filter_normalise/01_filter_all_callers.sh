#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# Stage 1 re-run: reproduce FILTERING (same filter as original), compare to documented counts.
# Filter rule (all tools): FILTER=PASS/.  AND  alt-depth > 3  AND  VAF > 0.02
# Per-tool expression adapted to each tool's VCF fields (verified 2026-07-12).
set -uo pipefail

THREADS=8   # per job; 10 jobs x 8 = 80 threads

A="${CALLS_DIR}"
B="${CALLS_DIR}"
OUT="${WORK_DIR}/01_filtered"
mkdir -p "$OUT"
LOG="$OUT/filter_log.txt"
: > "$LOG"
rm -f "$OUT"/.count_*.tsv

EXPR_STD='(FILTER="PASS" || FILTER=".") && FMT/AD[0:1] > 3 && FMT/AD[0:1] / FMT/DP[0] > 0.02'
EXPR_BCF='(FILTER="PASS" || FILTER=".") && FMT/AD[0:1] > 3 && FMT/AD[0:1] / INFO/DP > 0.02'
EXPR_VAR='(FILTER="PASS" || FILTER=".") && FMT/AD[0]   > 3 && FMT/AD[0]   / FMT/DP[0] > 0.02'

run() {
  local tool="$1" sample="$2" infile="$3" expr="$4" target="$5"
  local out="$OUT/${tool}_${sample}.filtered.vcf.gz"
  echo "[$(date +%H:%M:%S)] START $tool $sample" >> "$LOG"
  if [ ! -f "$infile" ]; then echo "[ERR] missing input: $infile" >> "$LOG"; return; fi
  bcftools view --threads $THREADS -i "$expr" "$infile" -Oz -o "$out" 2>>"$LOG"
  bcftools index --threads $THREADS "$out" 2>>"$LOG"
  local n; n=$(bcftools index -n "$out" 2>>"$LOG")
  printf "%s\t%s\t%s\t%s\n" "$tool" "$sample" "${n:-ERROR}" "$target" > "$OUT/.count_${tool}_${sample}.tsv"
  echo "[$(date +%H:%M:%S)] DONE  $tool $sample  count=${n:-ERROR} target=$target" >> "$LOG"
}

run BCF         normal   "$A/bcftools_out/normal.vcf.gz"                      "$EXPR_BCF" 14220913 &
run BCF         abnormal "$A/bcftools_out/abnormal.vcf.gz"                    "$EXPR_BCF" 14285440 &
run Deepvariant normal   "$A/DeepVaraints/CH-NORMS1_deepvariant_1_9.vcf.gz"  "$EXPR_STD" 13705037 &
run Deepvariant abnormal "$A/DeepVaraints/Ab-NormS2_deepvariant_1_9.vcf.gz"  "$EXPR_STD" 13767503 &
run Freebayes   normal   "$A/freebayes_out/normal.freebayes.vcf.gz"          "$EXPR_STD" 13239802 &
run Freebayes   abnormal "$A/freebayes_out/abnormal.freebayes.vcf.gz"        "$EXPR_STD" 13332001 &
run Gatk        normal   "$B/gatk_out/normal.gatk.genotyped.vcf.gz"          "$EXPR_STD" 14467559 &
run Gatk        abnormal "$B/gatk_out/abnormal.gatk.genotyped.vcf.gz"        "$EXPR_STD" 14522966 &
run Varscan     normal   "$A/varscan_out/normal.varscan.vcf.gz"              "$EXPR_VAR" 13753976 &
run Varscan     abnormal "$A/varscan_out/abnormal.varscan.vcf.gz"            "$EXPR_VAR" 13833482 &
wait

SUMMARY="$OUT/filter_counts.tsv"
printf "Tool\tSample\tNewCount\tTarget\tMatch\n" > "$SUMMARY"
cat "$OUT"/.count_*.tsv 2>/dev/null | while IFS=$'\t' read -r tool sample n target; do
  if [ "$n" = "$target" ]; then m="YES"; else m="NO(diff=$(( ${n:-0} - target )))"; fi
  printf "%s\t%s\t%s\t%s\t%s\n" "$tool" "$sample" "$n" "$target" "$m" >> "$SUMMARY"
done
rm -f "$OUT"/.count_*.tsv
echo "[$(date +%H:%M:%S)] ALL DONE" >> "$LOG"
