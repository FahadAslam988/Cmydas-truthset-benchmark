#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# Stage 5 (research standard): LEAVE-ONE-OUT benchmarking to remove circularity (Issue B1).
# For each tested caller X: truth = consensus of the OTHER 4 callers (>= K agree); score X with RTG vcfeval.
# Refs: Olson 2023 (Nat Rev Genet); Krusche 2019 (Nat Biotechnol); Zook 2019 (GIAB).
set -uo pipefail
THREADS=6
K=2   # independent support: >= 2 of the other 4 callers
REF_SDF="${REF_SDF}"
NSNP="${WORK_DIR}/04_norm_snps"
OUT="${WORK_DIR}/06_loo"
mkdir -p "$OUT"
LOG="$OUT/loo_log.txt"; : > "$LOG"
rm -f "$OUT"/.res_*.tsv

CALLERS=(BCF Deepvariant Freebayes Gatk Varscan)

run_loo() {
  local group="$1" test="$2"
  local others=()
  local c
  for c in "${CALLERS[@]}"; do [ "$c" != "$test" ] && others+=("$NSNP/${c}_${group}_SNPs.vcf.gz"); done
  local d="$OUT/${group}_leaveout_${test}"
  rm -rf "$d"; mkdir -p "$d"
  echo "[$(date +%H:%M:%S)] START $group leave-out=$test (truth = other 4, >=$K)" >> "$LOG"

  # 1) consensus of the other 4 at >=K  (compressed isec outputs)
  bcftools isec -n +$K -O z -p "$d/isec" "${others[@]}" 2>>"$LOG"
  # 2) merge per-input outputs into ONE truth VCF (dedup on POS+REF+ALT)
  local parts=( "$d"/isec/0*.vcf.gz )
  for p in "${parts[@]}"; do bcftools index -f "$p" 2>>"$LOG"; done
  local truth="$d/truth_${group}_wo_${test}.vcf.gz"
  bcftools concat -a "${parts[@]}" -Ou 2>>"$LOG" | bcftools sort -Ou 2>>"$LOG" | bcftools norm -d exact -Oz -o "$truth" 2>>"$LOG"
  bcftools index -f "$truth" 2>>"$LOG"
  local tsites; tsites=$(bcftools index -n "$truth" 2>>"$LOG")

  # 3) vcfeval: score the tested caller against that independent truth set
  local calls="$NSNP/${test}_${group}_SNPs.vcf.gz"
  local ev="$d/vcfeval"; rm -rf "$ev"
  rtg vcfeval -b "$truth" -c "$calls" -t "$REF_SDF" -o "$ev" --vcf-score-field QUAL -T $THREADS >>"$LOG" 2>&1

  # 4) parse the summary (last numeric row = overall / no-threshold)
  local row prec sens f1 tpb fp fn
  row=$(grep -vE 'Threshold|^-{3,}|^[[:space:]]*$' "$ev/summary.txt" 2>/dev/null | tail -1)
  # cols: Threshold TP-base TP-call FP FN Precision Sensitivity F-measure
  tpb=$(awk '{print $2}' <<<"$row"); fp=$(awk '{print $4}' <<<"$row"); fn=$(awk '{print $5}' <<<"$row")
  prec=$(awk '{print $6}' <<<"$row"); sens=$(awk '{print $7}' <<<"$row"); f1=$(awk '{print $8}' <<<"$row")
  printf "%s\t%s\t>=%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$group" "$test" "$K" "${tsites:-NA}" "${tpb:-NA}" "${fp:-NA}" "${fn:-NA}" "${prec:-NA}" "${sens:-NA}" "${f1:-NA}" \
    > "$OUT/.res_${group}_${test}.tsv"
  echo "[$(date +%H:%M:%S)] DONE  $group leave-out=$test  F1=${f1:-NA} (truth sites=$tsites)" >> "$LOG"
}

for g in normal abnormal; do
  for t in "${CALLERS[@]}"; do
    run_loo "$g" "$t" &
  done
done
wait

SUM="$OUT/loo_vcfeval_summary.tsv"
printf "Group\tTestedCaller\tThreshold\tTruthSites\tTP\tFP\tFN\tPrecision\tRecall\tF1\n" > "$SUM"
cat "$OUT"/.res_*.tsv 2>/dev/null | sort >> "$SUM"
rm -f "$OUT"/.res_*.tsv
echo "[$(date +%H:%M:%S)] ALL DONE" >> "$LOG"
