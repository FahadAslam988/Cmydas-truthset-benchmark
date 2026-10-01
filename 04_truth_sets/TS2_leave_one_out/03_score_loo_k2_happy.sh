#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# Stage 5 (hap.py version): Leave-One-Out benchmark, GA4GH-standard hap.py + vcfeval engine.
# Truth sets already built & integrity-checked (06_loo/<group>_leaveout_<caller>/truth_*.vcf.gz),
# so we DO NOT rebuild them. For each tested caller:
#   truth  = >=2 of the OTHER 4 callers (already built)
#   query  = that caller's own normalized SNPs (VarScan -> reheadered copy)
# Matching = ALLELE-LEVEL: hap.py --engine=vcfeval --set-gt hom forces all GTs to 1/1 in BOTH
#   truth and query, so a match depends only on the allele (genotype cannot cause a mismatch).
#   => same intent as rtg vcfeval --squash-ploidy, but reported as GA4GH-standard hap.py output.
# 0/0 (hom-ref) and ./. records are removed FIRST, else --set-gt hom would turn them into
#   spurious positives.
#
# Usage:
#   ./stage5_leaveoneout_happy.sh            # full run: all 10 (2 groups x 5 callers), 3 at a time
#   ./stage5_leaveoneout_happy.sh <contig>   # VALIDATION: only normal/BCF, restricted to <contig>
set -uo pipefail

# ---- config ----
THREADS=8                 # threads per hap.py job
MAXPAR=8                  # run up to 8 at once (8x8=64 threads < 80; ample RAM). Idempotent: done runs skipped.
IMG="jmcdani20/hap.py:v0.3.12"
FAHAD="${DATA_ROOT}/Fahad"
REFDIR="${REF_DIR}"
REF_FA_C="/ref/Reference.fasta"          # path INSIDE container
REF_SDF_C="/ref/Reference.sdf"           # path INSIDE container
RERUN="$FAHAD/03_Corrected_Work/rerun"
NSNP="$RERUN/04_norm_snps"
RH="$RERUN/04b_varscan_reheader"
OUT="$RERUN/06_loo"
HAPPY="$OUT/happy"                        # all hap.py outputs land here
mkdir -p "$HAPPY"
LOG="$HAPPY/happy_log.txt"; : > "$LOG"

REGION="${1:-}"                           # optional single contig => validation mode

CALLERS=(BCF Deepvariant Freebayes Gatk Varscan)
DSETS=(normal abnormal)

# host path to a caller's normalized SNPs (VarScan -> reheadered copy, others -> original)
vcf_host() { local c="$1" g="$2"; if [ "$c" = "Varscan" ]; then echo "$RH/Varscan_${g}.vcf.gz"; else echo "$NSNP/${c}_${g}_SNPs.vcf.gz"; fi; }
# translate a host path under $FAHAD to its in-container path (/data/fahad/...)
to_c() { echo "/data/fahad/${1#$FAHAD/}"; }

run_one() {
  local group="$1" test="$2"
  local pfx_done="$HAPPY/${group}_${test}.summary.csv"
  # IDEMPOTENT: if this run already produced a valid SNP/ALL summary row, skip it (don't recompute)
  if [ -f "$pfx_done" ] && awk -F',' 'NR>1 && $1=="SNP" && $2=="ALL"{ok=1} END{exit ok?0:1}' "$pfx_done"; then
    echo "[$(date +%H:%M:%S)] SKIP (already done) $group leave-out=$test" >> "$LOG"
    # re-emit its result row so the final summary still includes it
    awk -F',' -v g="$group" -v t="$test" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} $h["Type"]=="SNP"&&$h["Filter"]=="ALL"{
      printf "%s\t%s\t>=2\t%s\t%s\t%s\t%s\t%s\t%s\n",g,t,$h["TRUTH.TP"],$h["QUERY.FP"],$h["TRUTH.FN"],$h["METRIC.Precision"],$h["METRIC.Recall"],$h["METRIC.F1_Score"]}' \
      "$pfx_done" > "$HAPPY/.res_${group}_${test}.tsv"
    return 0
  fi
  local d="$OUT/${group}_leaveout_${test}"
  local truth_raw="$d/truth_${group}_wo_${test}.vcf.gz"
  local calls_raw; calls_raw="$(vcf_host "$test" "$group")"
  local pfx_host="$HAPPY/${group}_${test}"
  mkdir -p "$HAPPY/prep"

  # 1) keep only records carrying an ALT allele (drop 0/0 and ./.) in BOTH truth and query
  local truth_h="$HAPPY/prep/truth_${group}_wo_${test}.altonly.vcf.gz"
  local calls_h="$HAPPY/prep/calls_${group}_${test}.altonly.vcf.gz"
  bcftools view -i 'GT="alt"' "$truth_raw"  -Oz -o "$truth_h" 2>>"$LOG" && bcftools index -f -t "$truth_h" 2>>"$LOG"
  bcftools view -i 'GT="alt"' "$calls_raw"  -Oz -o "$calls_h" 2>>"$LOG" && bcftools index -f -t "$calls_h" 2>>"$LOG"

  # 2) hap.py (GA4GH) driving the vcfeval engine, allele-level via --set-gt hom
  local loc_opt=(); [ -n "$REGION" ] && loc_opt=(--location "$REGION")
  echo "[$(date +%H:%M:%S)] START hap.py $group leave-out=$test ${REGION:+(region=$REGION)}" >> "$LOG"
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/data/fahad \
    -v "$FAHAD":/data/fahad -v "$REFDIR":/ref \
    "$IMG" /opt/hap.py/bin/hap.py \
      "$(to_c "$truth_h")" "$(to_c "$calls_h")" \
      -r "$REF_FA_C" -o "$(to_c "$pfx_host")" \
      --engine=vcfeval --engine-vcfeval-template "$REF_SDF_C" \
      --set-gt hom --threads "$THREADS" --scratch-prefix "$(to_c "$HAPPY/prep")" \
      "${loc_opt[@]}" >>"$LOG" 2>&1
  local rc=$?
  echo "[$(date +%H:%M:%S)] END   hap.py $group leave-out=$test rc=$rc" >> "$LOG"

  # 3) pull the SNP / Filter=ALL row (Recall, Precision, F1) from <prefix>.summary.csv
  local sc="$pfx_host.summary.csv" row rec prec f1 tp fp fn
  if [ -f "$sc" ]; then
    row=$(awk -F',' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["Type"]=="SNP" && $h["Filter"]=="ALL"{
            print $h["TRUTH.TP"]","$h["QUERY.FP"]","$h["TRUTH.FN"]","$h["METRIC.Precision"]","$h["METRIC.Recall"]","$h["METRIC.F1_Score"]}' "$sc")
    tp=${row%%,*}; rest=${row#*,}; fp=${rest%%,*}; rest=${rest#*,}; fn=${rest%%,*}
    rest=${rest#*,}; prec=${rest%%,*}; rest=${rest#*,}; rec=${rest%%,*}; f1=${rest#*,}
  fi
  printf "%s\t%s\t>=2\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$group" "$test" "${tp:-NA}" "${fp:-NA}" "${fn:-NA}" "${prec:-NA}" "${rec:-NA}" "${f1:-NA}" \
    > "$HAPPY/.res_${group}_${test}.tsv"
}

# ---- VALIDATION MODE: one caller/group on one contig, then stop ----
if [ -n "$REGION" ]; then
  echo "VALIDATION run: normal / BCF on contig $REGION" | tee -a "$LOG"
  run_one normal BCF
  echo "=== validation summary.csv (SNP/ALL) ==="; cat "$HAPPY/.res_normal_BCF.tsv"
  exit 0
fi

# ---- FULL RUN: all 10, throttled to $MAXPAR at a time ----
for g in "${DSETS[@]}"; do
  for t in "${CALLERS[@]}"; do
    run_one "$g" "$t" &
    while [ "$(jobs -rp | wc -l)" -ge "$MAXPAR" ]; do wait -n; done
  done
done
wait

SUM="$HAPPY/loo_happy_summary.tsv"
printf "Group\tTestedCaller\tThreshold\tTP\tFP\tFN\tPrecision\tRecall\tF1\n" > "$SUM"
cat "$HAPPY"/.res_*.tsv 2>/dev/null | sort >> "$SUM"
rm -f "$HAPPY"/.res_*.tsv
echo "[$(date +%H:%M:%S)] ALL DONE" >> "$LOG"
echo "=== FINAL hap.py LOO SUMMARY ==="; column -t "$SUM"
