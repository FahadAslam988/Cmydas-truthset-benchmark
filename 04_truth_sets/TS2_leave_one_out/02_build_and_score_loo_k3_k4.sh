#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# LOO threshold SENSITIVITY sweep (research-standard k = 2, 3, 4 of the other 4 callers).
# k=2 is reused from the completed main LOO (06_loo). This script builds & benchmarks k=3 and k=4.
# For each tested caller X: truth = >=k of the OTHER 4 callers; query = X's own normalized SNPs.
# Matching = allele-level (hap.py --engine=vcfeval --set-gt hom); 0/0 & ./. removed first.
# Idempotent: a run with a valid summary.csv is skipped. Tuned for 80 threads / 640 GB RAM.
set -uo pipefail

THREADS=8                 # per hap.py job
MAXPAR=10                 # 10 x 8 = 80 threads (full machine); ~150 GB RAM peak
IMG="jmcdani20/hap.py:v0.3.12"
FAHAD="${DATA_ROOT}/Fahad"
REFDIR="${REF_DIR}"
REF_FA_C="/ref/Reference.fasta"; REF_SDF_C="/ref/Reference.sdf"
RERUN="$FAHAD/03_Corrected_Work/rerun"
NSNP="$RERUN/04_norm_snps"; RH="$RERUN/04b_varscan_reheader"
SENS="$RERUN/07_loo_sensitivity"
LOG="$SENS/summary/sensitivity_log.txt"; : > "$LOG"
log(){ echo "[$(date '+%H:%M:%S')] $*" >> "$LOG"; }

CALLERS=(BCF Deepvariant Freebayes Gatk Varscan)
DSETS=(normal abnormal)
KS=(3 4)                  # k=2 already reused from 06_loo

kdir(){ local k="$1"; case $k in 3) echo "$SENS/k3_ge3_of4";; 4) echo "$SENS/k4_eq4_of4";; esac; }
vcf_host(){ local c="$1" g="$2"; if [ "$c" = "Varscan" ]; then echo "$RH/Varscan_${g}.vcf.gz"; else echo "$NSNP/${c}_${g}_SNPs.vcf.gz"; fi; }
to_c(){ echo "/data/fahad/${1#$FAHAD/}"; }

run_one(){
  local k="$1" group="$2" test="$3"
  local d; d="$(kdir "$k")"
  local pfx="$d/happy/${group}_${test}"
  local sc="$pfx.summary.csv"
  if [ -f "$sc" ] && awk -F',' 'NR>1 && $1=="SNP" && $2=="ALL"{ok=1} END{exit ok?0:1}' "$sc"; then
    log "SKIP (done) k=$k $group wo=$test"; emit_res "$k" "$group" "$test" "$sc"; return 0
  fi
  # build truth = >=k of the OTHER 4 callers
  local others=() c; for c in "${CALLERS[@]}"; do [ "$c" != "$test" ] && others+=("$(vcf_host "$c" "$group")"); done
  local isec="$d/truth/isec_${group}_wo_${test}"; rm -rf "$isec"; mkdir -p "$isec"
  bcftools isec -n +$k -O z -p "$isec" "${others[@]}" 2>>"$LOG"
  local parts=( "$isec"/0*.vcf.gz ) p; for p in "${parts[@]}"; do bcftools index -f "$p" 2>>"$LOG"; done
  local truth="$d/truth/truth_${group}_wo_${test}.vcf.gz"
  bcftools concat -a "${parts[@]}" -Ou 2>>"$LOG" | bcftools sort -Ou 2>>"$LOG" | bcftools norm -d exact -Oz -o "$truth" 2>>"$LOG"
  bcftools index -f -t "$truth" 2>>"$LOG"
  rm -rf "$isec"    # keep only the merged truth (save disk)
  # allele-only: drop 0/0 & ./. from truth and query
  local th="$d/prep/truth_${group}_wo_${test}.altonly.vcf.gz"
  local qh="$d/prep/calls_${group}_${test}.altonly.vcf.gz"
  bcftools view -i 'GT="alt"' "$truth" -Oz -o "$th" 2>>"$LOG" && bcftools index -f -t "$th" 2>>"$LOG"
  bcftools view -i 'GT="alt"' "$(vcf_host "$test" "$group")" -Oz -o "$qh" 2>>"$LOG" && bcftools index -f -t "$qh" 2>>"$LOG"
  log "START hap.py k=$k $group wo=$test"
  docker run --rm -u "$(id -u):$(id -g)" -e HOME=/data/fahad -v "$FAHAD":/data/fahad -v "$REFDIR":/ref \
    "$IMG" /opt/hap.py/bin/hap.py "$(to_c "$th")" "$(to_c "$qh")" -r "$REF_FA_C" -o "$(to_c "$pfx")" \
    --engine=vcfeval --engine-vcfeval-template "$REF_SDF_C" --set-gt hom --threads "$THREADS" \
    --scratch-prefix "$(to_c "$d/prep")" >>"$LOG" 2>&1
  log "END   hap.py k=$k $group wo=$test rc=$?"
  rm -f "$th" "$qh" "${th}.tbi" "${qh}.tbi"    # prune big prep files after use
  emit_res "$k" "$group" "$test" "$sc"
}

emit_res(){   # k group test summary.csv -> one row in summary dir
  local k="$1" group="$2" test="$3" sc="$4"
  [ -f "$sc" ] || { printf "%s\t%s\t%s\tNA\tNA\tNA\tNA\tNA\n" "$k" "$group" "$test" > "$SENS/summary/.res_k${k}_${group}_${test}.tsv"; return; }
  awk -F',' -v k="$k" -v g="$group" -v t="$test" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i}
    $h["Type"]=="SNP"&&$h["Filter"]=="ALL"{printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",k,g,t,
      $h["TRUTH.TP"],$h["QUERY.FP"],$h["TRUTH.FN"],$h["METRIC.Precision"],$h["METRIC.Recall"],$h["METRIC.F1_Score"]}' \
    "$sc" > "$SENS/summary/.res_k${k}_${group}_${test}.tsv"
}

# k=2 rows from the reused results
for g in "${DSETS[@]}"; do for t in "${CALLERS[@]}"; do
  emit_res 2 "$g" "$t" "$SENS/k2_ge2_of4/happy/${g}_${t}.summary.csv"
done; done

# k=3 and k=4 runs, throttled to MAXPAR
for k in "${KS[@]}"; do for g in "${DSETS[@]}"; do for t in "${CALLERS[@]}"; do
  run_one "$k" "$g" "$t" &
  while [ "$(jobs -rp | wc -l)" -ge "$MAXPAR" ]; do wait -n; done
done; done; done
wait

# combined raw table
SUM="$SENS/summary/sensitivity_ranking.tsv"
printf "k\tGroup\tTestedCaller\tTP\tFP\tFN\tPrecision\tRecall\tF1\n" > "$SUM"
cat "$SENS/summary/".res_k*.tsv 2>/dev/null | sort -k1,1n -k2,2 -k3,3 >> "$SUM"
rm -f "$SENS/summary/".res_k*.tsv
log "ALL DONE"
echo "DONE -> $SUM"
