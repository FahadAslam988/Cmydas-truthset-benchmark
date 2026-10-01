#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# ============================================================================
# STAGE 6 — Benchmark all 6 callers vs the known truth (hap.py, default engine).
# For each caller: apply the SAME per-tool filter as the real pipeline
# (PASS/. + alt-depth>3 + VAF>2%, adapted to each tool's fields), normalize,
# then hap.py (genotype-aware) restricted to callable.bed. Collate SNP
# Precision/Recall/F1 -> 06_results/simulation_ranking.tsv. Guarded + logged.
# ============================================================================
set -uo pipefail
BASE="${SIM_DIR}"
REFDIR="${REF_DIR}"
REF="$REFDIR/Reference.fasta"
GDIR="$BASE/outputs/01_mutated_genome"
CALLS="$BASE/outputs/04_calls"
BENCH="$BASE/outputs/05_benchmark"
RES="$BASE/outputs/06_results"
TRUTH="$GDIR/truth.vcf.gz"; CALLABLE="$GDIR/callable.bed"
BCFT=bcftools
HAPIMG="jmcdani20/hap.py:v0.3.12"
UG="$(id -u):$(id -g)"; LOG="$BASE/run.log"
mkdir -p "$BENCH" "$RES"
log(){ echo "[$(date '+%F %T')] [s07] $*" | tee -a "$LOG"; }

# per-tool filter expressions (verbatim from real pipeline stage1_filter.sh)
EXPR_STD='(FILTER="PASS" || FILTER=".") && FMT/AD[0:1] > 3 && FMT/AD[0:1] / FMT/DP[0] > 0.02'
EXPR_BCF='(FILTER="PASS" || FILTER=".") && FMT/AD[0:1] > 3 && FMT/AD[0:1] / INFO/DP > 0.02'
EXPR_VAR='(FILTER="PASS" || FILTER=".") && FMT/AD[0]   > 3 && FMT/AD[0]   / FMT/DP[0] > 0.02'

# caller -> filter-expr map
declare -A EXPR=( [StockDV]="$EXPR_STD" [Gatk]="$EXPR_STD"
                  [Freebayes]="$EXPR_STD" [BCF]="$EXPR_BCF" [Varscan]="$EXPR_VAR" )
ORDER=(StockDV Gatk BCF Freebayes Varscan)

log "===== STAGE 6 benchmark START (filter -> normalize -> hap.py) ====="

# --- filter + normalize each caller ----------------------------------------
for c in "${ORDER[@]}"; do
  in="$CALLS/${c}_sim.vcf.gz"; out="$BENCH/${c}.filt.vcf.gz"
  [ -s "$out" ] && { log "$c filtered exists, skip"; continue; }
  raw=$("$BCFT" view -H "$in" 2>/dev/null | wc -l)
  "$BCFT" view -i "${EXPR[$c]}" "$in" 2>>"$LOG" \
    | "$BCFT" norm -f "$REF" -m -any 2>>"$LOG" \
    | "$BCFT" sort -Oz -o "$out" 2>>"$LOG"
  "$BCFT" index -t "$out" 2>>"$LOG"
  kept=$("$BCFT" view -H "$out" 2>/dev/null | wc -l)
  log "$c: filtered $raw -> $kept records"
  [ "$kept" -eq 0 ] && log "WARNING: $c filtered to 0 — check its AD/DP fields!"
done

# --- hap.py each caller vs truth (default engine, genotype-aware) -----------
for c in "${ORDER[@]}"; do
  q="$BENCH/${c}.filt.vcf.gz"
  [ -s "$BENCH/${c}.summary.csv" ] && { log "$c hap.py exists, skip"; continue; }
  [ -s "$q" ] || { log "$c: no filtered VCF, skip hap.py"; continue; }
  log "hap.py: $c ..."
  docker run --rm -u "$UG" -v "$BASE/outputs":/data -v "$REFDIR":/ref "$HAPIMG" \
    /opt/hap.py/bin/hap.py \
    /data/01_mutated_genome/truth.vcf.gz \
    /data/05_benchmark/${c}.filt.vcf.gz \
    -r /ref/Reference.fasta \
    -T /data/01_mutated_genome/callable.bed \
    -o /data/05_benchmark/${c} \
    --threads 16 2>&1 | tail -3 | tee -a "$LOG" || log "!!! hap.py FAILED for $c (continuing)"
done

# --- collate SNP metrics into the ranking table ----------------------------
RANK="$RES/simulation_ranking.tsv"
printf "Caller\tType\tTRUTH.TOTAL\tTP\tFN\tFP\tRecall\tPrecision\tF1\n" > "$RANK"
for c in "${ORDER[@]}"; do
  s="$BENCH/${c}.summary.csv"; [ -s "$s" ] || { printf "%s\tSNP\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n" "$c" >> "$RANK"; continue; }
  # summary.csv columns: Type,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.FP,...,METRIC.Recall,METRIC.Precision,...,METRIC.F1_Score
  awk -F',' -v C="$c" '
    NR==1{for(i=1;i<=NF;i++) h[$i]=i; next}
    $(h["Type"])=="SNP" && $(h["Filter"])=="PASS"{
      printf "%s\tSNP\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", C,
        $(h["TRUTH.TOTAL"]),$(h["TRUTH.TP"]),$(h["TRUTH.FN"]),$(h["QUERY.FP"]),
        $(h["METRIC.Recall"]),$(h["METRIC.Precision"]),$(h["METRIC.F1_Score"]) }' "$s" >> "$RANK" 2>/dev/null
done
log "===== STAGE 6 DONE -> $RANK ====="
column -t -s$'\t' "$RANK" | tee -a "$LOG"
