#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# ============================================================================
# checkpoint); only inputs/outputs point at the simulation. CPU image :latest
# (= v1.9.0, same as the real stock run). num_shards=80 (speed only; identical
# results). Runs via docker WITHOUT sudo. Resumable + logged.
# ============================================================================
set -uo pipefail
BASE="${SIM_DIR}"
REFDIR="${REF_DIR}"
BAMDIR="$BASE/outputs/03_aligned"
OUTDIR="$BASE/outputs/04_calls"
IMG="google/deepvariant:latest"          # = v1.9.0 CPU (same as real stock run)
LOG="$BASE/run.log"; PARAMS="$BASE/params.txt"
UG="$(id -u):$(id -g)"
mkdir -p "$OUTDIR"
log(){ echo "[$(date '+%F %T')] [s06a] $*" | tee -a "$LOG"; }

run_dv(){   # $1=label  $2=outname  $3=extra_docker(model mount or "")  $4=extra_arg(customized_model or "")
  local label="$1" out="$2" mnt="$3" marg="$4"
  if [ -s "$OUTDIR/${out}.vcf.gz" ]; then log "$label already done ($out.vcf.gz) — skip"; return 0; fi
  log "running $label (CPU, $IMG, num_shards=80)..."
  docker run --rm --user "$UG" \
    -v "$REFDIR":/ref:ro -v "$BAMDIR":/bam:ro -v "$OUTDIR":/output $mnt \
    "$IMG" /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS $marg \
    --ref=/ref/Reference.fasta \
    --reads=/bam/sim.dedup.bam \
    --output_vcf=/output/${out}.vcf.gz \
    --output_gvcf=/output/${out}.g.vcf.gz \
    --num_shards=80 \
    --intermediate_results_dir=/output/intermediates_${out} \
    --make_examples_extra_args="small_model_call_multiallelics=false" \
    2>&1 | tail -20 | tee -a "$LOG"
  if [ -s "$OUTDIR/${out}.vcf.gz" ]; then
    local n; n=$(bcftools view -H "$OUTDIR/${out}.vcf.gz" 2>/dev/null | wc -l)
    log "$label DONE -> $out.vcf.gz ($n records)"
    rm -rf "$OUTDIR/intermediates_${out}"
  else
    log "ERROR: $label produced no VCF"; return 3
  fi
}

run_dv "StockDV"     "StockDV_sim"     ""                          ""
{ echo "# ---- STAGE 5a DeepVariant ----"
  echo "DV_STOCK_VCF=$OUTDIR/StockDV_sim.vcf.gz"
  echo "# STAGE5a_DONE $(date '+%F %T')"; } >> "$PARAMS"
log "===== STAGE 5a DONE ====="
