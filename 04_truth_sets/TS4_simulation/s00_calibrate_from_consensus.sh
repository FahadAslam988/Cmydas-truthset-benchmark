#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# ============================================================================
# STAGE 0c — RECALIBRATE from a 5-tool CONSENSUS, for BOTH samples.
# Replaces the weak 2-tool/normal-only calibration. For each sample we merge the
# 5 callers, keep SNVs called by >=3/5 tools (consensus removes single-caller FPs),
# and measure density, Ti/Tv, and het:hom (majority genotype). Read-only on inputs.
# ============================================================================
set -uo pipefail
RERUN="${WORK_DIR}"
BASE="$RERUN/09_simulation_truthset/RUN_DeepVariant_fair_benchmark_simuG-ART"
NSNP_DIR="$RERUN/04_norm_snps"
PARAMS="$BASE/params.txt"; LOG="$BASE/run.log"
WORK="$BASE/outputs/recalib"; mkdir -p "$WORK"
BCFTOOLS=bcftools
GENOME_BP=2134000000
SUBSET_BP=$(grep -E '^SUBSET_BP=' "$PARAMS" | head -1 | cut -d= -f2 | awk '{print $1}')
TOOLS=(BCF Deepvariant Freebayes Gatk Varscan)
MINK=3          # >=3 of 5 tools = consensus
log(){ echo "[$(date '+%F %T')] [s00c] $*" | tee -a "$LOG"; }
log "===== STAGE 0c recalibrate (5-tool consensus, both samples) START ====="

# measure one sample -> prints: SAMPLE nSNV ts tv het hom
measure_sample(){
  local samp="$1"; local files=()
  for t in "${TOOLS[@]}"; do files+=("$NSNP_DIR/${t}_${samp}_SNPs.vcf.gz"); done
  "$BCFTOOLS" merge --force-samples -m none "${files[@]}" -Oz -o "$WORK/${samp}_merged.vcf.gz" 2>>"$LOG"
  "$BCFTOOLS" index -f "$WORK/${samp}_merged.vcf.gz" 2>>"$LOG"
  # per-site: 5 GTs; keep biallelic SNV called by >=MINK tools; majority het/hom; Ti/Tv
  "$BCFTOOLS" query -f '%REF\t%ALT[\t%GT]\n' "$WORK/${samp}_merged.vcf.gz" 2>>"$LOG" | awk -v mink="$MINK" '
    function isTs(r,a){ return ((r=="A"&&a=="G")||(r=="G"&&a=="A")||(r=="C"&&a=="T")||(r=="T"&&a=="C")) }
    {
      ref=$1; alt=$2; if(length(ref)!=1||length(alt)!=1) next;        # SNV only
      called=0; het=0; hom=0;
      for(i=3;i<=NF;i++){ g=$i; gsub(/\|/,"/",g);
        if(g ~ /1/){ called++; if(g=="1/1") hom++; else het++ } }
      if(called<mink) next;
      nSNV++;
      if(isTs(toupper(ref),toupper(alt))) ts++; else tv++;
      if(het>=hom) sHet++; else sHom++;
    }
    END{ printf "%d %d %d %d %d\n", nSNV, ts, tv, sHet, sHom }'
}

TOT_SNV=0; TOT_TS=0; TOT_TV=0; TOT_HET=0; TOT_HOM=0
for s in normal abnormal; do
  log "measuring $s (5-tool >=${MINK} consensus)..."
  read nsnv ts tv het hom < <(measure_sample "$s")
  dens=$(awk -v n=$nsnv -v g=$GENOME_BP 'BEGIN{printf "%.6f", n/g}')
  titv=$(awk -v a=$ts -v b=$tv 'BEGIN{if(b>0)printf "%.2f",a/b; else print "NA"}')
  hf=$(awk -v h=$het -v m=$hom 'BEGIN{if(h+m>0)printf "%.3f",h/(h+m); else print "NA"}')
  log "  $s consensus: SNV=$nsnv density=$dens (1/$(awk -v d=$dens 'BEGIN{printf "%.0f",1/d}')bp) Ti/Tv=$titv het_frac=$hf"
  TOT_SNV=$((TOT_SNV+nsnv)); TOT_TS=$((TOT_TS+ts)); TOT_TV=$((TOT_TV+tv)); TOT_HET=$((TOT_HET+het)); TOT_HOM=$((TOT_HOM+hom))
done

# combined (pooled both samples)
DENS=$(awk -v n=$TOT_SNV -v g=$((2*GENOME_BP)) 'BEGIN{printf "%.6f", n/g}')   # avg over 2 samples
TITV=$(awk -v a=$TOT_TS -v b=$TOT_TV 'BEGIN{printf "%.2f",a/b}')
HETF=$(awk -v h=$TOT_HET -v m=$TOT_HOM 'BEGIN{printf "%.3f",h/(h+m)}')
NEW_NSNP=$(awk -v d=$DENS -v s=$SUBSET_BP 'BEGIN{printf "%d", d*s}')
NEW_NIND=$(awk -v n=$NEW_NSNP 'BEGIN{printf "%d", n/8}')
log "COMBINED (5-tool consensus, both samples): density=$DENS (1/$(awk -v d=$DENS 'BEGIN{printf "%.0f",1/d}')bp) Ti/Tv=$TITV het_frac=$HETF"
log "-> subset counts: N_SNP=$NEW_NSNP N_INDEL=$NEW_NIND (8:1)"

{
  echo "# ==== STAGE 0c RECALIBRATION (5-tool >=3 consensus, normal+abnormal) $(date '+%F %T') ===="
  echo "# supersedes the earlier DV+GATK/normal-only values above."
  echo "RECAL_TITV=$TITV"
  echo "RECAL_HETF=$HETF"
  echo "RECAL_DENSITY=$DENS"
  echo "RECAL_N_SNP=$NEW_NSNP"
  echo "RECAL_N_INDEL=$NEW_NIND"
  echo "# STAGE0c_DONE $(date '+%F %T')"
} >> "$PARAMS"
rm -f "$WORK"/*_merged.vcf.gz* ; rmdir "$WORK" 2>/dev/null || true
log "===== STAGE 0c DONE (values in params.txt as RECAL_*) ====="
