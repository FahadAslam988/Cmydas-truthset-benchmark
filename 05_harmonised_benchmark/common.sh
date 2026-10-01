# Shared settings for the Paper 1 harmonised benchmark (sourced by every step script).
# ONE scoring method for every truth set, identical to the leave-one-out (TS2) benchmark
# (rerun/stage5_leaveoneout_happy.sh, rerun/07_loo_sensitivity/sensitivity_loo.sh):
#   hap.py v0.3.12, --engine=vcfeval (RTG 3.12.1), --set-gt hom (allele-level), 0/0 and ./. removed first.
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
DESK=${DATA_ROOT}
R=$DESK/Fahad/03_Corrected_Work/rerun
E=$R/10_Paper1_harmonised_benchmark_consensus_LOO_simulation
REFDIR=$DESK/DeepVariantTraining/ref                      # Reference.fasta (+.fai) and Reference.sdf
BCF=bcftools  # v1.22
IMG=jmcdani20/hap.py:v0.3.12
THREADS=8; MAXPAR=8
NSNP=$R/04_norm_snps                                      # SNPs from 03_normalized (filtered + normalised)
RH=$R/04b_varscan_reheader                                # VarScan copies used by the LOO hap.py runs
SIM=$R/09_simulation_truthset/RUN_DeepVariant_fair_benchmark_simuG-ART/outputs
CALLERS=(BCF Deepvariant Freebayes Gatk Varscan)
GROUPS_=(normal abnormal)
LOG=$E/run.log
mkdir -p "$E"
log(){ echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }
c(){ echo "/data/${1#$DESK/}"; }                         # host path -> path inside the hap.py container
query_snps(){ [ "$1" = Varscan ] && echo "$RH/Varscan_$2.vcf.gz" || echo "$NSNP/${1}_$2_SNPs.vcf.gz"; }
happy(){ # happy <truth> <query> <out-prefix> [extra args...]
  local t=$1 q=$2 o=$3; shift 3
  [ -s "$o.summary.csv" ] && { log "skip (exists) $(basename $o)"; return 0; }
  docker run --rm -u "$(id -u):$(id -g)" -e HOME="$(c $E)" -v "$DESK":/data "$IMG" /opt/hap.py/bin/hap.py \
    "$(c $t)" "$(c $q)" -r "$(c $REFDIR/Reference.fasta)" -o "$(c $o)" \
    --engine=vcfeval --engine-vcfeval-template "$(c $REFDIR/Reference.sdf)" --set-gt hom \
    --threads $THREADS --scratch-prefix "$(c $E/tmp)" "$@" >>"$E/happy_docker.log" 2>&1
  log "hap.py rc=$? $(basename $o)"; }
snp_row(){ # snp_row <summary.csv> -> TP FP FN Precision Recall F1
  awk -F',' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $h["Type"]=="SNP"&&$h["Filter"]=="ALL"{
    print $h["TRUTH.TP"]"\t"$h["QUERY.FP"]"\t"$h["TRUTH.FN"]"\t"$h["METRIC.Precision"]"\t"$h["METRIC.Recall"]"\t"$h["METRIC.F1_Score"]}' "$1"; }
