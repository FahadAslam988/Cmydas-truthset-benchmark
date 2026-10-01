#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# STEP 7 (robustness): repeat the truth-set comparison for INDELS (supplementary).
# Inputs: rerun/03_normalized/<Caller>_<group>.norm.vcf.gz (filtered + normalised, same as the SNP analysis) -> indels only.
# TS1 consensus (>=3, >=4, all 5) and TS2 leave-one-out (k = 2, 3, 4) built with bcftools isec (exact allele match,
# same as for SNPs); scored with the identical hap.py settings (vcfeval, --set-gt hom). Simulation indels: INDEL rows of
# the existing TS3 hap.py outputs (03_TS4_simulation_happy_allele).
set -uo pipefail; source "$(dirname "$0")/common.sh"
D=$E/07_indel_robustness; mkdir -p $D/calls $D/truth $D/happy $E/tmp
for g in "${GROUPS_[@]}"; do for c in "${CALLERS[@]}"; do
  o=$D/calls/${c}_${g}_INDELs.vcf.gz
  [ -s $o.tbi ] || { $BCF view -v indels -i 'GT="alt"' $R/03_normalized/${c}_${g}.norm.vcf.gz -Oz -o $o && $BCF index -f -t $o; }
done; done
mk_truth(){ # mk_truth <out.vcf.gz> <isec -n arg> <files...>
  local out=$1 n=$2; shift 2; [ -s $out.tbi ] && return 0
  local d=$out.isec; rm -rf $d; $BCF isec -n "$n" -p $d "$@" 2>>$LOG
  { echo '##fileformat=VCFv4.2'; awk '{print "##contig=<ID="$1",length="$2">"}' $REFDIR/Reference.fasta.fai
    echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'; printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\n'
    awk -v OFS='\t' '{print $1,$2,".",$3,$4,".","PASS",".","GT","1/1"}' $d/sites.txt; } | $BCF sort -Oz -o $out 2>>$LOG && $BCF index -f -t $out
  log "indel truth $(basename $out): $($BCF index -n $out) (sites $(wc -l < $d/sites.txt))"; rm -rf $d; }
for g in "${GROUPS_[@]}"; do
  all=(); for c in "${CALLERS[@]}"; do all+=($D/calls/${c}_${g}_INDELs.vcf.gz); done
  for p in "+3:3plus" "+4:4plus" "=5:5tools"; do mk_truth $D/truth/TS1_${p##*:}_${g}.vcf.gz "${p%%:*}" "${all[@]}"; done
  for t in "${CALLERS[@]}"; do oth=(); for c in "${CALLERS[@]}"; do [ $c != $t ] && oth+=($D/calls/${c}_${g}_INDELs.vcf.gz); done
    for k in 2 3 4; do n="+$k"; [ $k = 4 ] && n="=4"; mk_truth $D/truth/TS2_k${k}_${g}_wo_${t}.vcf.gz "$n" "${oth[@]}"; done; done
done
for g in "${GROUPS_[@]}"; do for c in "${CALLERS[@]}"; do
  for lv in 3plus 4plus 5tools; do happy $D/truth/TS1_${lv}_${g}.vcf.gz $D/calls/${c}_${g}_INDELs.vcf.gz $D/happy/TS1_${lv}_${g}_${c} &
    while [ "$(jobs -rp | wc -l)" -ge $MAXPAR ]; do wait -n; done; done
  for k in 2 3 4; do happy $D/truth/TS2_k${k}_${g}_wo_${c}.vcf.gz $D/calls/${c}_${g}_INDELs.vcf.gz $D/happy/TS2_k${k}_${g}_${c} &
    while [ "$(jobs -rp | wc -l)" -ge $MAXPAR ]; do wait -n; done; done
done; done; wait
O=$E/04_results/Paper1_INDEL_all_truthsets.tsv
irow(){ awk -F',' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $h["Type"]=="INDEL"&&$h["Filter"]=="ALL"{print $h["TRUTH.TP"]"\t"$h["QUERY.FP"]"\t"$h["TRUTH.FN"]"\t"$h["METRIC.Precision"]"\t"$h["METRIC.Recall"]"\t"$h["METRIC.F1_Score"]}' "$1"; }
printf "TruthSet\tLevel\tSample\tCaller\tTP\tFP\tFN\tPrecision\tRecall\tF1\n" > $O
for g in "${GROUPS_[@]}"; do for c in "${CALLERS[@]}"; do
  for lv in 3plus 4plus 5tools; do printf "TS1_consensus\t%s\t%s\t%s\t%s\n" $lv $g $c "$(irow $D/happy/TS1_${lv}_${g}_${c}.summary.csv)" >> $O; done
  for k in 2 3 4; do printf "TS2_leave_one_out\tk%s\t%s\t%s\t%s\n" $k $g $c "$(irow $D/happy/TS2_k${k}_${g}_${c}.summary.csv)" >> $O; done
done; done
for c in "${CALLERS[@]}"; do printf "TS4_simulation\tknown_truth\tsimulated\t%s\t%s\n" $c "$(irow $E/03_TS4_simulation_happy_allele/TS4_simulation_${c}.summary.csv)" >> $O; done
python3 - "$O" "$E/04_results/Paper1_INDEL_ranks_by_truthset.tsv" <<'PY'
import sys,csv,collections
g=collections.defaultdict(list)
for r in csv.DictReader(open(sys.argv[1]),delimiter='\t'): g[(r['TruthSet'],r['Level'],r['Sample'])].append(r)
with open(sys.argv[2],'w') as f:
    f.write('TruthSet\tLevel\tSample\tRanking(best->worst)\tSpread\n')
    for k,v in g.items():
        v=sorted(v,key=lambda r:-float(r['F1'] or 0))
        f.write('\t'.join(k)+'\t'+' > '.join(f"{r['Caller']}({float(r['F1']):.4f})" for r in v)+f"\t{float(v[0]['F1'])-float(v[-1]['F1']):.4f}\n")
PY
log "indel robustness done: $O"
