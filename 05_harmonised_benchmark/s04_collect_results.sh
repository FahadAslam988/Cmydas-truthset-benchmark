#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# STEP 4: one table with every caller x truth set x sample (SNPs), plus ranks and spread per truth set.
# TS2 (leave-one-out) is taken from 07_loo_sensitivity/summary/sensitivity_ranking.tsv: it was produced with the
# identical call sets and hap.py settings, so it is already harmonised (not re-run).
set -uo pipefail; source "$(dirname "$0")/common.sh"
O=$E/04_results; mkdir -p $O; T=$O/Paper1_all_truthsets_SNP_metrics.tsv
printf "TruthSet\tLevel\tSample\tCaller\tTP\tFP\tFN\tPrecision\tRecall\tF1\n" > $T
for g in "${GROUPS_[@]}"; do for lv in n3 3plus n4 4plus 5tools; do for cl in "${CALLERS[@]}"; do
  f=$E/02_TS1_consensus_happy/TS1_${lv}_${g}_${cl}.summary.csv
  [ -s $f ] && printf "TS1_consensus\t%s\t%s\t%s\t%s\n" $lv $g $cl "$(snp_row $f)" >> $T
done; done; done
awk -F'\t' -v OFS='\t' 'NR>1{print "TS2_leave_one_out","k"$1,$2,$3,$4,$5,$6,$7,$8,$9}' \
  $R/07_loo_sensitivity/summary/sensitivity_ranking.tsv >> $T
for cl in "${CALLERS[@]}"; do f=$E/03_TS4_simulation_happy_allele/TS4_simulation_${cl}.summary.csv
  [ -s $f ] && printf "TS4_simulation\tknown_truth\tsimulated\t%s\t%s\n" $cl "$(snp_row $f)" >> $T; done
# ranks (1 = highest F1) and spread within each truth set / level / sample
python3 - "$T" "$O" <<'PY'
import sys,csv,collections
rows=list(csv.DictReader(open(sys.argv[1]),delimiter='\t'))
grp=collections.defaultdict(list)
for r in rows: grp[(r['TruthSet'],r['Level'],r['Sample'])].append(r)
with open(sys.argv[2]+'/Paper1_ranks_by_truthset.tsv','w') as f:
    f.write('TruthSet\tLevel\tSample\tRanking(best->worst)\tBestF1\tWorstF1\tSpread\n')
    for k,v in grp.items():
        v=sorted(v,key=lambda r:-float(r['F1']))
        f.write('\t'.join(k)+'\t'+' > '.join(f"{r['Caller']}({float(r['F1']):.5f})" for r in v)
                +f"\t{float(v[0]['F1']):.5f}\t{float(v[-1]['F1']):.5f}\t{float(v[0]['F1'])-float(v[-1]['F1']):.5f}\n")
PY
column -t -s$'\t' $O/Paper1_ranks_by_truthset.tsv | tee -a $LOG
