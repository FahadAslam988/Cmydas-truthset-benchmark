#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# STEP 2 (TS1): score every caller against every consensus level (exactly 3, >=3, exactly 4, >=4, all 5), both samples: 50 hap.py runs.
set -uo pipefail; source "$(dirname "$0")/common.sh"
D=$E/02_TS1_consensus_happy; P=$E/00_query_altonly; mkdir -p $D $P $E/tmp
for g in "${GROUPS_[@]}"; do for cl in "${CALLERS[@]}"; do          # query: drop 0/0 and ./. (as in LOO)
  q=$P/${cl}_${g}.altonly.vcf.gz
  [ -s $q ] || { $BCF view -i 'GT="alt"' "$(query_snps $cl $g)" -Oz -o $q && $BCF index -f -t $q; }
done; done
for g in "${GROUPS_[@]}"; do for lv in n3 3plus n4 4plus 5tools; do for cl in "${CALLERS[@]}"; do
  happy $E/01_TS1_consensus_truth/TS1_consensus_${lv}_${g}.vcf.gz $P/${cl}_${g}.altonly.vcf.gz \
        $D/TS1_${lv}_${g}_${cl} &
  while [ "$(jobs -rp | wc -l)" -ge $MAXPAR ]; do wait -n; done
done; done; done; wait
