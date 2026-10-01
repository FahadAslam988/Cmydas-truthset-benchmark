#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# STEP 3 (TS4): re-score the EXISTING simulation calls (already filtered with the real-data rule and normalised,
# 05_benchmark/<caller>.filt.vcf.gz) against the known truth, restricted to callable.bed (72 contigs, 210.8 Mb),
# with the SAME allele-level vcfeval settings. The original genotype-aware result stays in 06_results/.
set -uo pipefail; source "$(dirname "$0")/common.sh"
D=$E/03_TS4_simulation_happy_allele; mkdir -p $D $E/tmp
declare -A SIMNAME=( [BCF]=BCF [Deepvariant]=StockDV [Freebayes]=Freebayes [Gatk]=Gatk [Varscan]=Varscan )
for cl in "${CALLERS[@]}"; do
  q=$D/${cl}_sim.altonly.vcf.gz
  [ -s $q ] || { $BCF view -i 'GT="alt"' $SIM/05_benchmark/${SIMNAME[$cl]}.filt.vcf.gz -Oz -o $q && $BCF index -f -t $q; }
  happy $SIM/01_mutated_genome/truth.vcf.gz $q $D/TS4_simulation_${cl} -T "$(c $SIM/01_mutated_genome/callable.bed)" &
done; wait
