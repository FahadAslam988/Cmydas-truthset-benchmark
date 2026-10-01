#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# STEP 5 (multi-caller test): is a consensus of several callers more accurate than a single caller?
# Uses the SIMULATION (known truth) only. The five filtered + normalised simulation call sets (same rule as the real data)
# -> SNPs, alt-only -> bcftools isec: SNPs called by >=2, >=3, >=4, all 5 callers -> one sites-only VCF each (GT 1/1)
# -> hap.py with the identical harmonised settings (vcfeval, --set-gt hom, callable.bed) -> compared with single callers.
set -uo pipefail; source "$(dirname "$0")/common.sh"
D=$E/05_TS4_multicaller_consensus; mkdir -p $D $E/tmp
declare -A SIMNAME=( [BCF]=BCF [Deepvariant]=StockDV [Freebayes]=Freebayes [Gatk]=Gatk [Varscan]=Varscan )
files=()
for cl in "${CALLERS[@]}"; do
  q=$D/${cl}_sim.SNP.altonly.vcf.gz
  [ -s $q.tbi ] || { $BCF view -v snps -i 'GT="alt"' $SIM/05_benchmark/${SIMNAME[$cl]}.filt.vcf.gz -Oz -o $q && $BCF index -f -t $q; }
  files+=($q)
done
for pair in "+2:ge2" "+3:ge3" "+4:ge4" "=5:all5"; do
  n=${pair%%:*}; lv=${pair##*:}; out=$D/multicaller_${lv}.vcf.gz
  if [ ! -s $out.tbi ]; then
    rm -rf $D/isec_$lv; $BCF isec -n "$n" -p $D/isec_$lv "${files[@]}" 2>>$LOG
    { echo '##fileformat=VCFv4.2'
      awk '{print "##contig=<ID="$1",length="$2">"}' $REFDIR/Reference.fasta.fai
      echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
      printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCONSENSUS\n'
      awk -v OFS='\t' '{print $1,$2,".",$3,$4,".","PASS",".","GT","1/1"}' $D/isec_$lv/sites.txt
    } | $BCF sort -Oz -o $out 2>>$LOG && $BCF index -f -t $out
    rm -f $D/isec_$lv/0*.vcf
  fi
  log "multi-caller $lv: records=$($BCF index -n $out) sites.txt=$(wc -l < $D/isec_$lv/sites.txt)"
  happy $SIM/01_mutated_genome/truth.vcf.gz $out $D/TS4_multicaller_${lv} -T "$(c $SIM/01_mutated_genome/callable.bed)" &
done; wait
O=$E/04_results/Paper1_multicaller_vs_single_simulation.tsv
printf "CallSet\tType\tTP\tFP\tFN\tPrecision\tRecall\tF1\n" > $O
for cl in "${CALLERS[@]}"; do printf "%s\tsingle_caller\t%s\n" $cl "$(snp_row $E/03_TS4_simulation_happy_allele/TS4_simulation_${cl}.summary.csv)" >> $O; done
for lv in ge2 ge3 ge4 all5; do printf "consensus_%s\tmulti_caller\t%s\n" $lv "$(snp_row $D/TS4_multicaller_${lv}.summary.csv)" >> $O; done
column -t $O | tee -a $LOG
