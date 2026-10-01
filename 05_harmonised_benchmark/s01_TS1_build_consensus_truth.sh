#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# STEP 1 (TS1): turn the existing consensus sites (rerun/05_consensus/isec_<group>_<level>/sites.txt, made by
# stage4_consensus.sh = bcftools isec on the 04_norm_snps files) into one truth VCF per group x level.
# Allele-level scoring uses only CHROM/POS/REF/ALT, so each consensus site is written once with GT=1/1.
# Check: records in the truth VCF must equal the lines in sites.txt (= consensus_counts.tsv).
set -uo pipefail; source "$(dirname "$0")/common.sh"
D=$E/01_TS1_consensus_truth; mkdir -p $D
for g in "${GROUPS_[@]}"; do for lv in 3plus 4plus 5tools; do
  s=$R/05_consensus/isec_${g}_${lv}/sites.txt; out=$D/TS1_consensus_${lv}_${g}.vcf.gz
  [ -s $out.tbi ] && { log "TS1 truth exists: $(basename $out)"; continue; }   # index written last = file complete
  { echo '##fileformat=VCFv4.2'
    awk '{print "##contig=<ID="$1",length="$2">"}' $REFDIR/Reference.fasta.fai
    echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\n'
    awk -v OFS='\t' '{print $1,$2,".",$3,$4,".","PASS",".","GT","1/1"}' $s
  } | $BCF sort -Oz -o $out 2>>$LOG && $BCF index -f -t $out
  n=$($BCF index -n $out); e=$(wc -l < $s)
  log "TS1 truth $g $lv: records=$n sites.txt=$e $([ "$n" = "$e" ] && echo MATCH || echo MISMATCH)"
done; done

# Exact strata (Fahad's original 5-level design): exactly 3 = (>=3) minus (>=4); exactly 4 = (>=4) minus (all 5).
# bcftools isec -C keeps records of the first file absent from the second (exact CHROM/POS/REF/ALT match).
# Check: record count must equal the difference of the two sites.txt counts.
for g in "${GROUPS_[@]}"; do for pair in "n3:3plus:4plus" "n4:4plus:5tools"; do
  IFS=: read lv a b <<< "$pair"; out=$D/TS1_consensus_${lv}_${g}.vcf.gz
  [ -s $out.tbi ] && { log "TS1 truth exists: $(basename $out)"; continue; }
  $BCF isec -C -w1 $D/TS1_consensus_${a}_${g}.vcf.gz $D/TS1_consensus_${b}_${g}.vcf.gz -Oz -o $out 2>>$LOG && $BCF index -f -t $out
  n=$($BCF index -n $out)
  e=$(( $(wc -l < $R/05_consensus/isec_${g}_${a}/sites.txt) - $(wc -l < $R/05_consensus/isec_${g}_${b}/sites.txt) ))
  log "TS1 truth $g $lv (exactly): records=$n expected=$e $([ "$n" = "$e" ] && echo MATCH || echo MISMATCH)"
done; done
