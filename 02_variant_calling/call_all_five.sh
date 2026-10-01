#!/usr/bin/env bash
# The five callers, with the options recorded in the header of each raw VCF (provenance in docs/PROVENANCE.md).
# Output names match those read by 03_filter_normalise/01_filter_all_callers.sh.
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
set -euo pipefail
mkdir -p "$CALLS_DIR"/{DeepVaraints,gatk_out,bcftools_out,freebayes_out,varscan_out}
for pair in "CH-NORMS1:normal" "Ab-NormS2:abnormal"; do s=${pair%%:*}; g=${pair##*:}; BAM="$BAM_DIR/$s.dedup.bam"

  # 1) DeepVariant 1.9.0, default human WGS model
  docker run --rm -v "$REF_DIR":/ref:ro -v "$BAM_DIR":/bam:ro -v "$CALLS_DIR/DeepVaraints":/output \
    google/deepvariant:1.9.0 /opt/deepvariant/bin/run_deepvariant --model_type=WGS \
    --ref=/ref/$(basename "$REF_FASTA") --reads=/bam/$s.dedup.bam \
    --output_vcf=/output/${s}_deepvariant_1_9.vcf.gz --output_gvcf=/output/${s}_deepvariant_1_9.g.vcf.gz \
    --num_shards=20 --make_examples_extra_args="small_model_call_multiallelics=false"

  # 2) GATK 4.6.2.0 HaplotypeCaller (GVCF, default settings) -> GenotypeGVCFs
  gatk HaplotypeCaller -R "$REF_FASTA" -I "$BAM" -O "$CALLS_DIR/gatk_out/$g.g.vcf.gz" --emit-ref-confidence GVCF
  gatk GenotypeGVCFs -R "$REF_FASTA" -V "$CALLS_DIR/gatk_out/$g.g.vcf.gz" -O "$CALLS_DIR/gatk_out/$g.gatk.genotyped.vcf.gz"

  # 3) BCFtools 1.22 mpileup | call -mv (defaults)
  bcftools mpileup --threads 16 -Ou -f "$REF_FASTA" "$BAM" | bcftools call --threads 16 -mv -Oz -o "$CALLS_DIR/bcftools_out/$g.vcf.gz"
  bcftools index "$CALLS_DIR/bcftools_out/$g.vcf.gz"

  # 4) FreeBayes 1.3.10 (options recorded in the VCF header ##commandline)
  freebayes -f "$REF_FASTA" --ploidy 2 --min-coverage 10 --min-base-quality 20 --min-mapping-quality 30 "$BAM" \
    | bgzip > "$CALLS_DIR/freebayes_out/$g.freebayes.vcf.gz"; tabix -p vcf "$CALLS_DIR/freebayes_out/$g.freebayes.vcf.gz"

  # 5) VarScan 2.4.6 mpileup2cns. VarScan writes no command line to the VCF header; options as in the author's
  #    released calling script (chelonia-snv-benchmarking, scripts/05_variant_calling/varscan.sh).
  samtools mpileup -f "$REF_FASTA" -q 30 -Q 20 -B "$BAM" \
    | varscan mpileup2cns --min-coverage 10 --min-avg-qual 20 --p-value 0.05 --variants --output-vcf 1 \
    | bgzip > "$CALLS_DIR/varscan_out/$g.varscan.vcf.gz"; tabix -p vcf "$CALLS_DIR/varscan_out/$g.varscan.vcf.gz"
done
