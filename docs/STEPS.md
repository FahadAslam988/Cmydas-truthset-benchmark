# Step-by-step workflow

Paths are variables from `config/paths.sh`. Samples: **CH-NORMS1** ("normal") and **Ab-NormS2** ("abnormal").

## Step 1: preprocessing (`01_preprocessing/`)

| Script | Input | Output | Settings |
|---|---|---|---|
| `01_fastp_trim.sh` | `$FASTQ_DIR/{W1-A,S1}_{1,2}.fastq.gz` | `$TRIM_DIR/*_trimmed.fastq.gz`, fastp html/json | fastp 0.23.4, `--detect_adapter_for_pe` |
| `02_bwa_align_markdup.sh` | trimmed reads, `$REF_FASTA` | `$BAM_DIR/<sample>.dedup.bam` (+ index, metrics) | BWA-MEM 0.7.17 default settings; read group LB lib1 / PL illumina / PU unit1; GATK 4.6.2.0 MarkDuplicates |

## Step 2: variant calling (`02_variant_calling/call_all_five.sh`)

| Caller | Version | Settings | Output |
|---|---|---|---|
| DeepVariant | 1.9.0 | `--model_type=WGS` (default model), `small_model_call_multiallelics=false` | `$CALLS_DIR/DeepVaraints/<sample>_deepvariant_1_9.vcf.gz` |
| GATK HaplotypeCaller | 4.6.2.0 | GVCF mode, defaults → GenotypeGVCFs | `$CALLS_DIR/gatk_out/<group>.gatk.genotyped.vcf.gz` |
| BCFtools | 1.22 | `mpileup -f REF` \| `call -mv` | `$CALLS_DIR/bcftools_out/<group>.vcf.gz` |
| FreeBayes | 1.3.10 | `--ploidy 2 --min-coverage 10 --min-base-quality 20 --min-mapping-quality 30` | `$CALLS_DIR/freebayes_out/<group>.freebayes.vcf.gz` |
| VarScan | 2.4.6 | `samtools mpileup -q 30 -Q 20 -B` \| `mpileup2cns --min-coverage 10 --min-avg-qual 20 --p-value 0.05 --variants` | `$CALLS_DIR/varscan_out/<group>.varscan.vcf.gz` |

## Step 3: filter, normalise, SNPs (`03_filter_normalise/`)

| Script | Does | Output |
|---|---|---|
| `01_filter_all_callers.sh` | keep FILTER PASS/"." AND alt-allele reads > 3 AND VAF > 0.02 (field names adapted per caller) | `$WORK_DIR/01_filtered/<Caller>_<group>.filtered.vcf.gz` |
| `02_extract_snps.sh` | SNPs from filtered VCFs (counts only; reporting) | `$WORK_DIR/02_snps/` |
| `03_normalise_and_extract_snps.sh` | `bcftools norm -f REF -m -any -c w`, then SNPs | `$WORK_DIR/03_normalized/`, **`$WORK_DIR/04_norm_snps/<Caller>_<group>_SNPs.vcf.gz`** (used for all scoring) |

## Step 4: truth sets (`04_truth_sets/`)

### TS1 consensus (`TS1_consensus/01_consensus_isec.sh`)
`bcftools isec -n +3`, `-n +4`, `-n =5` over the five `04_norm_snps` files → `$WORK_DIR/05_consensus/isec_<group>_<level>/sites.txt`.
The truth VCFs (≥3, ≥4, all 5, plus exactly 3 and exactly 4 by subtraction) are written in step 5.

### TS2 leave-one-out (`TS2_leave_one_out/`)
For each tested caller: `bcftools isec -n +k` over the **other four** callers, merged into one truth VCF
(`concat | sort | norm -d exact`); k = 2 (`01_build_loo_truth_k2.sh`), k = 3 and 4 (`02_build_and_score_loo_k3_k4.sh`).
Scored with hap.py (`03_score_loo_k2_happy.sh`, and within `02_...k3_k4.sh`) using the harmonised settings.

### TS4 simulation (`TS4_simulation/`)
| Script | Does |
|---|---|
| `s00_calibrate_from_consensus.sh` | density, Ti/Tv and het:hom from SNPs called by ≥3 of 5 callers in both samples |
| `s00a_measure_depth_insert_size.sh` | depth, insert size, read length from the real BAMs |
| `s01_select_contig_subset.sh` | 72 whole contigs (~211 Mb): all small scaffolds, one large and several mid-size chromosomes |
| `gen_truth_variants.py`, `s02_make_truth_and_haplotypes.sh` | known SNVs/indels (fixed seed) → `truth.vcf.gz`, two haplotype FASTAs (`bcftools consensus`), `callable.bed` |
| `s03a_build_art_error_profile.sh`, `s03b_simulate_reads_art.sh`, `art_chunk_worker.sh` | ART error profile learned from the real CH-NORMS1 reads; 2 × 151 bp reads, ~70× diploid |
| `s04_align_markdup.sh` | same alignment + duplicate marking as the real data |
| `s05a_call_deepvariant.sh`, `s05b_call_gatk_bcftools_freebayes_varscan.sh` | the five callers on the simulated BAM |
| `s06_filter_and_score_genotype_aware.sh` | the same filter as the real data, normalise, hap.py genotype-aware (supplementary) |

## Step 5: harmonised benchmark (`05_harmonised_benchmark/`)

One scoring method for every truth set: **hap.py 0.3.12, `--engine=vcfeval` (RTG 3.12.1), `--set-gt hom`
(allele-level), 0/0 and ./. removed first**; real data scored genome-wide, simulation within `callable.bed`.

| Script | Does |
|---|---|
| `s01_TS1_build_consensus_truth.sh` | consensus sites → truth VCFs (≥3, ≥4, all 5 from `sites.txt`; exactly 3 and exactly 4 by subtraction); every count checked |
| `s02_TS1_happy_consensus.sh` | 5 callers × 5 levels × 2 samples = 50 runs |
| `s03_TS4_happy_simulation_allele.sh` | 5 callers on the simulation |
| `s04_collect_results.sh` | `Paper1_all_truthsets_SNP_metrics.tsv` (TS1 + TS2 + TS4) and `Paper1_ranks_by_truthset.tsv` |
| `run_all.sh` | steps 1 to 4 |
| `s05_TS4_multicaller_consensus.sh` | multi-caller consensus (≥2, ≥3, ≥4, all 5) on the simulation vs single callers |
| `s07_indel_robustness.sh` | the same truth-set comparison for indels (supplementary) |
| `s06_bootstrap_ci.sh` | 95% CIs: paired block bootstrap (1,000 replicates, 1-Mb blocks, seed 20261001) of precision/recall/F1, pairwise F1 differences and rank stability |

## Results (`results/`)

Final tables of the paper and the source tables behind every number (raw, filtered and normalised counts, consensus
counts, leave-one-out, simulation genotype-aware check, simulation parameters).
