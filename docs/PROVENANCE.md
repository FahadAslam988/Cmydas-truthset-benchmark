# Provenance

Every script in `03_` to `05_` is a copy of the script that produced the results, with only these mechanical
changes (made by `build_repo.py`, kept outside the repository): personal absolute paths replaced by variables from
`config/paths.sh`; absolute tool paths replaced by tool names; a line sourcing `config/paths.sh` added after the
shebang; lines that ran a retrained DeepVariant model removed (that model belongs to a separate study).

## Script origin

| Repository file | Original script | Paths replaced | Retrained-model lines removed |
|---|---|---|---|
| `03_filter_normalise/01_filter_all_callers.sh` | `rerun/stage1_filter.sh` | 3 | 0 |
| `03_filter_normalise/02_extract_snps.sh` | `rerun/stage2_snps.sh` | 2 | 0 |
| `03_filter_normalise/03_normalise_and_extract_snps.sh` | `rerun/stage3_normalize_snps.sh` | 4 | 0 |
| `04_truth_sets/TS1_consensus/01_consensus_isec.sh` | `rerun/stage4_consensus.sh` | 2 | 0 |
| `04_truth_sets/TS2_leave_one_out/01_build_loo_truth_k2.sh` | `rerun/stage5_leaveoneout.sh` | 3 | 0 |
| `04_truth_sets/TS2_leave_one_out/02_build_and_score_loo_k3_k4.sh` | `rerun/07_loo_sensitivity/sensitivity_loo.sh` | 2 | 0 |
| `04_truth_sets/TS2_leave_one_out/03_score_loo_k2_happy.sh` | `rerun/stage5_leaveoneout_happy.sh` | 2 | 0 |
| `04_truth_sets/TS4_simulation/s00_calibrate_from_consensus.sh` | `rerun/09_simulation_truthset/RUN_DeepVariant_fair_benchmark_simuG-ART/scripts/s00c_recalibrate_consensus.sh` | 2 | 0 |
| `04_truth_sets/TS4_simulation/s00a_measure_depth_insert_size.sh` | `rerun/09_simulation_truthset/RUN_DeepVariant_fair_benchmark_simuG-ART/scripts/s00_calibrate.sh` | 5 | 0 |
| `04_truth_sets/TS4_simulation/s01_select_contig_subset.sh` | `rerun/09_simulation_truthset/RUN_DeepVariant_fair_benchmark_simuG-ART/scripts/s01_subset.sh` | 3 | 0 |
| `04_truth_sets/TS4_simulation/gen_truth_variants.py` | `rerun/09_simulation_truthset/RUN_DeepVariant_fair_benchmark_simuG-ART/scripts/gen_truth_variants.py` | 0 | 0 |
| `04_truth_sets/TS4_simulation/s02_make_truth_and_haplotypes.sh` | `rerun/09_simulation_truthset/RUN_DeepVariant_fair_benchmark_simuG-ART/scripts/s02b_make_truth_bcftools.sh` | 5 | 0 |
| `04_truth_sets/TS4_simulation/s03a_build_art_error_profile.sh` | `rerun/09_simulation_truthset/RUN_DeepVariant_fair_benchmark_simuG-ART/scripts/s03a_build_error_profile.sh` | 4 | 0 |
| `04_truth_sets/TS4_simulation/s03b_simulate_reads_art.sh` | `rerun/09_simulation_truthset/RUN_DeepVariant_fair_benchmark_simuG-ART/scripts/s03b_art_reads.sh` | 2 | 0 |
| `04_truth_sets/TS4_simulation/art_chunk_worker.sh` | `rerun/09_simulation_truthset/RUN_DeepVariant_fair_benchmark_simuG-ART/scripts/art_chunk_worker.sh` | 2 | 0 |
| `04_truth_sets/TS4_simulation/s04_align_markdup.sh` | `rerun/09_simulation_truthset/RUN_DeepVariant_fair_benchmark_simuG-ART/scripts/s05_align.sh` | 5 | 0 |
| `04_truth_sets/TS4_simulation/s05a_call_deepvariant.sh` | `rerun/09_simulation_truthset/RUN_DeepVariant_fair_benchmark_simuG-ART/scripts/s06a_deepvariant.sh` | 4 | 6 |
| `04_truth_sets/TS4_simulation/s05b_call_gatk_bcftools_freebayes_varscan.sh` | `rerun/09_simulation_truthset/RUN_DeepVariant_fair_benchmark_simuG-ART/scripts/s06b_callers.sh` | 6 | 0 |
| `04_truth_sets/TS4_simulation/s06_filter_and_score_genotype_aware.sh` | `rerun/09_simulation_truthset/RUN_DeepVariant_fair_benchmark_simuG-ART/scripts/s07_benchmark.sh` | 3 | 0 |
| `05_harmonised_benchmark/common.sh` | `rerun/10_Paper1_harmonised_benchmark_consensus_LOO_simulation/scripts/common.sh` | 2 | 0 |
| `05_harmonised_benchmark/s01_TS1_build_consensus_truth.sh` | `rerun/10_Paper1_harmonised_benchmark_consensus_LOO_simulation/scripts/s01_TS1_build_consensus_truth.sh` | 0 | 0 |
| `05_harmonised_benchmark/s02_TS1_happy_consensus.sh` | `rerun/10_Paper1_harmonised_benchmark_consensus_LOO_simulation/scripts/s02_TS1_happy_consensus.sh` | 0 | 0 |
| `05_harmonised_benchmark/s03_TS4_happy_simulation_allele.sh` | `rerun/10_Paper1_harmonised_benchmark_consensus_LOO_simulation/scripts/s03_TS4_happy_simulation_allele.sh` | 0 | 0 |
| `05_harmonised_benchmark/s04_collect_results.sh` | `rerun/10_Paper1_harmonised_benchmark_consensus_LOO_simulation/scripts/s04_collect_results.sh` | 0 | 0 |
| `05_harmonised_benchmark/run_all.sh` | `rerun/10_Paper1_harmonised_benchmark_consensus_LOO_simulation/scripts/run_all.sh` | 0 | 0 |
| `05_harmonised_benchmark/s07_indel_robustness.sh` | `rerun/10_Paper1_harmonised_benchmark_consensus_LOO_simulation/scripts/s07_indel_robustness.sh` | 0 | 0 |
| `05_harmonised_benchmark/s06_bootstrap_ci.sh` | `rerun/10_Paper1_harmonised_benchmark_consensus_LOO_simulation/scripts/s06_bootstrap_ci.sh` | 0 | 0 |
| `05_harmonised_benchmark/s05_TS4_multicaller_consensus.sh` | `rerun/10_Paper1_harmonised_benchmark_consensus_LOO_simulation/scripts/s05_TS4_multicaller_consensus.sh` | 0 | 0 |

`01_preprocessing/` and `02_variant_calling/` were written from the commands recorded in the files themselves:

| Step | Evidence |
|---|---|
| BWA-MEM | `@PG ID:bwa VN:0.7.17-r1188 CL:bwa mem -t 15 <ref> <R1> <R2>` in the dedup BAM header (no `-M`) |
| samtools | `@PG samtools 1.10 view -bS -`, `sort` |
| Read group | `@RG ID:CH-NORMS1 LB:lib1 PL:illumina SM:CH-NORMS1 PU:unit1` |
| MarkDuplicates | `@PG ID:MarkDuplicates VN:4.6.2.0` |
| GATK (the earlier repository script shows only HaplotypeCaller; the VCF header proves GenotypeGVCFs was also run) | `##GATKCommandLine` HaplotypeCaller `--emit-ref-confidence GVCF` (4.6.2.0, defaults) and GenotypeGVCFs (4.6.2.0, defaults) |
| BCFtools | `##bcftoolsCommand=mpileup --threads 16 -Ou -f <ref> <sample>.dedup.bam`; `call --threads 16 -mv` |
| FreeBayes (the earlier repository script shows defaults; the VCF header, which records the actual run, is followed) | `##commandline="freebayes -f <ref> --ploidy 2 --min-coverage 10 --min-base-quality 20 --min-mapping-quality 30 <sample>.dedup.bam"` |
| DeepVariant | `##DeepVariant_version=1.9.0`; run command from the project log (`--model_type=WGS`, `small_model_call_multiallelics=false`) |
| VarScan | VarScan writes no command line; options from the author's calling script `chelonia-snv-benchmarking/scripts/05_variant_calling/varscan.sh` (identical to the project command log) |
| fastp | project command log (`--detect_adapter_for_pe`) |
