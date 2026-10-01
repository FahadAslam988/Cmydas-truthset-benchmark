# Cmydas-truthset-benchmark

Code for the study **"Truth-set construction changes variant-caller rankings in a non-model vertebrate, the green
sea turtle (*Chelonia mydas*)"** (manuscript in preparation for *Scientific Reports*).

Five single-nucleotide variant (SNV) callers (DeepVariant, GATK HaplotypeCaller, BCFtools, FreeBayes, VarScan)
were run on whole-genome sequencing data from two green sea turtles and scored against three kinds of truth set
that differ in how independent they are of the callers themselves:

| Truth set | How it is built | Independent of the callers? |
|---|---|---|
| **TS1 consensus** | SNPs called by exactly 3, ≥3, exactly 4, ≥4 or all 5 callers (every site, whichever callers reported it) | no: the tested caller votes in its own truth |
| **TS2 leave-one-out** | for each tested caller, SNPs called by ≥k (k = 2, 3, 4) of the **other four** callers | partly: no self-vote, but still caller-derived |
| **TS4 simulation** | variants planted in the reference at a density and spectrum calibrated from the real data; reads simulated with an error profile learned from the real reads | yes: the answer is known by construction |

All truth sets are scored with the **same** method (hap.py 0.3.12, vcfeval engine, allele-level matching), so
differences in caller ranking between truth sets come from the truth sets alone.

## Main result

SNP F1, CH-NORMS1 (Ab-NormS2 within 0.001). Full tables: `results/`.

| Truth set | Ranking (best → worst) | Spread |
|---|---|---|
| Consensus ≥3 | BCFtools > VarScan > GATK > DeepVariant > FreeBayes | 0.059 |
| Consensus ≥4 | VarScan > DeepVariant > BCFtools > GATK > FreeBayes | 0.048 |
| Consensus, all 5 | FreeBayes > DeepVariant > VarScan > GATK > BCFtools | 0.069 |
| Leave-one-out k=2 | BCFtools > VarScan > GATK > DeepVariant > FreeBayes | 0.064 |
| Leave-one-out k=3 | VarScan > BCFtools > DeepVariant > GATK > FreeBayes | 0.052 |
| Leave-one-out k=4 | FreeBayes > DeepVariant > VarScan > GATK > BCFtools | 0.017 |
| Simulation (known truth) | FreeBayes > BCFtools > GATK > VarScan > DeepVariant | 0.001 |

---

## Data

| Data | Source |
|---|---|
| Raw reads (2 individuals, Illumina NovaSeq 6000, 2 × 151 bp) | NCBI SRA BioProject **PRJNA______** (to be added) |
| Reference genome rCheMyd1.pri.v2 | NCBI RefSeq **GCF_015237465.2** |
| Call sets, truth sets, hap.py outputs | Zenodo **10.5281/zenodo.______** (to be added) |

## Requirements

Linux, Docker, conda. Exact versions: `environment/TOOL_VERSIONS.md`; conda environments: `environment/*.yml`.

```
conda env create -f environment/bcftools-env.yml   # likewise gatk-env, freebayes-env, varscan-env, art-env
docker pull google/deepvariant:1.9.0
docker pull jmcdani20/hap.py:v0.3.12
```

## How to run

1. Copy `config/paths.sh`, set the paths, then `export CMYDAS_CONFIG=/full/path/to/config/paths.sh`.
2. Run the steps in order (inputs, outputs and settings of every step: `docs/STEPS.md`):

| Step | Folder | What it does |
|---|---|---|
| 1 | `01_preprocessing/` | fastp trimming; BWA-MEM alignment; read groups; MarkDuplicates |
| 2 | `02_variant_calling/` | the five callers |
| 3 | `03_filter_normalise/` | common filter (PASS, alt reads > 3, VAF > 0.02); `bcftools norm -m -any`; SNPs only |
| 4 | `04_truth_sets/` | TS1 consensus, TS2 leave-one-out, TS4 simulation (build + calls on simulated reads) |
| 5 | `05_harmonised_benchmark/` | scores every caller against every truth set with identical hap.py settings; collects one table |
| 6 | `06_figures/` | figures of the paper |

## Notes

- The analysis scripts are the ones that produced the published results; only absolute paths were replaced by the
  variables in `config/paths.sh` (`docs/PROVENANCE.md` gives the original file of every script and the evidence for
  every command).
- `04_truth_sets/TS2_leave_one_out/01_build_loo_truth_k2.sh` also contains an earlier RTG vcfeval scoring step;
  the reported leave-one-out results come from `03_score_loo_k2_happy.sh` and `02_build_and_score_loo_k3_k4.sh` (hap.py).
- `04_truth_sets/TS4_simulation/s06_filter_and_score_genotype_aware.sh` gives the genotype-aware simulation result
  (supplementary); the main allele-level result is produced in step 5.

## Citation

Aslam, F. et al. Truth-set construction changes variant-caller rankings in a non-model vertebrate, the green sea
turtle (*Chelonia mydas*). (in preparation). Code archive: Zenodo DOI to be added.

## License

MIT (see `LICENSE`).
