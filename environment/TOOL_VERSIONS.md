# Tool versions

| Tool | Version | Source |
|---|---|---|
| fastp | 0.23.4 | conda |
| BWA-MEM | 0.7.17-r1188 | BAM @PG |
| SAMtools | 1.10 (alignment), 1.22.1 (analysis) | BAM @PG, `bcftools-env.yml` |
| GATK | 4.6.2.0 | `gatk-env.yml`, VCF header |
| DeepVariant | 1.9.0 | Docker `google/deepvariant:1.9.0` |
| BCFtools / HTSlib | 1.22 / 1.22.1 | `bcftools-env.yml` |
| FreeBayes | 1.3.10 | `freebayes-env.yml`, VCF header |
| VarScan | 2.4.6 | `varscan-env.yml` |
| ART | 2016.06.05 (art_illumina, art_profiler_illumina) | `art-env.yml` |
| hap.py | 0.3.12 (RTG Tools 3.12.1 vcfeval engine) | Docker `jmcdani20/hap.py:v0.3.12` |
| Python | 3.10+ | |
