# ============================================================================
# Edit these paths for your system, then:  export CMYDAS_CONFIG=/full/path/to/config/paths.sh
# Every script in this repository sources this file.
# ============================================================================
export DATA_ROOT=/path/to/data                 # parent folder; must contain WORK_DIR and REF_DIR (mounted into Docker)
export WORK_DIR=$DATA_ROOT/work                # all intermediate and result folders are created here
export REF_DIR=$DATA_ROOT/ref                  # rCheMyd1.pri.v2 (RefSeq GCF_015237465.2)
export REF_FASTA=$REF_DIR/Reference.fasta      # + .fai, .dict, BWA index
export REF_SDF=$REF_DIR/Reference.sdf          # rtg format -i Reference.fasta -o Reference.sdf
export FASTQ_DIR=$DATA_ROOT/fastq              # raw reads from NCBI SRA (see README)
export TRIM_DIR=$DATA_ROOT/trimmed
export BAM_DIR=$DATA_ROOT/bam                  # <sample>.dedup.bam
export CALLS_DIR=$DATA_ROOT/calls              # raw caller VCFs (see 02_variant_calling)
export SIM_DIR=$WORK_DIR/simulation            # TS4 simulation working folder
export BENCH_DIR=$WORK_DIR/harmonised_benchmark
export CONDA_SH=$HOME/miniconda3/etc/profile.d/conda.sh
export THREADS=16
