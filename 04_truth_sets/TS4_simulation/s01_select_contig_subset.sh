#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# ============================================================================
# STAGE 1 (0.5) — Build a REPRESENTATIVE ~200 Mb subset of the reference.
# Rationale: benchmarking on a representative genome fraction yields statistically
# valid caller metrics at a fraction of the compute (standard GIAB-style practice;
# Olson 2023, Krusche 2019). We span the assembly's size range so the subset is
# NOT artificially easy:
#   (A) ALL small unplaced scaffolds (<1 Mb)  -> repeat-rich, hard-to-call sequence
#   (B) one large chromosome (smallest >=100 Mb) -> chromosome-scale euchromatin
#   (C) mid-size chromosomes (10-100 Mb), smallest-first, until total >= TARGET
# Whole contigs only => native coordinates => truth VCF stays exact (no lift-over).
# Resumable + logged.
# ============================================================================
set -euo pipefail

BASE="${SIM_DIR}"
REF="${REF_FASTA}"
SAMTOOLS="samtools"
FAI="$REF.fai"
OUT="$BASE/outputs/01_mutated_genome"
SUBFA="$OUT/subset.fa"
SUBBED="$OUT/subset.bed"
MAN="$OUT/subset_contigs.txt"
LOG="$BASE/run.log"
TARGET=200000000          # ~200 Mb

mkdir -p "$OUT"
if [ -s "$SUBFA" ] && [ -s "$SUBFA.fai" ]; then
  echo "[s01] subset already built -> $SUBFA (delete to re-run)"; exit 0
fi
log(){ echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }
log "===== STAGE 1 subset START (target ~$((TARGET/1000000)) Mb) ====="

# --- select contigs deterministically --------------------------------------
#  columns of .fai: name length ...
awk -v target="$TARGET" '
{
  name[NR]=$1; len[NR]=$2; n=NR
}
END{
  sum=0
  # (A) all scaffolds < 1 Mb
  for(i=1;i<=n;i++) if(len[i]<1000000){ pick[i]=1; sum+=len[i] }
  # (B) smallest chromosome >= 100 Mb
  big=0; bigidx=0
  for(i=1;i<=n;i++) if(len[i]>=100000000){ if(big==0||len[i]<big){big=len[i]; bigidx=i} }
  if(bigidx){ pick[bigidx]=1; sum+=len[bigidx] }
  # (C) mid chromosomes 10-100 Mb, smallest-first until sum>=target
  #   simple selection sort over mid class
  for(pass=1; sum<target; pass++){
    msmall=0; midx=0
    for(i=1;i<=n;i++){
      if(pick[i]) continue
      if(len[i]>=10000000 && len[i]<100000000){
        if(msmall==0||len[i]<msmall){ msmall=len[i]; midx=i }
      }
    }
    if(midx==0) break
    pick[midx]=1; sum+=len[midx]
  }
  # emit manifest (name \t length \t class)
  for(i=1;i<=n;i++) if(pick[i]){
    cls=(len[i]<1000000)?"scaffold_<1Mb":((len[i]>=100000000)?"chrom_>=100Mb":"chrom_10-100Mb")
    printf "%s\t%d\t%s\n", name[i], len[i], cls
  }
  printf "#SUBSET_TOTAL_BP\t%d\n", sum > "/dev/stderr"
}' "$FAI" > "$MAN" 2> "$OUT/.subtotal"

TOT=$(awk '{print $2}' "$OUT/.subtotal"); rm -f "$OUT/.subtotal"
NC=$(wc -l < "$MAN")
PCT=$(awk -v t="$TOT" 'BEGIN{printf "%.1f", 100*t/2134000000}')
log "selected $NC contigs, total $(awk -v t=$TOT 'BEGIN{printf "%.1f", t/1e6}') Mb (~${PCT}% of assembly)"
log "class breakdown:"; cut -f3 "$MAN" | sort | uniq -c | tee -a "$LOG"

# --- extract subset FASTA + index + BED ------------------------------------
log "extracting subset FASTA (samtools faidx)..."
cut -f1 "$MAN" | xargs "$SAMTOOLS" faidx "$REF" > "$SUBFA" 2>>"$LOG"
"$SAMTOOLS" faidx "$SUBFA"
awk 'BEGIN{OFS="\t"} {print $1,0,$2}' "$SUBFA.fai" > "$SUBBED"
log "subset.fa + subset.fa.fai + subset.bed written"

# --- record to params.txt ---------------------------------------------------
{
  echo "# ---- STAGE 1 subset ----"
  echo "SUBSET_FA=$SUBFA"
  echo "SUBSET_BED=$SUBBED"
  echo "SUBSET_BP=$TOT"
  echo "SUBSET_CONTIGS=$NC"
  echo "SUBSET_PCT_OF_GENOME=$PCT"
  echo "# STAGE1_DONE $(date '+%F %T')"
} >> "$BASE/params.txt"
log "===== STAGE 1 DONE -> $SUBFA ====="
echo ">>> Review $MAN ; next: s02 simuG (variant counts from measured SNV density)"
