#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# ============================================================================
# STAGE 2b (replaces slow simuG) — build the TRUTH SET with bcftools consensus.
# Method precedent: Hall et al. 2024, eLife (bcftools consensus applies a truthset
# of variants to the reference to make a mutated genome). Fast (minutes), exact,
# native coordinates, reproducible (fixed seed). Produces:
#   truth.vcf.gz (+ .tbi)         genotyped answer key (het 0/1, hom 1/1)
#   hap1.simseq.genome.fa         reference + ALL variants   (for ART)
#   hap2.simseq.genome.fa         reference + HOM variants   (for ART)
#   callable.bed                  non-N regions for scoring
# Resumable + logged.
# ============================================================================
set -uo pipefail

BASE="${SIM_DIR}"
PARAMS="$BASE/params.txt"
OUT="$BASE/outputs/01_mutated_genome"
LOG="$BASE/run.log"
SEED=20260713
BCFTOOLS=bcftools
BGZIP=bgzip
TABIX=tabix
SAMTOOLS=samtools

getp(){ grep -E "^$1=" "$PARAMS" | head -1 | cut -d= -f2- | awk '{print $1}'; }
SUBFA=$(getp SUBSET_FA); N_SNP=$(getp N_SNP); N_INDEL=$(getp N_INDEL)
TITV=$(getp SIM_TITV);   HETF=$(getp HET_FRACTION)
TRUTH="$OUT/truth.vcf.gz"; CALLABLE="$OUT/callable.bed"
log(){ echo "[$(date '+%F %T')] [s02b] $*" | tee -a "$LOG"; }

if [ -s "$TRUTH" ] && [ -s "$OUT/hap1.simseq.genome.fa" ] && [ -s "$OUT/hap2.simseq.genome.fa" ]; then
  log "already done -> $TRUTH (delete to re-run)"; exit 0
fi
log "===== STAGE 2b (bcftools truth set) START ====="
log "params: N_SNP=$N_SNP N_INDEL=$N_INDEL TiTv=$TITV HET_FRAC=$HETF seed=$SEED"

# --- 1. generate truth + hap1 + hap2 raw VCFs -------------------------------
log "generating variants (python)..."
python3 "$BASE/scripts/gen_truth_variants.py" "$SUBFA" "$N_SNP" "$N_INDEL" "$TITV" "$HETF" "$SEED" "$OUT" 2>&1 | tee -a "$LOG"
[ -s "$OUT/truth.raw.vcf" ] || { log "ERROR: generator produced no truth.raw.vcf"; exit 2; }

# --- 2. sort/normalize/compress/index each VCF ------------------------------
prep(){ # in.raw.vcf -> out.vcf.gz (+tbi), normalized against SUBFA
  local raw="$1" outgz="$2"
  "$BCFTOOLS" sort "$raw" -Oz -o "$outgz.tmp.gz" 2>>"$LOG"
  "$BCFTOOLS" norm -f "$SUBFA" -m -any "$outgz.tmp.gz" -Oz -o "$outgz" 2>>"$LOG"
  "$TABIX" -f -p vcf "$outgz" 2>>"$LOG"; rm -f "$outgz.tmp.gz"
}
log "sorting/normalizing/indexing VCFs..."
prep "$OUT/truth.raw.vcf" "$TRUTH"
prep "$OUT/hap1.raw.vcf"  "$OUT/hap1.vcf.gz"
prep "$OUT/hap2.raw.vcf"  "$OUT/hap2.vcf.gz"

# --- 3. build hap1 + hap2 genomes with bcftools consensus -------------------
log "building hap1 genome (bcftools consensus, ALL variants)..."
"$BCFTOOLS" consensus -f "$SUBFA" "$OUT/hap1.vcf.gz" > "$OUT/hap1.simseq.genome.fa" 2>>"$LOG"
log "building hap2 genome (bcftools consensus, HOM variants)..."
"$BCFTOOLS" consensus -f "$SUBFA" "$OUT/hap2.vcf.gz" > "$OUT/hap2.simseq.genome.fa" 2>>"$LOG"
"$SAMTOOLS" faidx "$OUT/hap1.simseq.genome.fa" 2>>"$LOG"
"$SAMTOOLS" faidx "$OUT/hap2.simseq.genome.fa" 2>>"$LOG"

# --- 4. callable BED = non-N regions of the subset --------------------------
log "building callable.bed (non-N regions)..."
python3 - "$SUBFA" "$CALLABLE" <<'PY' 2>>"$LOG"
import sys, re
fa, out = sys.argv[1], sys.argv[2]
with open(fa) as f, open(out, "w") as o:
    name=None; seq=[]
    def flush(n,s):
        for m in re.finditer(r'[ACGTacgt]+', s): o.write(f"{n}\t{m.start()}\t{m.end()}\n")
    for line in f:
        if line[0]=='>':
            if name is not None: flush(name,''.join(seq))
            name=line[1:].split()[0]; seq=[]
        else: seq.append(line.strip())
    if name is not None: flush(name,''.join(seq))
PY

# --- 5. verify --------------------------------------------------------------
log "verifying truth set..."
NT=$("$BCFTOOLS" view -H "$TRUTH" | wc -l)
NSNP=$("$BCFTOOLS" view -H -v snps "$TRUTH" | wc -l)
NIND=$("$BCFTOOLS" view -H -v indels "$TRUTH" | wc -l)
TITV_OUT=$("$BCFTOOLS" stats "$TRUTH" 2>/dev/null | awk -F'\t' '/^TSTV/{print $5}')
read HET HOM < <("$BCFTOOLS" stats -s - "$TRUTH" 2>/dev/null | awk -F'\t' '/^PSC/{print $6, $5}')
CBP=$(awk '{s+=$3-$2}END{print s}' "$CALLABLE")
log "TRUTH: total=$NT SNP=$NSNP INDEL=$NIND | Ti/Tv=$TITV_OUT | het=$HET hom=$HOM | callable=$(awk -v b=$CBP 'BEGIN{printf "%.1f",b/1e6}')Mb"
log "hap1 genome: $(grep -c '^>' $OUT/hap1.simseq.genome.fa) contigs; hap2: $(grep -c '^>' $OUT/hap2.simseq.genome.fa)"

# --- record + cleanup -------------------------------------------------------
{
  echo "# ---- STAGE 2b truth set (bcftools consensus; replaces simuG) ----"
  echo "TRUTH_VCF=$TRUTH"; echo "CALLABLE_BED=$CALLABLE"
  echo "TRUTH_TOTAL=$NT"; echo "TRUTH_SNP=$NSNP"; echo "TRUTH_INDEL=$NIND"
  echo "TRUTH_TITV=$TITV_OUT"; echo "TRUTH_HET=$HET"; echo "TRUTH_HOM=$HOM"
  echo "HAP1_GENOME=$OUT/hap1.simseq.genome.fa"; echo "HAP2_GENOME=$OUT/hap2.simseq.genome.fa"
  echo "# STAGE2b_DONE $(date '+%F %T')"
} >> "$PARAMS"
rm -f "$OUT"/truth.raw.vcf "$OUT"/hap1.raw.vcf "$OUT"/hap2.raw.vcf
log "===== STAGE 2b DONE -> $TRUTH ====="
