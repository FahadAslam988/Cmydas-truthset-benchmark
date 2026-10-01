#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# ============================================================================
# STAGE 5b — Traditional callers on the simulated BAM, best-practice + parallel.
# Runs ONE BY ONE: GATK -> BCFtools -> FreeBayes -> VarScan. Each is GUARDED:
# if a tool fails, it is logged and we SKIP to the next tool (no disturbance).
# Versions matched via conda envs: GATK 4.6.2.0, BCFtools 1.22 (samtools 1.22.1),
# FreeBayes 1.3.10, VarScan 2.4.6. Each caller parallelized across the 72 contigs.
# Resumable (skips a caller whose output VCF already exists). Logged.
# ============================================================================
set -uo pipefail
BASE="${SIM_DIR}"
REF="${REF_FASTA}"
BAM="$BASE/outputs/03_aligned/sim.dedup.bam"
OUT="$BASE/outputs/04_calls"
GDIR="$BASE/outputs/01_mutated_genome"
LOG="$BASE/run.log"; PARAMS="$BASE/params.txt"
PZ=48                                  # parallel jobs per caller (leaves headroom)
BCFT=bcftools
SAMT=samtools   # samtools 1.22.1
mapfile -t CONTIGS < <(cut -f1 "$GDIR/subset.fa.fai")       # 72 read-bearing contigs
export REF BAM
log(){ echo "[$(date '+%F %T')] [s06b] $*" | tee -a "$LOG"; }
source ${CONDA_SH}

# guard: run a caller function; on failure log + CONTINUE to next tool
guard(){ local name="$1"; shift; log ">>> START $name";
  if "$@"; then log "<<< DONE $name"; else log "!!! FAILED $name (rc=$?) -> skipping to next tool"; fi; }

concat_sorted(){ # $1=out.vcf.gz  $2..=per-contig vcf.gz (in contig order)
  local o="$1"; shift
  "$BCFT" concat "$@" 2>>"$LOG" | "$BCFT" sort -Oz -o "$o" 2>>"$LOG" && "$BCFT" index -t "$o" 2>>"$LOG"
}

# ---------------- GATK: scatter HaplotypeCaller(GVCF)->GenotypeGVCFs -> gather
run_gatk(){
  local o="$OUT/Gatk_sim.vcf.gz"; [ -s "$o" ] && { log "Gatk exists, skip"; return 0; }
  set +u; conda activate gatk-env; set -u
  local wd="$OUT/.gatk_tmp"; mkdir -p "$wd"; : > "$wd/fail"
  printf '%s\n' "${CONTIGS[@]}" | xargs -P "$PZ" -I{} bash -c '
    c="$1"; wd="'"$wd"'"
    gatk --java-options "-Xmx8g" HaplotypeCaller -R "$REF" -I "$BAM" -L "$c" -ERC GVCF \
         -O "$wd/$c.g.vcf.gz" >/dev/null 2>>"$wd/$c.err" \
    && gatk --java-options "-Xmx8g" GenotypeGVCFs -R "$REF" -V "$wd/$c.g.vcf.gz" -L "$c" \
         -O "$wd/$c.gt.vcf.gz" >/dev/null 2>>"$wd/$c.err" || echo "$c" >>"$wd/fail"
  ' _ {}
  [ -s "$wd/fail" ] && log "GATK: failed contigs: $(tr '\n' ' ' <"$wd/fail")"
  local files=(); for c in "${CONTIGS[@]}"; do [ -s "$wd/$c.gt.vcf.gz" ] && files+=("$wd/$c.gt.vcf.gz"); done
  concat_sorted "$o" "${files[@]}"
  set +u; conda deactivate; set -u; rm -rf "$wd"; [ -s "$o" ]
}

# ---------------- BCFtools: mpileup(AD,DP)->call -mv, per contig -------------
run_bcf(){
  local o="$OUT/BCF_sim.vcf.gz"; [ -s "$o" ] && { log "BCF exists, skip"; return 0; }
  local wd="$OUT/.bcf_tmp"; mkdir -p "$wd"
  printf '%s\n' "${CONTIGS[@]}" | xargs -P "$PZ" -I{} bash -c '
    c="$1"; wd="'"$wd"'"; B="'"$BCFT"'"
    "$B" mpileup -f "$REF" -a FORMAT/AD,FORMAT/DP -r "$c" "$BAM" 2>/dev/null \
      | "$B" call -mv -Oz -o "$wd/$c.vcf.gz" 2>/dev/null
  ' _ {}
  local files=(); for c in "${CONTIGS[@]}"; do [ -s "$wd/$c.vcf.gz" ] && files+=("$wd/$c.vcf.gz"); done
  concat_sorted "$o" "${files[@]}"
  rm -rf "$wd"; [ -s "$o" ]
}

# ---------------- FreeBayes: freebayes-parallel across region chunks ---------
run_fb(){
  local o="$OUT/Freebayes_sim.vcf.gz"; [ -s "$o" ] && { log "FreeBayes exists, skip"; return 0; }
  set +u; conda activate freebayes-env; set -u
  local reg="$OUT/.fb_regions"
  fasta_generate_regions.py "$GDIR/subset.fa.fai" 2000000 > "$reg" 2>>"$LOG"
  freebayes-parallel "$reg" "$PZ" -f "$REF" "$BAM" > "$OUT/.fb.vcf" 2>>"$LOG"
  set +u; conda deactivate; set -u
  "$BCFT" sort "$OUT/.fb.vcf" -Oz -o "$o" 2>>"$LOG" && "$BCFT" index -t "$o" 2>>"$LOG"
  rm -f "$OUT/.fb.vcf" "$reg"; [ -s "$o" ]
}

# ---------------- VarScan: samtools mpileup | varscan mpileup2cns ------------
run_varscan(){
  local o="$OUT/Varscan_sim.vcf.gz"; [ -s "$o" ] && { log "VarScan exists, skip"; return 0; }
  set +u; conda activate varscan-env; set -u
  local wd="$OUT/.vs_tmp"; mkdir -p "$wd"
  # VarScan omits ##contig headers -> build them from the reference for the merge
  awk '{print "##contig=<ID="$1",length="$2">"}' "$REF.fai" > "$wd/contigs.hdr"
  printf '%s\n' "${CONTIGS[@]}" | xargs -P "$PZ" -I{} bash -c '
    c="$1"; wd="'"$wd"'"; S="'"$SAMT"'"
    "$S" mpileup -f "$REF" -r "$c" "$BAM" 2>/dev/null \
      | varscan mpileup2cns --variants 1 --output-vcf 1 > "$wd/$c.vcf" 2>/dev/null
    /usr/bin/bgzip -f "$wd/$c.vcf" 2>/dev/null && /usr/bin/tabix -f -p vcf "$wd/$c.vcf.gz" 2>/dev/null
  ' _ {}
  set +u; conda deactivate; set -u
  local files=(); for c in "${CONTIGS[@]}"; do [ -s "$wd/$c.vcf.gz" ] && files+=("$wd/$c.vcf.gz"); done
  [ ${#files[@]} -eq 0 ] && { log "VarScan: no per-contig VCFs produced"; rm -rf "$wd"; return 3; }
  # concat (indexed) -> inject ##contig headers -> sort -> index  (the fix)
  "$BCFT" concat -a "${files[@]}" 2>>"$LOG" \
    | "$BCFT" annotate -h "$wd/contigs.hdr" 2>>"$LOG" \
    | "$BCFT" sort -Oz -o "$o" 2>>"$LOG"
  "$BCFT" index -t "$o" 2>>"$LOG"
  rm -rf "$wd"; [ -s "$o" ]
}

log "===== STAGE 5b traditional callers START (one-by-one, guarded, PZ=$PZ) ====="
guard "GATK"      run_gatk
guard "BCFtools"  run_bcf
guard "FreeBayes" run_fb
guard "VarScan"   run_varscan
# summary
log "caller outputs:"; for t in Gatk BCF Freebayes Varscan; do
  f="$OUT/${t}_sim.vcf.gz"
  if [ -s "$f" ]; then log "  $t: $(bcftools view -H "$f" 2>/dev/null | wc -l) records"; else log "  $t: MISSING"; fi
done
echo "# STAGE5b_DONE $(date '+%F %T')" >> "$PARAMS"
log "===== STAGE 5b DONE ====="
