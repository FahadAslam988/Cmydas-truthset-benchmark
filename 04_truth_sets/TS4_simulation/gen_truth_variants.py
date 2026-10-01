#!/usr/bin/env python3
"""
Generate a KNOWN truth set of variants for the C. mydas simulation benchmark.
Fully specified + reproducible (fixed seed); every distribution calibrated to the
real data. Emits three VCFs (native reference coordinates):
  truth.vcf  : all variants with REAL genotypes (het 0/1, hom 1/1)  <- the answer key
  hap1.vcf   : ALL variants as 1/1  (-> bcftools consensus builds haplotype 1)
  hap2.vcf   : HOM variants as 1/1  (-> bcftools consensus builds haplotype 2)

Method (matches Hall et al. 2024, eLife, who apply a truthset with bcftools consensus):
  * positions sampled uniformly across NON-N sequence, one per length/k bin so
    variants are sorted and never overlap (spacing >> max indel length);
  * SNP alt allele respects the measured Ti/Tv ratio;
  * indels: ins:del = 1:1, short-length spectrum (geometric, capped);
  * genotype het:hom set to the measured HET_FRACTION.
Usage: gen_truth_variants.py SUBSET_FA N_SNP N_INDEL TITV HET_FRACTION SEED OUTDIR
"""
import sys, random

SUBSET_FA, N_SNP, N_INDEL, TITV, HETF, SEED, OUTDIR = (
    sys.argv[1], int(sys.argv[2]), int(sys.argv[3]),
    float(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6]), sys.argv[7])

random.seed(SEED)
BASES = "ACGT"
TS = {"A": "G", "G": "A", "C": "T", "T": "C"}          # transitions
TV = {"A": ["C", "T"], "G": ["C", "T"], "C": ["A", "G"], "T": ["A", "G"]}
P_TS = TITV / (TITV + 1.0)                               # P(transition)
P_INDEL = N_INDEL / float(N_SNP + N_INDEL)
MAX_INDEL = 10
MARGIN = MAX_INDEL + 2                                   # keep indels inside intervals

# ---- load reference (uppercase) -------------------------------------------
def load_fa(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line[0] == ">":
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
        if name is not None:
            seqs[name] = "".join(buf).upper()
    return seqs

sys.stderr.write("loading reference...\n")
seqs = load_fa(SUBSET_FA)
order = list(seqs.keys())

# ---- non-N callable intervals per contig (trim MARGIN at ends) ------------
def callable_intervals(seq):
    ivs, i, n = [], 0, len(seq)
    while i < n:
        if seq[i] in "ACGT":
            j = i
            while j < n and seq[j] in "ACGT":
                j += 1
            if j - i > 2 * MARGIN:
                ivs.append((i + MARGIN, j - MARGIN))     # [start, end)
            i = j
        else:
            i += 1
    return ivs

intervals = {c: callable_intervals(seqs[c]) for c in order}
clen = {c: sum(e - s for s, e in intervals[c]) for c in order}
total_callable = sum(clen.values())
sys.stderr.write(f"callable (non-N, trimmed) bp: {total_callable}\n")

# ---- allocate variant counts per contig, proportional to callable length --
TOTAL = N_SNP + N_INDEL
def alloc(c):
    return int(round(TOTAL * clen[c] / total_callable)) if total_callable else 0

# map a callable-offset (0..clen-1) to a genomic 0-based position
def offset_to_pos(ivs, off):
    for s, e in ivs:
        L = e - s
        if off < L:
            return s + off
        off -= L
    return None

# ---- generate ------------------------------------------------------------
truth, hap1, hap2 = [], [], []      # each entry: (chrom,pos,ref,alt,gt)
n_snp = n_ind = n_het = n_hom = 0
for c in order:
    seq = seqs[c]; ivs = intervals[c]; L = clen[c]
    k = alloc(c)
    if k == 0 or L == 0:
        continue
    binsz = L / k
    last = -10
    for i in range(k):
        lo = int(i * binsz); hi = max(lo + 1, int((i + 1) * binsz))
        off = random.randrange(lo, min(hi, L))
        pos0 = offset_to_pos(ivs, off)
        if pos0 is None or pos0 <= last + MARGIN:
            continue
        rbase = seq[pos0]
        if rbase not in "ACGT":
            continue
        is_indel = (random.random() < P_INDEL)
        if not is_indel:
            # SNP honoring Ti/Tv
            if random.random() < P_TS:
                alt = TS[rbase]
            else:
                alt = random.choice(TV[rbase])
            ref = rbase; n_snp += 1
        else:
            # indel: 1:1 ins:del, short geometric length (1..MAX_INDEL)
            ln = 1
            while random.random() < 0.5 and ln < MAX_INDEL:
                ln += 1
            if random.random() < 0.5:                    # insertion
                ins = "".join(random.choice(BASES) for _ in range(ln))
                ref = rbase; alt = rbase + ins
            else:                                        # deletion
                span = seq[pos0:pos0 + ln + 1]
                if len(span) < ln + 1 or any(b not in "ACGT" for b in span):
                    continue
                ref = span; alt = rbase
            n_ind += 1
        # genotype
        if random.random() < HETF:
            gt = "0/1"; n_het += 1
        else:
            gt = "1/1"; n_hom += 1
        last = pos0
        vpos = pos0 + 1                                  # VCF is 1-based
        truth.append((c, vpos, ref, alt, gt))
        hap1.append((c, vpos, ref, alt, "1/1"))          # hap1 carries ALL
        if gt == "1/1":
            hap2.append((c, vpos, ref, alt, "1/1"))      # hap2 carries HOM only

# ---- write VCFs ----------------------------------------------------------
def header(fh):
    fh.write("##fileformat=VCFv4.2\n")
    fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    for c in order:
        fh.write(f"##contig=<ID={c},length={len(seqs[c])}>\n")
    fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSIMSAMPLE\n")

def write_vcf(path, rows):
    with open(path, "w") as fh:
        header(fh)
        for c, p, r, a, gt in rows:
            fh.write(f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\t.\tGT\t{gt}\n")

write_vcf(f"{OUTDIR}/truth.raw.vcf", truth)
write_vcf(f"{OUTDIR}/hap1.raw.vcf",  hap1)
write_vcf(f"{OUTDIR}/hap2.raw.vcf",  hap2)

sys.stderr.write(
    f"DONE: variants={len(truth)} SNP={n_snp} INDEL={n_ind} "
    f"het={n_het} hom={n_hom} het_frac={n_het/max(1,len(truth)):.3f}\n")
