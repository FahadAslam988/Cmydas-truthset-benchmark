#!/usr/bin/env bash
source "${CMYDAS_CONFIG:?export CMYDAS_CONFIG=/path/to/config/paths.sh}"
# STEP 6: 95% confidence intervals by PAIRED BLOCK BOOTSTRAP (1,000 replicates, 1-Mb blocks, seed 20261001).
# Input: the per-variant hap.py output VCFs (FORMAT/BD = TP/FP/FN for TRUTH and QUERY; SNPs only).
# A) per-block counts per run; B) resample blocks (same blocks for all callers of one truth set) ->
#    CI of precision/recall/F1, CI of pairwise F1 differences, rank stability.
set -uo pipefail; source "$(dirname "$0")/common.sh"
D=$E/06_bootstrap_ci; mkdir -p $D/blocks
L=$D/runs.tsv; : > $L
for g in normal abnormal; do for cl in "${CALLERS[@]}"; do
  for lv in 3plus 4plus 5tools; do echo -e "TS1_consensus\t$lv\t$g\t$cl\t$E/02_TS1_consensus_happy/TS1_${lv}_${g}_${cl}.vcf.gz" >> $L; done
  echo -e "TS2_leave_one_out\tk2\t$g\t$cl\t$R/06_loo/happy/${g}_${cl}.vcf.gz" >> $L
  echo -e "TS2_leave_one_out\tk3\t$g\t$cl\t$R/07_loo_sensitivity/k3_ge3_of4/happy/${g}_${cl}.vcf.gz" >> $L
  echo -e "TS2_leave_one_out\tk4\t$g\t$cl\t$R/07_loo_sensitivity/k4_eq4_of4/happy/${g}_${cl}.vcf.gz" >> $L
done; done
for cl in "${CALLERS[@]}"; do echo -e "TS4_simulation\tknown_truth\tsimulated\t$cl\t$E/03_TS4_simulation_happy_allele/TS4_simulation_${cl}.vcf.gz" >> $L; done
for lv in ge2 ge3 ge4 all5; do echo -e "TS4_simulation\tknown_truth\tsimulated\tconsensus_$lv\t$E/05_TS4_multicaller_consensus/TS4_multicaller_${lv}.vcf.gz" >> $L; done
# A) per-block counts: block, truthTP, FN, queryTP, FP
export BCF D
cut -f5 $L | xargs -P 16 -I{} bash -c 'f={}; o=$D/blocks/$(echo "$f" | md5sum | cut -c1-16).tsv; [ -s $o ] && exit 0
  $BCF query -f "%CHROM\t%POS[\t%BD\t%BVT]\n" "$f" | awk -F"\t" -v OFS="\t" "{b=\$1\":\"int(\$2/1000000)
    if(\$4==\"SNP\"){if(\$3==\"TP\")tt[b]++; else if(\$3==\"FN\")fn[b]++}
    if(\$6==\"SNP\"){if(\$5==\"TP\")qt[b]++; else if(\$5==\"FP\")fp[b]++}; k[b]=1}
    END{for(b in k)print b,tt[b]+0,fn[b]+0,qt[b]+0,fp[b]+0}" > $o.tmp && mv $o.tmp $o'
# B) bootstrap
python3 - "$L" "$D" <<'PY'
import sys,csv,hashlib,itertools,numpy as np,collections
L,D=sys.argv[1],sys.argv[2]; rng=np.random.default_rng(20261001); NB=1000
runs=collections.defaultdict(dict)
for ts,lv,s,cl,f in csv.reader(open(L),delimiter='\t'):
    h=hashlib.md5((f+'\n').encode()).hexdigest()[:16]
    runs[(ts,lv,s)][cl]={r[0]:np.array(r[1:],float) for r in csv.reader(open(f'{D}/blocks/{h}.tsv'),delimiter='\t')}
def f1(c):  # c[...,4] = tt,fn,qt,fp
    p=c[...,2]/(c[...,2]+c[...,3]); r=c[...,0]/(c[...,0]+c[...,1]); return p,r,2*p*r/(p+r)
o1=open(f'{D}/bootstrap_F1_CI.tsv','w'); o1.write('TruthSet\tLevel\tSample\tCaller\tPrecision\tP_lo\tP_hi\tRecall\tR_lo\tR_hi\tF1\tF1_lo\tF1_hi\tTP\tFN\tFP\n')
o2=open(f'{D}/bootstrap_pairwise_F1_diff.tsv','w'); o2.write('TruthSet\tLevel\tSample\tCallerA\tCallerB\tF1_A_minus_B\tlo\thi\tShare_A_better\n')
o3=open(f'{D}/bootstrap_rank_stability.tsv','w'); o3.write('TruthSet\tLevel\tSample\tCaller\tRank\tShare_rank1\tMean_rank\tRank_lo\tRank_hi\n')
for key,cl in runs.items():
    names=list(cl); blocks=sorted(set().union(*[set(v) for v in cl.values()]))
    M=np.array([[cl[n].get(b,np.zeros(4)) for b in blocks] for n in names])      # callers x blocks x 4
    idx=rng.integers(0,len(blocks),(NB,len(blocks)))
    W=np.stack([np.bincount(i,minlength=len(blocks)) for i in idx])            # NB x blocks
    tot=M.sum(1); bs=np.einsum('kb,nbc->knc',W,M)                               # NB x callers x 4
    P,Rr,F=f1(tot); bp,br,bf=f1(bs)
    for j,n in enumerate(names):
        q=lambda a:np.percentile(a[:,j],[2.5,97.5])
        o1.write('\t'.join(key)+f'\t{n}\t{P[j]:.6f}\t'+'\t'.join(f'{x:.6f}' for x in q(bp))+f'\t{Rr[j]:.6f}\t'+'\t'.join(f'{x:.6f}' for x in q(br))
                 +f'\t{F[j]:.6f}\t'+'\t'.join(f'{x:.6f}' for x in q(bf))+f'\t{int(tot[j,0])}\t{int(tot[j,1])}\t{int(tot[j,3])}\n')
    for a,b in itertools.combinations(range(len(names)),2):
        d=bf[:,a]-bf[:,b]; lo,hi=np.percentile(d,[2.5,97.5])
        o2.write('\t'.join(key)+f'\t{names[a]}\t{names[b]}\t{F[a]-F[b]:.6f}\t{lo:.6f}\t{hi:.6f}\t{(d>0).mean():.3f}\n')
    sing=[j for j,n in enumerate(names) if not n.startswith('consensus_')]
    rk=(-bf[:,sing]).argsort(1).argsort(1)+1; r0=(-F[sing]).argsort().argsort()+1
    for k,j in enumerate(sing):
        o3.write('\t'.join(key)+f'\t{names[j]}\t{r0[k]}\t{(rk[:,k]==1).mean():.3f}\t{rk[:,k].mean():.2f}\t'+'\t'.join(str(int(x)) for x in np.percentile(rk[:,k],[2.5,97.5]))+'\n')
print('groups:',len(runs))
PY
log "bootstrap done: $D"
