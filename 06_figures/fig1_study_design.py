#!/usr/bin/env python3
"""Figure 1: study design. (a) workflow from reads to scoring; (b) how the three truth sets are built,
each as a short flow of steps ending in the truth set used for scoring.
usage: python3 fig1_study_design.py [output_dir]"""
import sys, os, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle, Polygon
from matplotlib.colors import LinearSegmentedColormap
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import NAVY, TEAL, PALE, CALLERS, LABEL, COL, setup, save
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
setup()
TXT, GREY, EDGE, OFF = '#1F2D36', '#5F6B7A', '#B8C2CE', '#E4E9F0'
W_MM, Y0, Y1 = 180, 0, 128
fig = plt.figure(figsize=(W_MM / 25.4, (Y1 - Y0) / 25.4)); ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W_MM); ax.set_ylim(Y0, Y1); ax.set_aspect('equal'); ax.axis('off')

def rbox(x0, y0, w, h, fc, ec='none', r=1.6, lw=0.7, z=1):
    ax.add_patch(FancyBboxPatch((x0, y0), w, h, boxstyle=f'round,pad=0,rounding_size={r}', fc=fc, ec=ec, lw=lw, zorder=z))
def t(x, y, s, fs=6.5, c=TXT, w='normal', ha='center', va='center', **k):
    ax.text(x, y, s, fontsize=fs, color=c, fontweight=w, ha=ha, va=va, zorder=5, **k)
def arrow(x1, y1, x2, y2, lw=0.9): ax.annotate('', xy=(x2, y2), xytext=(x1, y1), zorder=4,
    arrowprops=dict(arrowstyle='-|>', color=NAVY, lw=lw, mutation_scale=7, shrinkA=0, shrinkB=0))
def letter(x, y, s): t(x, y, s, fs=9, w='bold', ha='left', va='top', c='black')

# ================= (a) workflow =================
letter(1, Y1 - 0.5, 'a'); t(90, Y1 - 3, 'Identical processing and scoring for every caller and every truth set', fs=7.2, c=NAVY, w='bold')
steps = [('Reads', 'real: 2 turtles\nsimulated: 72 contigs'), ('Align', 'fastp, BWA-MEM,\nMarkDuplicates'),
         ('Call', 'BCFtools, DeepVariant,\nFreeBayes, GATK,\nVarScan'), ('Filter', 'PASS, ≥ 4 alt. reads,\nVAF > 2%'),
         ('Normalise', 'split, left-align,\nkeep SNPs'), ('Score', 'hap.py vs TS1–TS3\nF1 with 95% CI')]
BW, GAP, YB, BH = 25.5, 4.0, 104.5, 15.5
xs = [3.5 + BW / 2 + i * (BW + GAP) for i in range(6)]
for i, (x, (ti, su)) in enumerate(zip(xs, steps)):
    last = i == 5
    rbox(x - BW / 2, YB, BW, BH, NAVY if last else 'white', ec=NAVY if last else EDGE, lw=0.8)
    ax.add_patch(Circle((x - BW / 2 + 3.6, YB + BH - 3.6), 2.2, fc='white' if last else NAVY, ec='none', zorder=4))
    t(x - BW / 2 + 3.6, YB + BH - 3.6, str(i + 1), fs=6.3, w='bold', c=NAVY if last else 'white')
    t(x + 1.6, YB + BH - 3.6, ti, fs=7.2, w='bold', c='white' if last else NAVY)
    t(x, YB + 5.4, su, fs=5.9, c='white' if last else TXT, linespacing=1.25)
for a, b in zip(xs[:-1], xs[1:]): arrow(a + BW / 2 + 0.4, YB + BH / 2, b - BW / 2 - 0.4, YB + BH / 2)
t(90, YB - 3.2, 'Reference: rCheMyd1.pri.v2 (RefSeq GCF_015237465.2, 2.13 Gb).  Real reads: NovaSeq 6000, 2 × 151 bp; simulated reads: ART, 2 × 151 bp, 70×',
  fs=5.8, c=GREY)

# ================= (b) truth sets =================
letter(1, 96.5, 'b'); t(90, 94, 'How the three truth sets were built', fs=7.2, c=NAVY, w='bold')
CW, CG, CY0, CY1 = 55.0, 3.5, 22, 90
cx = [3.5 + CW / 2 + i * (CW + CG) for i in range(3)]
heads = [('TS1  Consensus', 'every caller votes, including the one scored'),
         ('TS2  Leave-one-out', 'the caller being scored has no vote'),
         ('TS3  Simulation', 'variants known before any caller is run')]
for c, (h1, h2) in zip(cx, heads):
    rbox(c - CW / 2, CY0, CW, CY1 - CY0, '#F7F9FB', ec=EDGE, lw=0.8)
    rbox(c - CW / 2, CY1 - 8.5, CW, 8.5, TEAL, r=1.6, z=2)
    ax.add_patch(Rectangle((c - CW / 2, CY1 - 8.5), CW, 2, fc=TEAL, ec='none', zorder=2))   # square lower corners
    t(c, CY1 - 4.25, h1, fs=7.5, w='bold', c='white')
    t(c, CY1 - 11.8, h2, fs=6, c=GREY, style='italic')

SY = [72.0, 63.0, 54.0, 45.0]; SW, SH = 47, 7.0
def dots(c, y, out=None):                                   # five caller dots; the left-out caller faded and crossed
    for j, k in enumerate(CALLERS):
        x = c + (j - 2) * 5.2
        ax.add_patch(Circle((x, y), 1.25, fc=COL[k], ec='none', alpha=0.25 if k == out else 1, zorder=4))
        if k == out:
            for d in (1, -1): ax.plot([x - 1.3, x + 1.3], [y - 1.3 * d, y + 1.3 * d], color='0.35', lw=0.8, zorder=5)
def flow(c, steps, final, note, out=None):
    for n, (y, (a1, a2)) in enumerate(zip(SY, steps)):
        rbox(c - SW / 2, y - SH / 2, SW, SH, 'white', ec=EDGE, lw=0.7, r=1.2)
        t(c, y + 1.55, a1, fs=6.4, w='bold', c=NAVY)
        if a2 == 'DOTS': dots(c, y - 1.75, out)
        else: t(c, y - 1.75, a2, fs=5.9)
    for a, b in zip(SY[:-1], SY[1:]): arrow(c, a - SH / 2 - 0.2, c, b + SH / 2 + 0.2, lw=0.8)
    rbox(c - SW / 2, 33.2, SW, 6.0, NAVY, r=1.2); t(c, 36.2, final, fs=6.4, w='bold', c='white')
    arrow(c, SY[-1] - SH / 2 - 0.2, c, 39.4, lw=0.8)
    t(c, CY0 + 4.2, note, fs=6.0, c=TXT, linespacing=1.3)
flow(cx[0], [('Five filtered SNP call sets', 'DOTS'), ('Count callers per position', 'bcftools isec over all five'),
             ('Keep agreed positions', '≥ 3, ≥ 4 or all 5 callers'), ('Truth sizes (CH-NORMS1)', '12.40 M, 11.88 M, 10.28 M SNPs')],
     'Consensus truth (3 levels)', 'True = reported by enough callers; the\ncaller being scored also contributes')
flow(cx[1], [('Leave out the tested caller', 'DOTS'), ('Count the other four', 'bcftools isec over four call sets'),
             ('Keep agreed positions', '≥ k of the other four (k = 2, 3, 4)'), ('Repeat for every caller', '5 callers × 3 levels')],
     'Leave-one-out truth (per caller)', 'True = reported by the other callers; the\ncaller being scored has no influence', out='Freebayes')
flow(cx[2], [('Reference subset', '72 whole contigs, 211.7 Mb'), ('Plant known variants', '1,227,026 SNVs, 153,429 indels'),
             ('Simulate reads (ART)', 'real error profile, 70×'), ('Align, call and filter', 'exactly as for the real genomes')],
     'Known truth (callable regions)', 'True = planted in the genome; does not\ndepend on any caller')

# independence gradient arrow
cmap = LinearSegmentedColormap.from_list('b', ['#DCE8F5', NAVY])
ax.imshow(np.linspace(0, 1, 256)[None, :], extent=(14, 166, 14.2, 16.6), aspect='auto', cmap=cmap, zorder=2)
ax.add_patch(Polygon([[166, 13.2], [170.5, 15.4], [166, 17.6]], closed=True, fc=NAVY, ec='none', zorder=2))
t(11.5, 15.4, 'TS1', fs=6.2, w='bold', c=NAVY, ha='right'); t(172, 15.4, 'TS3', fs=6.2, w='bold', c=NAVY, ha='left')
t(90, 11.0, 'Increasing independence of the truth from the callers being compared', fs=6.2, c=NAVY, style='italic')
ax.set_xlim(0, W_MM); ax.set_ylim(8.5, Y1); fig.set_size_inches(W_MM / 25.4, (Y1 - 8.5) / 25.4)
save(fig, OUT, 'Fig1_study_design')
