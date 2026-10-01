#!/usr/bin/env python3
"""Figure 2: (a) caller rank under each truth definition (both samples); (b, c) SNP F1 with 95% bootstrap CIs.
Input: bootstrap_F1_CI.tsv and bootstrap_rank_stability.tsv (paired block bootstrap, 1,000 replicates)."""
import csv, sys, os, matplotlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import *
matplotlib.use('Agg'); import matplotlib.pyplot as plt
setup()
from matplotlib import font_manager
# usage: python3 fig2_rank_changes_and_F1_CI.py [tables_dir] [output_dir]
T = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TD = [('TS1_consensus', '3plus', 'Consensus\n≥3'), ('TS1_consensus', '4plus', 'Consensus\n≥4'), ('TS1_consensus', '5tools', 'Consensus\nall 5'),
      ('TS2_leave_one_out', 'k2', 'Leave-one-out\nk = 2'), ('TS2_leave_one_out', 'k3', 'Leave-one-out\nk = 3'),
      ('TS2_leave_one_out', 'k4', 'Leave-one-out\nk = 4'), ('TS4_simulation', 'known_truth', 'Simulation\n(known truth)')]
rank = {(r['TruthSet'], r['Level'], r['Sample'], r['Caller']): int(r['Rank']) for r in csv.DictReader(open(f'{T}/bootstrap_rank_stability.tsv'), delimiter='\t')}
ci = {(r['TruthSet'], r['Level'], r['Sample'], r['Caller']): (float(r['F1']), float(r['F1_lo']), float(r['F1_hi'])) for r in csv.DictReader(open(f'{T}/bootstrap_F1_CI.tsv'), delimiter='\t')}
smp = lambda ts, s: 'simulated' if ts == 'TS4_simulation' else s

fig = plt.figure(figsize=(180 / 25.4, 150 / 25.4))
gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.05], hspace=0.42)
gb = gs[1].subgridspec(1, 2, width_ratios=[6, 1.25], wspace=0.55)
ax = fig.add_subplot(gs[0])
import numpy as np
SHADE = RANK_SHADE
for i, c in enumerate(CALLERS):
    for k, (ts, lv, _) in enumerate(TD):
        rN = rank[(ts, lv, smp(ts, 'normal'), c)]; rA = rank[(ts, lv, smp(ts, 'abnormal'), c)]
        f = ci[(ts, lv, smp(ts, 'normal'), c)][0]
        ax.add_patch(plt.Rectangle((k, i), 1, 1, fc=SHADE[rN], ec='white', lw=1.5))
        txt = 'white' if rN <= 2 else '#1B1B1B'
        ax.text(k + 0.5, i + 0.36, f'#{rN}' + ('' if rA == rN else f' ({rA})'), ha='center', va='center', color=txt, fontsize=7.5, fontweight='bold')
        ax.text(k + 0.5, i + 0.70, f'{f:.4f}' if ts == 'TS4_simulation' else f'{f:.3f}', ha='center', va='center', color=txt, fontsize=6)
ax.set_xlim(0, len(TD)); ax.set_ylim(len(CALLERS), 0)
ax.set_yticks([i + 0.5 for i in range(5)]); ax.set_yticklabels([LABEL[c] for c in CALLERS])
ax.set_xticks([k + 0.5 for k in range(len(TD))]); ax.set_xticklabels([t.replace('Simulation\n(known truth)', 'Known\ntruth').replace('Leave-one-out\n', '').replace('Consensus\n', '') for *_, t in TD])
ax.xaxis.tick_top(); ax.tick_params(length=0); [sp.set_visible(False) for sp in ax.spines.values()]
for xm, t in ((1.5, 'Consensus (TS1)'), (4.5, 'Leave-one-out (TS2)'), (6.5, 'Simulation (TS3)')):
    ax.text(xm, -0.95, t, ha='center', va='bottom', fontsize=7, fontweight='bold', color=NAVY)
for b_ in (3, 6): ax.axvline(b_, color='white', lw=4)
ax.text(7, 5.25, 'Not separable (95% CI of the F1 difference includes 0): FreeBayes vs DeepVariant at k = 4; GATK vs VarScan in the simulation.\nAb-NormS2 ranks were identical in every column.', ha='right', va='top', fontsize=5.8, color='0.35')
ax.text(-0.13, 1.20, 'a', transform=ax.transAxes, fontweight='bold', fontsize=9)

def dots(ax, items, ylab=None):
    for k, (ts, lv, _) in enumerate(items):
        for j, c in enumerate(CALLERS):
            f, lo, hi = ci[(ts, lv, smp(ts, 'normal'), c)]
            xx = k + (j - 2) * 0.14
            ax.errorbar(xx, f, yerr=[[f - lo], [hi - f]], fmt=MK[c], color=COL[c], ms=3.5, elinewidth=0.9, capsize=1.5, mew=0)
    ax.set_xticks(range(len(items))); ax.set_xticklabels([t for *_, t in items]); ax.set_xlim(-0.6, len(items) - 0.4)
    ax.spines[['top', 'right']].set_visible(False)
    if ylab: ax.set_ylabel(ylab)
axb = fig.add_subplot(gb[0]); dots(axb, TD[:6], 'SNP F1 (95% CI)'); axb.set_ylim(0.88, 1.0)
axb.set_xticklabels(['≥3', '≥4', 'all 5', 'k = 2', 'k = 3', 'k = 4']); axb.axvline(2.5, color='0.75', lw=0.6, ls=':')
for xm, t in ((1, 'Consensus (TS1)'), (4, 'Leave-one-out (TS2)')):
    axb.text(xm, -0.13, t, transform=axb.get_xaxis_transform(), ha='center', va='top', fontsize=7, fontweight='bold', color=NAVY)
axb.text(-0.09, 1.05, 'b', transform=axb.transAxes, fontweight='bold', fontsize=9)
axb.set_title('Caller-derived truth sets', fontsize=7, pad=14)
hb = [plt.Line2D([], [], color=COL[c], marker=MK[c], ls='none', ms=4) for c in CALLERS]
axb.legend(hb, [LABEL[c] for c in CALLERS], ncol=5, frameon=False, loc='lower center', bbox_to_anchor=(0.5, 0.99), handletextpad=0.2, columnspacing=1.2, fontsize=6.5)
axc = fig.add_subplot(gb[1]); dots(axc, TD[6:]); axc.set_ylim(0.9960, 0.9995)
axc.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.4f'))
axc.set_title('Known truth\n(expanded scale)', fontsize=7)
axc.text(-0.75, 1.05, 'c', transform=axc.transAxes, fontweight='bold', fontsize=9)
save(fig, OUT, 'Fig2_rank_changes_and_F1_CI')
