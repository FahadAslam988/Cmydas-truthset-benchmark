#!/usr/bin/env python3
"""Figure 3: why rankings change. Recall (a, c) and precision (b, d) of each caller as the caller-derived truth
becomes stricter (consensus >=3 -> >=4 -> all 5; leave-one-out k = 2 -> 3 -> 4), 95% bootstrap CIs;
(e) share of disputed sites (called by exactly 4 or exactly 3 of the 5 callers) reported by each caller.
usage: python3 fig3_mechanism.py [tables_dir] [output_dir]"""
import csv, sys, os, matplotlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import *
matplotlib.use('Agg'); import matplotlib.pyplot as plt
setup()
from matplotlib import font_manager
T = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SMP = os.environ.get('FIG_SAMPLE', 'normal'); NAME = 'Fig3_mechanism' if SMP == 'normal' else 'FigS1_mechanism_Ab-NormS2'
ci = {(r['TruthSet'], r['Level'], r['Sample'], r['Caller']): r for r in csv.DictReader(open(f'{T}/bootstrap_F1_CI.tsv'), delimiter='\t')}
met = {(r['TruthSet'], r['Level'], r['Sample'], r['Caller']): r for r in csv.DictReader(open(f'{T}/Paper1_all_truthsets_SNP_metrics.tsv'), delimiter='\t')}
FAM = [('TS1_consensus', ['3plus', '4plus', '5tools'], ['≥3', '≥4', 'All 5'], 'Consensus truth (callers agreeing)'),
       ('TS2_leave_one_out', ['k2', 'k3', 'k4'], ['k = 2', 'k = 3', 'k = 4'], 'Leave-one-out truth (other callers agreeing)')]
MET = [('Recall', 'R_lo', 'R_hi', 'Recall'), ('Precision', 'P_lo', 'P_hi', 'Precision')]

fig = plt.figure(figsize=(180 / 25.4, 112 / 25.4))
outer = fig.add_gridspec(1, 2, width_ratios=[2.1, 1], wspace=0.32)
left = outer[0].subgridspec(2, 2, hspace=0.62, wspace=0.42)
letters = iter('abcd')
for row, (ts, lvls, xl, title) in enumerate(FAM):
    for col, (m, lo, hi, ylab) in enumerate(MET):
        ax = fig.add_subplot(left[row, col])
        for j, c in enumerate(CALLERS):
            xs = [i + (j - 2) * 0.06 for i in range(3)]
            v = [float(ci[(ts, l, SMP, c)][m]) for l in lvls]
            e = [[v[i] - float(ci[(ts, l, SMP, c)][lo]) for i, l in enumerate(lvls)],
                 [float(ci[(ts, l, SMP, c)][hi]) - v[i] for i, l in enumerate(lvls)]]
            ax.errorbar(xs, v, yerr=e, color=COL[c], marker=MK[c], ms=3.6, lw=1.2, elinewidth=0.7, capsize=0, label=LABEL[c])
        ax.set_xticks(range(3)); ax.set_xticklabels(xl); ax.set_xlim(-0.35, 2.35)
        ax.set_ylim(0.80, 1.005); ax.set_yticks([0.80, 0.85, 0.90, 0.95, 1.00]); ax.set_ylabel(ylab)
        ax.spines[['top', 'right']].set_visible(False); ax.grid(axis='y', color='0.92', lw=0.6); ax.set_axisbelow(True)
        ax.text(-0.36, 1.08, next(letters), transform=ax.transAxes, fontweight='bold', fontsize=9)
        if col == 0:
            ax.text(1.18, 1.17, title, transform=ax.transAxes, ha='center', fontsize=7, fontweight='bold', color=NAVY)
        if row == 1: ax.set_xlabel('Stricter truth  →')
ax5 = fig.add_subplot(outer[1])
w = 0.38
for j, (lv, hat) in enumerate((('n4', None), ('n3', '////'))):
    for i, c in enumerate(CALLERS):
        v = 100 * float(met[('TS1_consensus', lv, SMP, c)]['Recall'])
        ax5.bar(i + (j - 0.5) * w, v, w * 0.95, color=COL[c] if j == 0 else 'white', edgecolor=COL[c], hatch=hat, lw=0.8)
ax5.set_xticks(range(5)); ax5.set_xticklabels([LABEL[c] for c in CALLERS], rotation=40, ha='right')
ax5.set_ylim(0, 112); ax5.set_yticks(range(0, 101, 20)); ax5.set_ylabel('Disputed positions called (%)')
ax5.spines[['top', 'right']].set_visible(False); ax5.grid(axis='y', color='0.92', lw=0.6); ax5.set_axisbelow(True)
ax5.text(0.5, 1.035, 'Disputed SNP positions', transform=ax5.transAxes, ha='center', fontsize=7, fontweight='bold', color=NAVY)
ax5.text(-0.30, 1.035, 'e', transform=ax5.transAxes, fontweight='bold', fontsize=9)
hs = [plt.Rectangle((0, 0), 1, 1, fc='0.5', ec='0.5'), plt.Rectangle((0, 0), 1, 1, fc='white', ec='0.5', hatch='////')]
ax5.legend(hs, ['Called by exactly 4 of 5', 'Called by exactly 3 of 5'], frameon=False, loc='upper right', bbox_to_anchor=(1.03, 1.03), ncol=1, fontsize=6.5, handlelength=1.3)
h, l = fig.axes[0].get_legend_handles_labels()
fig.legend(h, l, ncol=5, frameon=False, loc='upper center', bbox_to_anchor=(0.5, 1.03), handlelength=1.8, columnspacing=1.2)
save(fig, OUT, NAME)
