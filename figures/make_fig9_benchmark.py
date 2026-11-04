# Fig. 9: known-signal benchmark in ten simulated datasets.
# Same style as the other multi-panel figures of the paper: Times New Roman, 12 pt, bold,
# drawn on a 13-inch-wide canvas (scaled to the page width in the manuscript).
# Usage: python make_fig9_benchmark.py per_ranking_metrics.csv
import sys
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.lines import Line2D

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman', 'Liberation Serif', 'Tinos', 'DejaVu Serif']
rcParams['font.size'] = 12
rcParams['font.weight'] = 'bold'
rcParams['axes.labelweight'] = 'bold'
rcParams['axes.titleweight'] = 'bold'
rcParams['pdf.fonttype'] = 42

d = pd.read_csv(sys.argv[1])
ref = d[d.condition == 'Reference']; nr = d[d.condition != 'Reference'].copy()
nr['group'] = np.where(nr.condition == 'Observed', nr.generator, nr.condition)
COL = {'ARF': '#F8766D', 'CTGAN': '#00BFC4', 'TVAE': '#E5A500',
       'Control A': '#BDBDBD', 'Control B': '#7F7F7F', 'Control C': '#B07CC6'}
order = ['ARF', 'CTGAN', 'TVAE', 'Control A', 'Control B', 'Control C']
xlabels = ['ARF', 'CTGAN', 'TVAE', 'Control\nA', 'Control\nB', 'Control\nC']
fig, axs = plt.subplots(1, 2, figsize=(13, 5.8), gridspec_kw={'width_ratios': [1, 1]})
rng = np.random.default_rng(3)

ax = axs[0]
for i, g in enumerate(order):
    y = nr[nr.group == g].tp50_main.values; c = COL[g]
    ax.boxplot([y], positions=[i], widths=0.62, whis=(0, 100), showfliers=False, patch_artist=True,
               boxprops=dict(facecolor=c, edgecolor='black', linewidth=0.8, alpha=0.85),
               medianprops=dict(color='black', linewidth=1.6),
               whiskerprops=dict(color='black', linewidth=0.8), capprops=dict(color='black', linewidth=0.8))
    if g in ('ARF', 'CTGAN', 'TVAE', 'Control C'):
        ax.scatter(i + rng.uniform(-0.17, 0.17, len(y)), y, s=14, color='black', alpha=0.55, linewidths=0, zorder=3)
rm = ref.tp50_main.mean()
ax.axhline(rm, ls='--', color='black', lw=1)
ax.text(-0.55, rm + 0.5, 'reference mean', ha='left', va='bottom', fontsize=12, fontweight='bold')
ax.set_xticks(range(6)); ax.set_xticklabels(xlabels, fontsize=12, fontweight='bold')
ax.set_xlim(-0.6, 5.6)
ax.set_ylabel('True main-effect genes in the Top-50', fontsize=12, fontweight='bold')
ax.set_title('(a)', fontsize=12, fontweight='bold')

ax = axs[1]
# Unit of every point: one RF ranking (not averaged within simulations).
# Drawing order: Control B (300 rankings), Control A (30), Control C (30), then ARF/CTGAN/TVAE (10 each) on top.
sty = {'Control B': dict(marker='o', s=16, alpha=0.75, edgecolors='none', linewidths=0, zorder=2),
       'Control A': dict(marker='s', s=40, alpha=1.0, edgecolors='black', linewidths=0.8, zorder=3),
       'Control C': dict(marker='o', s=46, alpha=0.95, edgecolors='black', linewidths=0.7, zorder=4),
       'ARF': dict(marker='o', s=46, alpha=0.95, edgecolors='black', linewidths=0.7, zorder=5),
       'CTGAN': dict(marker='o', s=46, alpha=0.95, edgecolors='black', linewidths=0.7, zorder=5),
       'TVAE': dict(marker='o', s=46, alpha=0.95, edgecolors='black', linewidths=0.7, zorder=5)}
for g in ['Control B', 'Control A', 'Control C', 'ARF', 'CTGAN', 'TVAE']:
    s = nr[nr.group == g]
    ax.scatter(s.mean_abs_shift, s.loss_tp, color=COL[g], **sty[g])
ax.axhline(0, ls=':', color='black', lw=1, zorder=1)   # 0 = same number of true genes as in the reference Top-50
ax.set_xlabel('Mean absolute rank shift of reference Top-50 genes', fontsize=12, fontweight='bold')
ax.set_ylabel('Loss of true main-effect genes from the Top-50', fontsize=12, fontweight='bold')
ax.set_title('(b)', fontsize=12, fontweight='bold')
# boxed colour legend at the top right, as in the other figures of the paper
# legend symbols follow panel (b): a square for Control A, circles for all other groups
handles = [Line2D([0], [0], marker='s' if g == 'Control A' else 'o', ls='', markersize=9,
                  markerfacecolor=COL[g], markeredgecolor='black', markeredgewidth=0.8, label=g) for g in order]
fig.legend(handles=handles, loc='upper right', ncol=6, frameon=True, facecolor='white', edgecolor='black',
           fancybox=False, prop={'weight': 'bold', 'size': 12}, handlelength=1.0, handletextpad=0.3, columnspacing=1.2,
           bbox_to_anchor=(0.995, 0.995))
for a in axs:
    for tl in a.get_xticklabels() + a.get_yticklabels(): tl.set_fontsize(12); tl.set_fontweight('bold')
    a.grid(axis='y', color='#e5e5e5', lw=0.8); a.set_axisbelow(True)
    for s_ in ('top', 'right'): a.spines[s_].set_visible(False)
fig.tight_layout(w_pad=2.5, rect=[0, 0, 1, 0.92])
fig.savefig('Fig9_benchmark.png', dpi=300); fig.savefig('Fig9_benchmark.pdf')
print('ok')
