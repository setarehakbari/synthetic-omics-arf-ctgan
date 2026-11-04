# Fig. 8: dilution summaries across repeated synthetic draws and controls.
# Same style as the other multi-panel figures of the paper: Times New Roman, 12 pt, bold,
# drawn on a 13-inch-wide canvas (the figure is scaled to the page width in the manuscript,
# exactly as for the R figures such as the utility and RF-importance figures).
# Usage: python make_fig8_repeated_draws.py per_draw_all_metrics.csv
import sys
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman', 'Liberation Serif', 'Tinos', 'DejaVu Serif']
rcParams['font.size'] = 12
rcParams['font.weight'] = 'bold'
rcParams['axes.labelweight'] = 'bold'
rcParams['axes.titleweight'] = 'bold'
rcParams['pdf.fonttype'] = 42

d = pd.read_csv(sys.argv[1]).rename(columns={'rbo_p0.95': 'rbo95'})
COL = {'ARF': '#F8766D', 'CTGAN': '#00BFC4', 'TVAE': '#E5A500', 'CTGAN-long': '#C77CFF',
       'Control A': '#BDBDBD', 'Control B': '#7F7F7F'}
groups = [('ARF', 'Observed', 'ARF'), ('CTGAN', 'Observed', 'CTGAN'), ('TVAE', 'Observed', 'TVAE'),
          ('CTGAN-long', 'Observed', 'CTGAN_long'), ('Control A', 'Control A', None), ('Control B', 'Control B', None)]
xlabels = ['ARF', 'CTGAN', 'TVAE', 'CTGAN-\nlong', 'Control\nA', 'Control\nB']
panels = [('mean_abs_shift', 'Mean absolute rank shift (Top-50)', '(a)'),
          ('dropout_pct', 'Top-50 dropout (%)', '(b)'),
          ('rbo95', 'RBO (p = 0.95)', '(c)'),
          ('jaccard_top50', 'Top-50 Jaccard', '(d)')]
fig, axs = plt.subplots(2, 2, figsize=(13, 10))
rng = np.random.default_rng(1)
for ax, (m, ylab, tag) in zip(axs.ravel(), panels):
    for i, (lab, cond, gen) in enumerate(groups):
        sub = d[(d.condition == cond) & (d.manuscript_draw == False)]
        if gen: sub = sub[sub.generator == gen]
        y = sub[m].values; c = COL[lab]
        # whiskers span the full range (minimum to maximum)
        ax.boxplot([y], positions=[i], widths=0.62, whis=(0, 100), showfliers=False, patch_artist=True,
                   boxprops=dict(facecolor=c, edgecolor='black', linewidth=0.8, alpha=0.85),
                   medianprops=dict(color='black', linewidth=1.6),
                   whiskerprops=dict(color='black', linewidth=0.8), capprops=dict(color='black', linewidth=0.8))
        if gen:  # individual draws are shown for the generators only; controls are summarized by the box
            ax.scatter(i + rng.uniform(-0.17, 0.17, len(y)), y, s=14, color='black', alpha=0.55, linewidths=0, zorder=3)
            ms = d[(d.condition == 'Observed') & (d.generator == gen) & (d.manuscript_draw == True)]
            if len(ms):
                ax.scatter([i], ms[m].values, marker='D', s=95, facecolor='white', edgecolor='black', linewidths=1.6, zorder=5)
    ax.set_xticks(range(6)); ax.set_xticklabels(xlabels, fontsize=12, fontweight='bold')
    ax.set_xlim(-0.6, 5.6)
    ax.set_ylabel(ylab, fontsize=12, fontweight='bold')
    ax.set_title(tag, fontsize=12, fontweight='bold')
    for tl in ax.get_yticklabels(): tl.set_fontsize(12); tl.set_fontweight('bold')
    ax.grid(axis='y', color='#e5e5e5', linewidth=0.8); ax.set_axisbelow(True)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
# boxed colour legend at the top right, as in the other figures of the paper
handles = [Patch(facecolor=COL[k], edgecolor='black', linewidth=0.8, label=k)
           for k in ['ARF', 'CTGAN', 'TVAE', 'CTGAN-long', 'Control A', 'Control B']]
handles.append(Line2D([0], [0], marker='D', ls='', markersize=8, markerfacecolor='white', markeredgecolor='black',
                      markeredgewidth=1.6, label='Main-analysis dataset'))
fig.legend(handles=handles, loc='upper right', ncol=7, frameon=True, facecolor='white', edgecolor='black',
           fancybox=False, prop={'weight': 'bold', 'size': 12}, handlelength=1.2, handletextpad=0.5, columnspacing=1.2,
           bbox_to_anchor=(0.995, 0.995))
fig.tight_layout(h_pad=2.0, w_pad=2.5, rect=[0, 0, 1, 0.95])
fig.savefig('Fig8_repeated_draws.png', dpi=300); fig.savefig('Fig8_repeated_draws.pdf')
print('ok')
