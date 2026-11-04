# Fig. 4: domain-separator AUC (random-forest classifier), real = positive class.
# Style: Times New Roman, 12 pt, bold (same as the other figures in the paper).
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman', 'Liberation Serif', 'Tinos', 'DejaVu Serif']
rcParams['font.size'] = 12
rcParams['font.weight'] = 'bold'
rcParams['axes.labelweight'] = 'bold'
rcParams['axes.titleweight'] = 'bold'
rcParams['pdf.fonttype'] = 42

gens = ['ARF', 'CTGAN', 'TVAE']
mean = [0.701, 1.000, 1.000]
sd   = [0.0726, 0.0, 0.000156]
mn   = [0.583, 1.000, 0.999]
mx   = [0.821, 1.000, 1.000]
col  = ['#F8766D', '#00BFC4', '#E5A500']

# Figure size equals the printed size in the manuscript (about 6.3 x 5.3 in),
# so that 12 pt text is 12 pt on the page.
fig, ax = plt.subplots(figsize=(6.33, 5.35), dpi=300)
x = range(3)
ax.bar(x, mean, width=0.6, color=col, zorder=2)
for i in x:
    ax.vlines(i, mn[i], mx[i], color='#333333', lw=1.2, zorder=3)
    ax.hlines([mn[i], mx[i]], i - 0.07, i + 0.07, color='#333333', lw=1.2, zorder=3)
    ax.errorbar(i, mean[i], yerr=sd[i], color='black', lw=3, capsize=0, zorder=4)
lab = ['0.701\n(SD 0.073)', '1.000\n(SD 0)', '1.000\n(SD 0.0002)']
for i in x:
    top = max(mx[i], mean[i]) + 0.02
    ax.text(i, top if i == 0 else 1.02, lab[i], ha='center', va='bottom', fontsize=12, fontweight='bold')
ax.axhline(0.5, color='black', ls='--', lw=1.2, zorder=1)
ax.text(2.36, 0.52, 'chance', ha='left', va='bottom', fontsize=12, fontweight='bold')
ax.set_xlim(-0.5, 2.95)
ax.set_xticks(list(x)); ax.set_xticklabels(gens, fontsize=12, fontweight='bold')
ax.set_ylim(0, 1.2); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
for t in ax.get_yticklabels(): t.set_fontsize(12); t.set_fontweight('bold')
ax.set_ylabel('Domain-separator AUC (real = positive class)', fontsize=12, fontweight='bold')
ax.grid(axis='y', color='#e8e8e8', zorder=0)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig('Fig4_privacy.png', dpi=300)
fig.savefig('Fig4_privacy.pdf')
print('ok', matplotlib.font_manager.findfont(rcParams['font.family'][0]))
