#!/usr/bin/env python3
"""Plot the theorem's coefficient and verified analytic core convergence.

This visualises formulas and the recorded core calculations. It does not
claim that large completed frames were materialised or factorised.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from verify import sharp_constant

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'figures'
DEST.mkdir(exist_ok=True)
record = json.loads((ROOT / 'results' / 'verification.json').read_text())
palette = ['#75436C', '#AD702A', '#327B73']
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 10,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.labelcolor': '#302D35', 'text.color': '#302D35',
    'xtick.color': '#49434D', 'ytick.color': '#49434D',
    'axes.edgecolor': '#A9A1AB', 'savefig.facecolor': '#FFFDFA',
    'figure.facecolor': '#FFFDFA', 'axes.facecolor': '#FFFDFA',
})

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
h_values = np.geomspace(1, 16, 400)
c_values = np.array([sharp_constant(float(h)) for h in h_values])
axes[0].plot(h_values, c_values, color=palette[0], lw=2.5)
points = np.array([1., 2., 4., 8., 16.])
coeffs = [sharp_constant(float(h)) for h in points]
axes[0].scatter(points, coeffs, color=palette[0], s=28, zorder=5)
for h, c in zip(points[:3], coeffs[:3]):
    axes[0].annotate(f'{c:.3f}', (h, c), xytext=(5, 9),
                     textcoords='offset points', fontsize=9)
axes[0].set_xscale('log', base=2)
axes[0].set_xticks(points, [str(int(h)) for h in points])
axes[0].set_xlabel('Coherence cap H')
axes[0].set_ylabel(r'Sharp coefficient $C(H)$')
axes[0].set_ylim(1.92, max(coeffs)+.12)
axes[0].set_title('A  Exact leading coefficient', loc='left', weight='bold', pad=13)
axes[0].grid(axis='y', color='#E8E0E7', linewidth=.7)

curves = []
for label, colour in zip(['1.01', '1.5', '8.0'], palette):
    rows = [r for r in record['large_rank_core_asymptotics']['records']
            if r['effective_profile_coherence'] == label]
    ranks = np.array([r['rank'] for r in rows])
    values = np.array([float(r['log_dk_inverse_entry_over_sqrt_rank']) for r in rows])
    limit = float(rows[0]['limit_constant_80_digits'])
    axes[1].plot(ranks, values, '-o', color=colour, lw=1.8, ms=4,
                 label=rf'$\bar H={float(label):g}$; limit {limit:.3f}')
    axes[1].axhline(limit, color=colour, ls='--', lw=1.1, alpha=.7)
    curves.append({'profile_parameter':label,'ranks':ranks.tolist(),
                   'normalised_core_log_product':values.tolist(),'limit':limit})
axes[1].set_xscale('log')
axes[1].set_xlabel('Rank k')
axes[1].set_ylabel(r'$\log[d_k |(A^{-1})_{1k}|]/\sqrt{k}$')
axes[1].set_title('B  Analytic core approaches the limit', loc='left', weight='bold', pad=13)
axes[1].grid(axis='y', color='#E8E0E7', linewidth=.7)
axes[1].legend(frameon=False, loc='lower right', fontsize=8.8)
fig.subplots_adjust(left=.07, right=.98, top=.90, bottom=.12, wspace=.32)
fig.savefig(DEST/'sharp-growth.png', dpi=200)
fig.savefig(DEST/'sharp-growth.pdf')
(DEST/'plot-data.json').write_text(json.dumps({
    'coefficient_curve':{'H':h_values.tolist(),'C_H':c_values.tolist()},
    'selected_coefficients':dict(zip(map(str, points), coeffs)),
    'core_curves':curves,
    'scope':'Analytic formulas and recorded triangular-core calculations; no large full-frame experiment.',
}, indent=2)+'\n')
print('Created figures/sharp-growth.png, sharp-growth.pdf and plot-data.json')
