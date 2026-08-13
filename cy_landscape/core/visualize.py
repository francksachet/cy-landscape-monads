"""
visualize.py -- Generation des graphiques de synthese.

Produit 4 PNG a partir du fichier results.json :
  - hodge_diagram.png    : diagramme de Hodge (h11 vs h21)
  - statistics.png       : scores, methodes, groupes de jauge
  - top_candidates.png   : fiches des 5 meilleurs
  - best_spectrum.png    : spectre du meilleur candidat
"""
import os, json
import numpy as np
from collections import Counter, defaultdict


def generate_all_plots(results_path, output_dir=None):
    """Genere tous les graphiques depuis le fichier results.json."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch
    except ImportError:
        print("  matplotlib non installe, graphiques ignores")
        print("  pip install matplotlib")
        return

    with open(results_path) as f:
        data = json.load(f)

    results = data.get('results', [])
    if not results:
        print("  Aucun resultat a visualiser")
        return

    if output_dir is None:
        output_dir = os.path.dirname(results_path)

    stable = [r for r in results if r.get('stable', False)]
    if not stable:
        stable = results  # Fallback

    # ─── 1. DIAGRAMME DE HODGE ───
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    cicys_with_stable = defaultdict(int)
    all_h = set()
    for r in stable:
        cicys_with_stable[(r['h11'], r['h21'])] += 1
        all_h.add((r['h11'], r['h21']))

    h11_s = [k[0] for k in cicys_with_stable]
    h21_s = [k[1] for k in cicys_with_stable]
    counts = [cicys_with_stable[k] for k in cicys_with_stable]

    scatter = ax.scatter(h11_s, h21_s, c=counts, cmap='hot_r', s=80,
                         edgecolors='black', linewidth=0.5, zorder=5)
    plt.colorbar(scatter, ax=ax, label='Fibres stables')
    mx = max(max(h11_s, default=1), max(h21_s, default=1)) + 2
    ax.plot([0, mx], [0, mx], 'k--', alpha=0.2)
    ax.set_xlabel('h¹¹'); ax.set_ylabel('h²¹')
    ax.set_title('Diagramme de Hodge — CICYs fertiles')

    ax = axes[1]
    colors_type = {'pos_monad': '#e74c3c', 'extension': '#2ecc71', 'monad': '#3498db'}
    for r in stable:
        ax.scatter(r['h11'], r['h21'],
                   c=colors_type.get(r['type'], 'grey'),
                   s=max(20, r['score'] * 2), alpha=0.7,
                   edgecolors='black', linewidth=0.3)
    for t, c, l in [('pos_monad','#e74c3c','Monade positive'),
                     ('extension','#2ecc71','Extension'),
                     ('monad','#3498db','Monade classique')]:
        if any(r['type'] == t for r in stable):
            ax.scatter([], [], c=c, s=60, label=l, edgecolors='black', linewidth=0.3)
    ax.legend(fontsize=8)
    ax.set_xlabel('h¹¹'); ax.set_ylabel('h²¹')
    ax.set_title('Par methode (taille ~ score)')

    plt.tight_layout()
    path = os.path.join(output_dir, 'hodge_diagram.png')
    plt.savefig(path, dpi=150); plt.close()
    print(f"  → {path}")

    # ─── 2. STATISTIQUES ───
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    ax = axes[0]
    scores = [r['score'] for r in stable]
    ax.hist(scores, bins=20, color='#3498db', edgecolor='white', alpha=0.8)
    ax.axvline(np.mean(scores), color='red', linestyle='--',
               label=f'Moyenne : {np.mean(scores):.1f}')
    ax.legend(fontsize=8)
    ax.set_xlabel('Score'); ax.set_ylabel('Nombre')
    ax.set_title('Distribution des scores')

    ax = axes[1]
    mc = Counter(r['type'] for r in stable)
    labels = {'pos_monad': 'Monades\npositives', 'extension': 'Extensions',
              'monad': 'Monades\nclassiques'}
    names = [labels.get(k, k) for k in mc]
    vals = list(mc.values())
    cols = [colors_type.get(k, 'grey') for k in mc]
    bars = ax.bar(names, vals, color=cols, edgecolor='white')
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                str(val), ha='center', fontsize=11, fontweight='bold')
    ax.set_ylabel('Fibres stables')
    ax.set_title('Par methode')

    ax = axes[2]
    gc = Counter(r['gauge'] for r in stable)
    colors_g = {'SO(10)': '#e67e22', 'SU(5)': '#9b59b6', 'E₆': '#1abc9c'}
    bars = ax.bar(list(gc.keys()), list(gc.values()),
                  color=[colors_g.get(k, 'grey') for k in gc], edgecolor='white')
    for bar, val in zip(bars, gc.values()):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                str(val), ha='center', fontsize=11, fontweight='bold')
    ax.set_ylabel('Fibres stables')
    ax.set_title('Par groupe de jauge')

    plt.tight_layout()
    path = os.path.join(output_dir, 'statistics.png')
    plt.savefig(path, dpi=150); plt.close()
    print(f"  → {path}")

    # ─── 3. FICHES TOP 5 ───
    top_n = min(5, len(stable))
    fig, axes = plt.subplots(1, top_n, figsize=(4 * top_n, 5))
    if top_n == 1: axes = [axes]

    for idx, (ax, r) in enumerate(zip(axes, stable[:top_n])):
        ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.set_aspect('equal'); ax.axis('off')
        sf = min(r['score'] / 110, 1.0)
        ax.add_patch(plt.Rectangle((0, 0), 10, 10, facecolor=plt.cm.RdYlGn(sf),
                                    alpha=0.15, edgecolor='black', linewidth=2))
        ax.text(5, 9.2, f"#{idx+1}", ha='center', fontsize=16, fontweight='bold')
        ax.text(5, 8.3, f"CICY #{r['cicy']}", ha='center', fontsize=12)
        ax.text(5, 7.5, r.get('ambient', ''), ha='center', fontsize=9, color='grey')

        hi = r.get('cohomology', [0,0,0,0])
        lines = [
            f"h = ({r['h11']}, {r['h21']})   chi = {r['chi']}",
            "", f"{r['type']}  rk {r['rank_V']}", f"Jauge : {r['gauge']}",
            "", f"Generations : {r['n_gen']}",
            f"Higgs : {r['higgs']}    Exo : {r['exotics']}",
            f"H^i = {hi}",
        ]
        for i, line in enumerate(lines):
            ax.text(5, 6.5 - i * 0.7, line, ha='center', fontsize=9, family='monospace')
        ax.text(5, 0.8, f"Score : {r['score']}", ha='center', fontsize=14,
                fontweight='bold', color='darkgreen' if r['score'] > 90 else 'darkorange')

    plt.tight_layout()
    path = os.path.join(output_dir, 'top_candidates.png')
    plt.savefig(path, dpi=150); plt.close()
    print(f"  → {path}")

    # ─── 4. SPECTRE DU MEILLEUR ───
    best = stable[0]
    reps = best.get('reps', {})
    if reps:
        fig, ax = plt.subplots(figsize=(8, 5))
        rep_names = list(reps.keys())
        rep_vals = [reps[k] for k in rep_names]
        colors_rep = []
        for name in rep_names:
            if '16' in name or '27' in name: colors_rep.append('#e74c3c')
            elif '10' in name or '5' in name: colors_rep.append('#2ecc71')
            else: colors_rep.append('#95a5a6')
        bars = ax.bar(rep_names, rep_vals, color=colors_rep, edgecolor='white', width=0.6)
        for bar, val in zip(bars, rep_vals):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                    str(val), ha='center', fontsize=12, fontweight='bold')
        ax.set_ylabel('Multiplicite')
        ax.set_title(f"Spectre #{best['cicy']} — {best['gauge']}")
        ax.legend(handles=[Patch(facecolor='#e74c3c', label='Fermions'),
                           Patch(facecolor='#2ecc71', label='Higgs'),
                           Patch(facecolor='#95a5a6', label='Singlets')], fontsize=9)
        plt.tight_layout()
        path = os.path.join(output_dir, 'best_spectrum.png')
        plt.savefig(path, dpi=150); plt.close()
        print(f"  → {path}")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        generate_all_plots(sys.argv[1])
    else:
        print("Usage: python -m cy_landscape.core.visualize results.json")
