# %% [markdown]
# # CY Landscape Explorer
#
# Exploration du paysage des variétés de Calabi-Yau pour la théorie des cordes hétérotique.
#
# **Deux modes d'utilisation** (à choisir dans la cellule de configuration ci-dessous) :
# - `MODE = "load"` : charge un `results.json` déjà produit (par `main_optimized.py` par exemple) et affiche directement les graphiques. Rapide, ne relance aucun calcul.
# - `MODE = "scan"` : relance le pipeline complet dans le notebook (étages 1 à 4), puis affiche les graphiques. Utile pour tester sur un petit échantillon ou explorer pas à pas.

# %% [markdown]
# ## Configuration
#
# C'est la seule cellule à modifier.

# %%
import os, sys

# === MODE D'UTILISATION ===
MODE = "load"   # "load" : charger un results.json existant (rapide)
                # "scan" : relancer le scan complet dans le notebook (lent)

# --- Utilise si MODE == "load" ---
RESULTS_JSON_PATH = "output_optimized/results.json"   # Chemin vers le fichier a charger

# --- Utilise si MODE == "scan" ---
CICYLIST_PATH = "cicylist.txt"  # Chemin vers le fichier Oxford (7890 CICYs)
MAX_PS = 5                      # Nombre max d'espaces projectifs (5=rapide, 6=standard)
N_RANDOM = 150                  # Nombre de fibres aleatoires par rang
SEED = 42                       # Graine aleatoire (reproductibilite)
# ===========================

# Ajouter le bon dossier au path
_found = False
for _p in ['.', '..', os.path.dirname(os.path.abspath(''))]:
    if os.path.isdir(os.path.join(_p, 'cy_landscape', 'core')):
        sys.path.insert(0, os.path.abspath(_p))
        _found = True
        break
if not _found:
    print("ERREUR : dossier cy_landscape/ introuvable.")
    print("Placer ce notebook A COTE du dossier cy_landscape/, pas dedans.")
    print("Structure attendue :")
    print("  mon_projet/")
    print("    cy_landscape/")
    print("    CY_Landscape_Explorer.ipynb  <-- ici")
    print("    cicylist.txt")

print(f"Mode selectionne : {MODE!r}")

# %% [markdown]
# ## Chargement des résultats
#
# Cette cellule s'adapte automatiquement au `MODE` choisi ci-dessus.
# Elle produit dans les deux cas une liste `all_results_list` de dicts,
# au même format, utilisée par toutes les cellules de visualisation qui suivent.

# %%
import json
import numpy as np

all_results_list = []
n_cicys_total = None
scan_params = {}

if MODE == "load":
    # --- Chargement direct d'un results.json existant ---
    if not os.path.exists(RESULTS_JSON_PATH):
        print(f"ERREUR : fichier introuvable : {RESULTS_JSON_PATH}")
        print("Verifie le chemin, ou lance d'abord :")
        print("  python -m cy_landscape.main_optimized cicylist.txt --max-ps 6 -j 8")
    else:
        with open(RESULTS_JSON_PATH) as f:
            data = json.load(f)
        all_results_list = data.get('results', [])
        n_cicys_total = data.get('n_cicys', None)
        scan_params = data.get('parameters', {})

        print(f"✓ Charge depuis {RESULTS_JSON_PATH}")
        print(f"  CICYs analysees      : {n_cicys_total}")
        print(f"  Fibres Hoppe-stables : {data.get('n_stable', len(all_results_list))}")
        if scan_params:
            print(f"  Parametres du scan   : {scan_params}")

elif MODE == "scan":
    # --- Relance du pipeline complet dans le notebook ---
    import time
    from cy_landscape.core.intersection import (
        compute_intersection_numbers, compute_euler_from_intersection, compute_c2_tangent)
    from cy_landscape.core.bundles import CICYGeometry, GAUGE_GROUP_TABLE
    from cy_landscape.core.monads import (
        MonadBundle, compute_monad_cohomology, generate_monads, check_map_exists)
    from cy_landscape.core.positive_monads import generate_positive_monads
    from cy_landscape.core.extensions import (
        ExtensionBundle, check_extension_exists, compute_extension_cohomology, generate_extensions)
    from cy_landscape.core.monad_wedge import cohomology_wedge2_V
    from cy_landscape.core.stability_full import hoppe_criterion
    from cy_landscape.core.cohomology import extract_spectrum_su5, extract_spectrum_so10, extract_spectrum_e6

    rng = np.random.RandomState(SEED)

    # Chargement des CICYs
    if os.path.exists(CICYLIST_PATH):
        from cy_landscape.data.parse_oxford import load_oxford_file
        all_entries = load_oxford_file(CICYLIST_PATH)
        print(f"✓ {len(all_entries)} CICYs chargees depuis {CICYLIST_PATH}")
    else:
        from cy_landscape.data.oxford_cicys import get_all_oxford
        all_entries = get_all_oxford()
        print(f"⚠ Fichier Oxford non trouve, {len(all_entries)} CICYs embarquees utilisees")

    if MAX_PS:
        all_entries = [e for e in all_entries if len(e['ambient']) <= MAX_PS]
        print(f"  Filtre max_ps <= {MAX_PS} : {len(all_entries)} CICYs")

    # Validation
    valid_cicys = []
    for c in all_entries:
        try:
            d = compute_intersection_numbers(c['ambient'], c['config'])
            chi = compute_euler_from_intersection(c['ambient'], c['config'], d)
            if chi != c['chi']:
                continue
            c2 = compute_c2_tangent(c['ambient'], c['config'], d)
            geom = CICYGeometry(
                ambient_dims=c['ambient'], config_matrix=c['config'],
                h11=c['h11'], h21=c['h21'],
                intersection_numbers=d, c2_tangent=c2)
            valid_cicys.append((c, geom))
        except Exception:
            pass
    print(f"  {len(valid_cicys)} CICYs validees")
    n_cicys_total = len(valid_cicys)

    def _spectrum(gauge, cohom_fmt):
        if gauge == "SU(5)": return extract_spectrum_su5(cohom_fmt)
        elif gauge == "SO(10)": return extract_spectrum_so10(cohom_fmt)
        else: return extract_spectrum_e6(cohom_fmt)

    def scan_all_methods(c, geom, rng, n_random):
        m = len(geom.ambient_dims)
        candidates = []
        for rank_V in [3, 4, 5]:
            for monad in generate_monads(m, rank_V, max_charge=3, n_random=n_random, rng=rng):
                if not monad.c1_vanishes: continue
                map_ok, _ = check_map_exists(monad, geom.ambient_dims, geom.config_matrix)
                if not map_ok: continue
                cohom = compute_monad_cohomology(monad, geom.ambient_dims, geom.config_matrix)
                if cohom and abs(cohom[1] - cohom[2]) == 3:
                    candidates.append(('monad', rank_V, monad, cohom))
            for monad in generate_positive_monads(m, rank_V, max_charge=4,
                                                   n_systematic=n_random, rng=rng):
                map_ok, _ = check_map_exists(monad, geom.ambient_dims, geom.config_matrix)
                if not map_ok: continue
                cohom = compute_monad_cohomology(monad, geom.ambient_dims, geom.config_matrix)
                if cohom and abs(cohom[1] - cohom[2]) == 3:
                    candidates.append(('pos_monad', rank_V, monad, cohom))
            for ext in generate_extensions(m, rank_V, max_charge=3, n_random=n_random, rng=rng):
                exists, _ = check_extension_exists(ext, geom.ambient_dims, geom.config_matrix)
                if not exists: continue
                cohom = compute_extension_cohomology(ext, geom.ambient_dims, geom.config_matrix)
                if cohom and abs(cohom[1] - cohom[2]) == 3:
                    pseudo = MonadBundle(ext.f1_charges + ext.f2_charges, ext.f2_charges)
                    candidates.append(('extension', rank_V, pseudo, cohom))
        return candidates

    print("Scan en cours...")
    t0 = time.time()
    all_candidates = []
    for i, (c, geom) in enumerate(valid_cicys):
        cands = scan_all_methods(c, geom, rng, N_RANDOM)
        if cands:
            all_candidates.extend([(c, geom, *cand) for cand in cands])
            print(f"  CICY #{c['num']:>4} ({c['h11']},{c['h21']}): {len(cands)} candidats a 3 gen")
        if (i + 1) % 50 == 0:
            print(f"    ... {i+1}/{len(valid_cicys)} ({time.time()-t0:.0f}s)")
    print(f"\nTotal : {len(all_candidates)} fibres a 3 generations ({time.time()-t0:.0f}s)")

    print("\nTest de Hoppe (rangs 1 a rk-1) en cours...")
    t0 = time.time()
    for c, geom, kind, rank_V, monad, cohom in all_candidates:
        hoppe = hoppe_criterion(geom.ambient_dims, geom.config_matrix, monad, max_H=1)
        gauge = GAUGE_GROUP_TABLE.get(rank_V, {}).get("group", "?")
        if monad.rank_C == 1:
            try:
                w2 = cohomology_wedge2_V(geom.ambient_dims, geom.config_matrix,
                                          monad.b_charges, monad.c_charges)
                w2V = w2['wedge2V']
            except Exception:
                rV = monad.rank_V
                w2V = {0:0, 1: max(0, cohom[1]*(rV-1)//2), 2: max(0, cohom[2]*(rV-1)//2), 3:0}
        else:
            rV = monad.rank_V
            w2V = {0:0, 1: max(0, cohom[1]*(rV-1)//2), 2: max(0, cohom[2]*(rV-1)//2), 3:0}
        end_V = {0:1, 1: max(1, monad.rank_V**2-1), 2: max(1, monad.rank_V**2-1), 3:1}
        fmt = {"V": {i: cohom[i] for i in range(4)},
               "V_dual": {i: cohom[3-i] for i in range(4)},
               "wedge2V": w2V, "end_V": end_V}
        sp = _spectrum(gauge, fmt)

        score = sp.sm_compatibility * 0.5 + 30
        if hoppe.get('stable'): score += 25
        if kind == 'pos_monad': score += 3
        if kind == 'extension': score += 2

        result = {
            'type': kind, 'stable': hoppe.get('stable', False),
            'cicy': c['num'], 'h11': c['h11'], 'h21': c['h21'], 'chi': c['chi'],
            'ambient': "x".join(f"P{n}" for n in c['ambient']),
            'rank_V': rank_V, 'gauge': gauge,
            'n_gen': sp.n_generations, 'higgs': sp.n_higgs_candidates,
            'exotics': sp.n_exotics, 'score': round(score, 1),
            'cohomology': [int(cohom[i]) for i in range(4)],
            'reps': {k: v for k, v in sp.representations.items() if v > 0},
            'b_charges': [list(b) for b in monad.b_charges],
            'c_charges': [list(cc) for cc in monad.c_charges],
            'hoppe': hoppe.get('reason', ''),
        }
        if result['stable']:
            all_results_list.append(result)

    all_results_list.sort(key=lambda r: r['score'], reverse=True)
    print(f"\nHoppe-stables : {len(all_results_list)}  ({time.time()-t0:.0f}s)")

else:
    print(f"MODE invalide : {MODE!r}. Utiliser 'load' ou 'scan'.")

print(f"\n{len(all_results_list)} resultats disponibles pour la visualisation.")

# %% [markdown]
# ## Export (uniquement si MODE == "scan")
#
# Sauvegarde les résultats fraîchement calculés, pour pouvoir les recharger
# plus tard en mode `"load"` sans tout relancer.

# %%
output_dir = "output_notebook"
os.makedirs(output_dir, exist_ok=True)

if MODE == "scan" and all_results_list:
    class NpEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.integer,)): return int(obj)
            if isinstance(obj, (np.floating,)): return float(obj)
            if isinstance(obj, np.ndarray): return obj.tolist()
            return super().default(obj)

    export = {
        'n_cicys': n_cicys_total,
        'n_stable': len(all_results_list),
        'parameters': {'max_ps': MAX_PS, 'n_random': N_RANDOM, 'seed': SEED},
        'results': all_results_list[:500],
    }
    path = os.path.join(output_dir, 'results.json')
    with open(path, 'w') as f:
        json.dump(export, f, indent=2, cls=NpEncoder)
    print(f"Resultats exportes dans {path}")
    print(f"Tu pourras les recharger plus tard avec MODE='load' et")
    print(f"RESULTS_JSON_PATH = {path!r}")
else:
    print("Rien a exporter (MODE='load' : les resultats viennent deja d'un fichier).")

# %% [markdown]
# ## Visualisation des résultats
#
# Ces cellules fonctionnent à l'identique, que les résultats viennent d'un
# scan fraîchement calculé ou d'un `results.json` chargé.

# %%
if not all_results_list:
    print("Aucun resultat a visualiser. Verifie MODE et le chemin des donnees ci-dessus.")
else:
    import matplotlib
    # matplotlib.use('Agg')  # Decommenter uniquement hors Jupyter
    import matplotlib.pyplot as plt
    from collections import Counter, defaultdict

    stable_results = [r for r in all_results_list if r.get('stable', True)]
    print(f"{len(stable_results)} fibres stables a visualiser")

# %% [markdown]
# ### Diagramme de Hodge

# %%
if stable_results:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    cicys_with_stable = defaultdict(int)
    for r in stable_results:
        cicys_with_stable[(r['h11'], r['h21'])] += 1

    h11_s = [k[0] for k in cicys_with_stable]
    h21_s = [k[1] for k in cicys_with_stable]
    counts = [cicys_with_stable[k] for k in cicys_with_stable]
    scatter = ax.scatter(h11_s, h21_s, c=counts, cmap='hot_r', s=80,
                         edgecolors='black', linewidth=0.5, zorder=5)
    plt.colorbar(scatter, ax=ax, label='Fibrés stables')
    mx = max(max(h11_s, default=1), max(h21_s, default=1)) + 2
    ax.plot([0, mx], [0, mx], 'k--', alpha=0.2)
    ax.set_xlabel('h¹¹'); ax.set_ylabel('h²¹')
    ax.set_title('Diagramme de Hodge')

    ax = axes[1]
    colors_type = {'pos_monad': '#e74c3c', 'extension': '#2ecc71', 'monad': '#3498db'}
    for r in stable_results:
        ax.scatter(r['h11'], r['h21'], c=colors_type.get(r['type'], 'grey'),
                   s=max(20, r['score']*2), alpha=0.7, edgecolors='black', linewidth=0.3)
    for t, c, l in [('pos_monad','#e74c3c','Monade positive'),
                     ('extension','#2ecc71','Extension'),
                     ('monad','#3498db','Monade classique')]:
        if any(r['type']==t for r in stable_results):
            ax.scatter([], [], c=c, s=60, label=l, edgecolors='black', linewidth=0.3)
    ax.legend(fontsize=8)
    ax.set_xlabel('h¹¹'); ax.set_ylabel('h²¹')
    ax.set_title('Zone fertile (taille ∝ score)')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'hodge_diagram.png'), dpi=150)
    plt.show()
    print("→ hodge_diagram.png")

# %% [markdown]
# ### Statistiques : scores, méthodes, groupes de jauge

# %%
if stable_results:
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    ax = axes[0]
    scores = [r['score'] for r in stable_results]
    ax.hist(scores, bins=20, color='#3498db', edgecolor='white', alpha=0.8)
    ax.axvline(np.mean(scores), color='red', linestyle='--', label=f'Moyenne : {np.mean(scores):.1f}')
    ax.legend(fontsize=8)
    ax.set_xlabel('Score'); ax.set_ylabel('Nombre de fibrés')
    ax.set_title('Distribution des scores')

    ax = axes[1]
    method_counts = Counter(r['type'] for r in stable_results)
    labels = {'pos_monad': 'Monades\npositives', 'extension': 'Extensions', 'monad': 'Monades\nclassiques'}
    names = [labels.get(k, k) for k in method_counts]
    vals = list(method_counts.values())
    cols = [colors_type.get(k, 'grey') for k in method_counts]
    bars = ax.bar(names, vals, color=cols, edgecolor='white')
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                str(val), ha='center', fontsize=11, fontweight='bold')
    ax.set_ylabel('Fibrés Hoppe-stables')
    ax.set_title('Par méthode de construction')

    ax = axes[2]
    gauge_counts = Counter(r['gauge'] for r in stable_results)
    colors_g = {'SO(10)': '#e67e22', 'SU(5)': '#9b59b6', 'E₆': '#1abc9c'}
    names_g = list(gauge_counts.keys())
    vals_g = list(gauge_counts.values())
    cols_g = [colors_g.get(k, 'grey') for k in names_g]
    bars = ax.bar(names_g, vals_g, color=cols_g, edgecolor='white')
    for bar, val in zip(bars, vals_g):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                str(val), ha='center', fontsize=11, fontweight='bold')
    ax.set_ylabel('Fibrés Hoppe-stables')
    ax.set_title('Par groupe de jauge GUT')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'statistics.png'), dpi=150)
    plt.show()
    print("→ statistics.png")

# %% [markdown]
# ### Fiches des meilleurs candidats

# %%
if stable_results:
    top_n = min(5, len(stable_results))
    fig, axes = plt.subplots(1, top_n, figsize=(4*top_n, 5))
    if top_n == 1: axes = [axes]

    for idx, (ax, r) in enumerate(zip(axes, stable_results[:top_n])):
        ax.set_xlim(0, 10); ax.set_ylim(0, 10)
        ax.set_aspect('equal')
        ax.axis('off')
        score_frac = min(r['score'] / 110, 1.0)
        bg_color = plt.cm.RdYlGn(score_frac)
        ax.add_patch(plt.Rectangle((0, 0), 10, 10, facecolor=bg_color,
                                    alpha=0.15, edgecolor='black', linewidth=2))
        ax.text(5, 9.2, f"#{idx+1}", ha='center', fontsize=16, fontweight='bold')
        ax.text(5, 8.3, f"CICY #{r['cicy']}", ha='center', fontsize=12)
        ax.text(5, 7.5, f"{r.get('ambient','')}", ha='center', fontsize=9, color='grey')

        hi = r.get('cohomology', [0,0,0,0])
        lines = [
            f"h = ({r['h11']}, {r['h21']})   χ = {r['chi']}",
            f"", f"{r['type']}  rk {r['rank_V']}", f"Jauge : {r['gauge']}",
            f"", f"Générations : {r['n_gen']}",
            f"Higgs : {r['higgs']}    Exo : {r['exotics']}",
            f"H^i = {hi}",
        ]
        for i, line in enumerate(lines):
            ax.text(5, 6.5 - i*0.7, line, ha='center', fontsize=9, family='monospace')

        ax.text(5, 0.8, f"Score : {r['score']}", ha='center',
                fontsize=14, fontweight='bold',
                color='darkgreen' if r['score'] > 90 else 'darkorange')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'top_candidates.png'), dpi=150)
    plt.show()
    print("→ top_candidates.png")

# %% [markdown]
# ### Spectre du meilleur candidat

# %%
if stable_results:
    best = stable_results[0]
    reps = best.get('reps', {})

    if reps:
        fig, ax = plt.subplots(figsize=(8, 5))
        rep_names = list(reps.keys())
        rep_vals = [reps[k] for k in rep_names]
        colors_rep = []
        for name in rep_names:
            if '16' in name or '27' in name:
                colors_rep.append('#e74c3c')
            elif '10' in name or '5' in name:
                colors_rep.append('#2ecc71')
            else:
                colors_rep.append('#95a5a6')

        bars = ax.bar(rep_names, rep_vals, color=colors_rep, edgecolor='white', width=0.6)
        for bar, val in zip(bars, rep_vals):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                    str(val), ha='center', fontsize=12, fontweight='bold')

        ax.set_ylabel('Multiplicité')
        ax.set_title(f"Spectre du candidat #1 — CICY #{best['cicy']} — {best['gauge']}")

        from matplotlib.patches import Patch
        legend = [Patch(facecolor='#e74c3c', label='Fermions (générations)'),
                  Patch(facecolor='#2ecc71', label='Higgs'),
                  Patch(facecolor='#95a5a6', label='Singlets (moduli)')]
        ax.legend(handles=legend, fontsize=9)

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'best_spectrum.png'), dpi=150)
        plt.show()
        print("→ best_spectrum.png")

# %% [markdown]
# ### Résumé final

# %%
if stable_results:
    print("═"*60)
    print("  RÉSUMÉ")
    print("═"*60)
    print(f"""
  Mode utilisé          : {MODE}
  CICYs analysées       : {n_cicys_total}
  Fibrés Hoppe-stables  : {len(stable_results)}

  Par méthode :""")
    for t in ['pos_monad', 'extension', 'monad']:
        n = sum(1 for r in stable_results if r['type'] == t)
        if n: print(f"    {t:<15}: {n}")
    print(f"""
  Par groupe de jauge :""")
    for g in ['SO(10)', 'SU(5)', 'E₆']:
        n = sum(1 for r in stable_results if r['gauge'] == g)
        if n: print(f"    {g:<15}: {n}")

    b = stable_results[0]
    print(f"""
  Meilleur candidat :
    CICY #{b['cicy']} ({b.get('ambient','')})
    {b['type']} {b['gauge']} rk {b['rank_V']}
    {b['n_gen']} générations, {b['higgs']} Higgs, {b['exotics']} exotiques
    Score : {b['score']}
""")
    print(f"  Images : {output_dir}/")
