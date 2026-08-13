# %% [markdown]
# # CY Landscape Explorer
# 
# Exploration du paysage des variétés de Calabi-Yau pour la théorie des cordes hétérotique.
# 
# **Pipeline :** Topologie CY → Fibrés vectoriels → Spectre de particules → Couplages de Yukawa

# %% [markdown]
# ## Configuration
# 
# Modifier `CICYLIST_PATH` pour pointer vers le fichier téléchargé depuis Oxford.
# Si le fichier n'est pas disponible, les 14 CICYs embarquées seront utilisées.

# %%
import numpy as np
import os, sys, json, time

# Ajouter le bon dossier au path
# Cas 1 : le notebook est DANS cy_landscape/ (cd .. necessaire)
# Cas 2 : le notebook est A COTE de cy_landscape/ (deja bon)
# Cas 3 : chemin manuel
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

# === PARAMETRES ===
CICYLIST_PATH = "cicylist.txt"  # Chemin vers le fichier Oxford (7890 CICYs)
MAX_PS = 6                      # Nombre max d'espaces projectifs (5=rapide, 6=standard, None=tout)
N_RANDOM = 150                  # Nombre de fibrés aléatoires par rang
SEED = 42                       # Graine aléatoire (reproductibilité)
# ==================

rng = np.random.RandomState(SEED)

# %% [markdown]
# ## Étage 1 — Chargement et validation des CICYs

# %%
from cy_landscape.core.intersection import (
    compute_intersection_numbers, compute_euler_from_intersection, compute_c2_tangent)
from cy_landscape.core.bundles import CICYGeometry

# Chargement
if os.path.exists(CICYLIST_PATH):
    from cy_landscape.data.parse_oxford import load_oxford_file
    all_entries = load_oxford_file(CICYLIST_PATH)
    print(f"✓ {len(all_entries)} CICYs chargées depuis {CICYLIST_PATH}")
else:
    from cy_landscape.data.oxford_cicys import get_all_oxford
    all_entries = get_all_oxford()
    print(f"⚠ Fichier Oxford non trouvé, {len(all_entries)} CICYs embarquées utilisées")
    print(f"  Télécharger depuis: https://www-thphys.physics.ox.ac.uk/projects/CalabiYau/cicylist/cicylist.txt")

# Filtrage
if MAX_PS:
    all_entries = [e for e in all_entries if len(e['ambient']) <= MAX_PS]
    print(f"  Filtre max_ps ≤ {MAX_PS} : {len(all_entries)} CICYs")

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
    except:
        pass

print(f"  {len(valid_cicys)} CICYs validées (intersection numbers + χ)")

# %% [markdown]
# ### Distribution des caractéristiques d'Euler

# %%
from collections import Counter

chi_dist = Counter(c['chi'] for c, _ in valid_cicys)
print(f"Distribution de χ (Euler) :")
print(f"{'χ':>6} {'|h11-h21|':>10} {'count':>6}")
print(f"{'─'*6} {'─'*10} {'─'*6}")
for chi in sorted(chi_dist.keys(), key=abs)[:20]:
    dh = abs(chi) // 2
    print(f"{chi:>6} {dh:>10} {chi_dist[chi]:>6}")

has_chi6 = -6 in chi_dist or 6 in chi_dist
print(f"\nχ = ±6 (3 générations topologiques) : {'PRESENT' if has_chi6 else 'ABSENT'}")
if not has_chi6:
    print("→ Les 3 générations du Modèle Standard DOIVENT venir du fibré vectoriel")

# %% [markdown]
# ## Étage 2 — Construction de fibrés vectoriels
# 
# Trois méthodes :
# - **Monades positives** : 0 → V → B → C → 0 avec c_j − b_i > 0
# - **Extensions** : 0 → F₁ → V → F₂ → 0
# - **Monades classiques** : 0 → V → B → C → 0 (aléatoires)

# %%
from cy_landscape.core.monads import (
    MonadBundle, compute_monad_cohomology, generate_monads, check_map_exists)
from cy_landscape.core.positive_monads import generate_positive_monads
from cy_landscape.core.extensions import (
    ExtensionBundle, check_extension_exists, compute_extension_cohomology, generate_extensions)
from cy_landscape.core.bundles import GAUGE_GROUP_TABLE

def scan_all_methods(c, geom, rng, n_random):
    """Scan les 3 méthodes sur une CICY."""
    m = len(geom.ambient_dims)
    candidates = []

    for rank_V in [3, 4, 5]:
        # Monades classiques
        for monad in generate_monads(m, rank_V, max_charge=3, n_random=n_random, rng=rng):
            if not monad.c1_vanishes: continue
            map_ok, _ = check_map_exists(monad, geom.ambient_dims, geom.config_matrix)
            if not map_ok: continue
            cohom = compute_monad_cohomology(monad, geom.ambient_dims, geom.config_matrix)
            if cohom and abs(cohom[1] - cohom[2]) == 3:
                candidates.append(('monad', rank_V, monad, cohom))

        # Monades positives
        for monad in generate_positive_monads(m, rank_V, max_charge=4,
                                               n_systematic=n_random, rng=rng):
            map_ok, _ = check_map_exists(monad, geom.ambient_dims, geom.config_matrix)
            if not map_ok: continue
            cohom = compute_monad_cohomology(monad, geom.ambient_dims, geom.config_matrix)
            if cohom and abs(cohom[1] - cohom[2]) == 3:
                candidates.append(('pos_monad', rank_V, monad, cohom))

        # Extensions
        for ext in generate_extensions(m, rank_V, max_charge=3, n_random=n_random, rng=rng):
            exists, _ = check_extension_exists(ext, geom.ambient_dims, geom.config_matrix)
            if not exists: continue
            cohom = compute_extension_cohomology(ext, geom.ambient_dims, geom.config_matrix)
            if cohom and abs(cohom[1] - cohom[2]) == 3:
                pseudo = MonadBundle(ext.f1_charges + ext.f2_charges, ext.f2_charges)
                candidates.append(('extension', rank_V, pseudo, cohom))

    return candidates

# Lancement du scan
print("Scan en cours...")
t0 = time.time()
all_candidates = []

for i, (c, geom) in enumerate(valid_cicys):
    cands = scan_all_methods(c, geom, rng, N_RANDOM)
    if cands:
        all_candidates.extend([(c, geom, *cand) for cand in cands])
        print(f"  CICY #{c['num']:>4} ({c['h11']},{c['h21']}): {len(cands)} candidats à 3 gen")
    if (i + 1) % 50 == 0:
        print(f"    ... {i+1}/{len(valid_cicys)} ({time.time()-t0:.0f}s)")

print(f"\nTotal : {len(all_candidates)} fibrés à 3 générations ({time.time()-t0:.0f}s)")

# %% [markdown]
# ## Étage 3 — Test de stabilité (Hoppe complet)

# %%
from cy_landscape.core.stability_full import hoppe_criterion
from cy_landscape.core.monad_wedge import cohomology_wedge2_V, cohomology_end_V_approx
from cy_landscape.core.cohomology import extract_spectrum_su5, extract_spectrum_so10, extract_spectrum_e6

def extract_spectrum(gauge, cohom, monad, ambient, config):
    """Calcule le spectre complet."""
    if monad.rank_C == 1:
        try:
            w2 = cohomology_wedge2_V(ambient, config, monad.b_charges, monad.c_charges)
            w2V = w2['wedge2V']
        except:
            rV = monad.rank_V
            w2V = {0:0, 1: max(0, cohom[1]*(rV-1)//2), 2: max(0, cohom[2]*(rV-1)//2), 3:0}
    else:
        rV = monad.rank_V
        w2V = {0:0, 1: max(0, cohom[1]*(rV-1)//2), 2: max(0, cohom[2]*(rV-1)//2), 3:0}

    end_V = {0:1, 1: max(1, monad.rank_V**2-1), 2: max(1, monad.rank_V**2-1), 3:1}
    fmt = {"V": {i: cohom[i] for i in range(4)},
           "V_dual": {i: cohom[3-i] for i in range(4)},
           "wedge2V": w2V, "end_V": end_V}
    if gauge == "SU(5)": return extract_spectrum_su5(fmt)
    elif gauge == "SO(10)": return extract_spectrum_so10(fmt)
    return extract_spectrum_e6(fmt)

# Test de Hoppe + spectre
print("Test de Hoppe (rangs 1 à rk-1) en cours...")
t0 = time.time()
stable_results = []

for c, geom, kind, rank_V, monad, cohom in all_candidates:
    hoppe = hoppe_criterion(geom.ambient_dims, geom.config_matrix, monad, max_H=1)
    gauge = GAUGE_GROUP_TABLE.get(rank_V, {}).get("group", "?")
    sp = extract_spectrum(gauge, cohom, monad, geom.ambient_dims, geom.config_matrix)

    score = sp.sm_compatibility * 0.5 + 30
    if hoppe.get('stable'): score += 25
    if kind == 'pos_monad': score += 3
    if kind == 'extension': score += 2

    result = {
        'type': kind, 'stable': hoppe.get('stable', False),
        'cicy': c['num'], 'h11': c['h11'], 'h21': c['h21'], 'chi': c['chi'],
        'ambient': "×".join(f"P{n}" for n in c['ambient']),
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
        stable_results.append(result)

stable_results.sort(key=lambda r: r['score'], reverse=True)
n_stable = len(stable_results)
print(f"\n{'─'*60}")
print(f"  Candidats à 3 générations  : {len(all_candidates)}")
print(f"  HOPPE-STABLES              : {n_stable}")
print(f"  Temps                      : {time.time()-t0:.0f}s")
print(f"{'─'*60}")

# %% [markdown]
# ### Tableau des meilleurs candidats stables

# %%
if stable_results:
    print(f"{'#':>2} {'Type':<12} {'CICY':>5} {'Hodge':>8} {'Rk':>2} {'Jauge':>7} "
          f"{'H^i(V)':<16} {'Hig':>3} {'Exo':>3} {'Score':>5}")
    print(f"{'─'*2} {'─'*12} {'─'*5} {'─'*8} {'─'*2} {'─'*7} {'─'*16} {'─'*3} {'─'*3} {'─'*5}")
    for i, r in enumerate(stable_results[:20]):
        hi = r['cohomology']
        print(f"{i+1:>2} {r['type']:<12} #{r['cicy']:>4} ({r['h11']:>2},{r['h21']:>2}) "
              f"{r['rank_V']:>2} {r['gauge']:>7} [{hi[0]},{hi[1]},{hi[2]},{hi[3]}]"
              f"{'':>4} {r['higgs']:>3} {r['exotics']:>3} {r['score']:>5}")
else:
    print("Aucun fibré Hoppe-stable trouvé. Essayer avec --n-random plus élevé.")

# %% [markdown]
# ## Étage 4 — Couplages de Yukawa

# %%
from cy_landscape.core.yukawa import compute_yukawa_texture

if stable_results:
    print("Analyse des couplages de Yukawa :\n")
    cicys_dict = {c['num']: c for c, _ in valid_cicys}

    for i, r in enumerate(stable_results[:5]):
        c = cicys_dict.get(r['cicy'])
        if c is None: continue
        monad = MonadBundle(r['b_charges'], r['c_charges'])
        yukawa = compute_yukawa_texture(
            c['ambient'], c['config'], monad, r['gauge'])

        h1, h2 = r['cohomology'][1], r['cohomology'][2]
        gen_source = "H¹(V)" if h1 >= h2 else "H²(V)=H¹(V*)"

        print(f"  #{i+1} CICY #{r['cicy']} {r['type']} {r['gauge']}")
        print(f"     Générations dans {gen_source}")
        if yukawa.total_rank > 0:
            print(f"     Yukawa rang {yukawa.total_rank}/{yukawa.n_gen} → "
                  f"{yukawa.total_rank} massive(s)")
            if yukawa.eigenvalues is not None:
                print(f"     Masses normalisées : "
                      f"{', '.join(f'{ev:.4f}' for ev in yukawa.eigenvalues[:3])}")
        else:
            print(f"     Yukawa tree-level = 0 → masses via instantons")
        print()

# %% [markdown]
# ## Vérification de la map générique

# %%
from cy_landscape.core.explicit_map import verify_generic_map

if stable_results:
    print("Vérification de la map générique (rang matriciel) :\n")
    n_ok = 0
    for r in stable_results[:10]:
        if r['type'] == 'extension': continue
        c = cicys_dict.get(r['cicy'])
        if c is None: continue

        result = verify_generic_map(
            c['ambient'], c['config'],
            r['b_charges'], r['c_charges'],
            n_trials=5, rng=rng)

        status = "✓" if result['verified'] else "✗"
        if result['verified']: n_ok += 1
        print(f"  {status} CICY #{r['cicy']} {r['type']} {r['gauge']} : "
              f"rang {result.get('rank_actual_max','?')}/{result.get('rank_expected','?')}")

    print(f"\n  {n_ok} vérifiées sur {min(len(stable_results), 10)} testées")

# %% [markdown]
# ## Export des résultats

# %%
output_dir = "output_notebook"
os.makedirs(output_dir, exist_ok=True)

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.integer,)): return int(obj)
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        return super().default(obj)

export = {
    'n_cicys': len(valid_cicys),
    'n_candidates': len(all_candidates),
    'n_stable': len(stable_results),
    'parameters': {'max_ps': MAX_PS, 'n_random': N_RANDOM, 'seed': SEED},
    'results': stable_results[:100],
}

path = os.path.join(output_dir, 'results.json')
with open(path, 'w') as f:
    json.dump(export, f, indent=2, cls=NpEncoder)

print(f"Résultats exportés dans {path}")
print(f"  {export['n_cicys']} CICYs, {export['n_stable']} fibrés Hoppe-stables")

# %% [markdown]
# ## Étage 5 — Visualisation des résultats
# 
# Ces graphiques prennent tout leur sens sur un scan massif (7 890 CICYs).

# %%
import matplotlib
# matplotlib.use('Agg')  # Décommenter uniquement hors Jupyter
import matplotlib.pyplot as plt
from collections import Counter, defaultdict

# Préparer les données
all_results_list = stable_results  # Résultats Hoppe-stables

# === 1. DIAGRAMME DE HODGE ===
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# 1a. Toutes les CICYs validées
ax = axes[0]
h11_all = [c['h11'] for c, _ in valid_cicys]
h21_all = [c['h21'] for c, _ in valid_cicys]
ax.scatter(h11_all, h21_all, c='lightgrey', s=15, alpha=0.6, label='Toutes CICYs')

# CICYs avec fibrés stables
if all_results_list:
    stable_cicys = defaultdict(int)
    for r in all_results_list:
        stable_cicys[(r['h11'], r['h21'])] += 1

    h11_s = [k[0] for k in stable_cicys]
    h21_s = [k[1] for k in stable_cicys]
    counts = [stable_cicys[k] for k in stable_cicys]
    scatter = ax.scatter(h11_s, h21_s, c=counts, cmap='hot_r', s=80,
                         edgecolors='black', linewidth=0.5, zorder=5)
    plt.colorbar(scatter, ax=ax, label='Fibrés stables')

ax.plot([0, max(h11_all+h21_all)], [0, max(h11_all+h21_all)],
        'k--', alpha=0.2, label='h¹¹ = h²¹ (miroir)')
ax.set_xlabel('h¹¹'); ax.set_ylabel('h²¹')
ax.set_title('Diagramme de Hodge')
ax.legend(fontsize=8)

# 1b. Zoom sur la zone fertile
ax = axes[1]
if all_results_list:
    for r in all_results_list:
        color = {'pos_monad': '#e74c3c', 'extension': '#2ecc71', 'monad': '#3498db'}
        ax.scatter(r['h11'], r['h21'], c=color.get(r['type'], 'grey'),
                   s=r['score']*2, alpha=0.7, edgecolors='black', linewidth=0.3)
    # Légende manuelle
    for t, c, l in [('pos_monad','#e74c3c','Monade positive'),
                     ('extension','#2ecc71','Extension'),
                     ('monad','#3498db','Monade classique')]:
        if any(r['type']==t for r in all_results_list):
            ax.scatter([], [], c=c, s=60, label=l, edgecolors='black', linewidth=0.3)
    ax.legend(fontsize=8)
ax.set_xlabel('h¹¹'); ax.set_ylabel('h²¹')
ax.set_title('Zone fertile (taille ∝ score)')

plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'hodge_diagram.png'), dpi=150)
plt.show()
print("→ hodge_diagram.png")

# %% [markdown]
# ### Distribution des scores et comparaison des méthodes

# %%
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# 2. DISTRIBUTION DES SCORES
ax = axes[0]
if all_results_list:
    scores = [r['score'] for r in all_results_list]
    ax.hist(scores, bins=20, color='#3498db', edgecolor='white', alpha=0.8)
    ax.axvline(np.mean(scores), color='red', linestyle='--', label=f'Moyenne : {np.mean(scores):.1f}')
    ax.legend(fontsize=8)
ax.set_xlabel('Score'); ax.set_ylabel('Nombre de fibrés')
ax.set_title('Distribution des scores')

# 3. COMPARAISON DES METHODES
ax = axes[1]
if all_results_list:
    method_counts = Counter(r['type'] for r in all_results_list)
    labels = {'pos_monad': 'Monades\npositives', 'extension': 'Extensions', 'monad': 'Monades\nclassiques'}
    colors = {'pos_monad': '#e74c3c', 'extension': '#2ecc71', 'monad': '#3498db'}
    names = [labels.get(k, k) for k in method_counts]
    vals = list(method_counts.values())
    cols = [colors.get(k, 'grey') for k in method_counts]
    bars = ax.bar(names, vals, color=cols, edgecolor='white')
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                str(val), ha='center', fontsize=11, fontweight='bold')
ax.set_ylabel('Fibrés Hoppe-stables')
ax.set_title('Par méthode de construction')

# 4. GROUPES DE JAUGE
ax = axes[2]
if all_results_list:
    gauge_counts = Counter(r['gauge'] for r in all_results_list)
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
if all_results_list:
    top_n = min(5, len(all_results_list))
    fig, axes = plt.subplots(1, top_n, figsize=(4*top_n, 5))
    if top_n == 1: axes = [axes]

    for idx, (ax, r) in enumerate(zip(axes, all_results_list[:top_n])):
        ax.set_xlim(0, 10); ax.set_ylim(0, 10)
        ax.set_aspect('equal')
        ax.axis('off')

        # Fond coloré selon score
        score_frac = min(r['score'] / 110, 1.0)
        bg_color = plt.cm.RdYlGn(score_frac)
        ax.add_patch(plt.Rectangle((0, 0), 10, 10, facecolor=bg_color,
                                    alpha=0.15, edgecolor='black', linewidth=2))

        # Titre
        ax.text(5, 9.2, f"#{idx+1}", ha='center', fontsize=16, fontweight='bold')
        ax.text(5, 8.3, f"CICY #{r['cicy']}", ha='center', fontsize=12)
        ax.text(5, 7.5, f"{r['ambient']}", ha='center', fontsize=9, color='grey')

        # Infos
        lines = [
            f"h = ({r['h11']}, {r['h21']})   χ = {r['chi']}",
            f"",
            f"{r['type']}  rk {r['rank_V']}",
            f"Jauge : {r['gauge']}",
            f"",
            f"Générations : {r['n_gen']}",
            f"Higgs : {r['higgs']}    Exo : {r['exotics']}",
            f"H^i = {r['cohomology']}",
        ]
        for i, line in enumerate(lines):
            ax.text(5, 6.5 - i*0.7, line, ha='center', fontsize=9,
                    family='monospace')

        # Score
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
if all_results_list:
    best = all_results_list[0]
    reps = best.get('reps', {})

    if reps:
        fig, ax = plt.subplots(figsize=(8, 5))

        # Barres pour chaque représentation
        rep_names = list(reps.keys())
        rep_vals = [reps[k] for k in rep_names]

        # Couleurs par type physique
        colors_rep = []
        for name in rep_names:
            if '16' in name or '27' in name:
                colors_rep.append('#e74c3c')   # Fermions (rouge)
            elif '10' in name or '5' in name:
                colors_rep.append('#2ecc71')   # Higgs (vert)
            else:
                colors_rep.append('#95a5a6')   # Singlets (gris)

        bars = ax.bar(rep_names, rep_vals, color=colors_rep, edgecolor='white', width=0.6)
        for bar, val in zip(bars, rep_vals):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                    str(val), ha='center', fontsize=12, fontweight='bold')

        ax.set_ylabel('Multiplicité')
        ax.set_title(f"Spectre du candidat #1 — CICY #{best['cicy']} — {best['gauge']}")

        # Légende
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
print("═"*60)
print("  RÉSUMÉ DU SCAN")
print("═"*60)
print(f"""
  CICYs analysées       : {len(valid_cicys)}
  Fibrés à 3 générations: {len(all_candidates)}
  Hoppe-stables         : {len(stable_results)}

  Par méthode :""")
for t in ['pos_monad', 'extension', 'monad']:
    n = sum(1 for r in stable_results if r['type'] == t)
    if n: print(f"    {t:<15}: {n}")
print(f"""
  Par groupe de jauge :""")
for g in ['SO(10)', 'SU(5)', 'E₆']:
    n = sum(1 for r in stable_results if r['gauge'] == g)
    if n: print(f"    {g:<15}: {n}")

if stable_results:
    b = stable_results[0]
    print(f"""
  Meilleur candidat :
    CICY #{b['cicy']} ({b['ambient']})
    {b['type']} {b['gauge']} rk {b['rank_V']}
    {b['n_gen']} générations, {b['higgs']} Higgs, {b['exotics']} exotiques
    Score : {b['score']}
""")

print(f"  Images : {output_dir}/")
print(f"  Données : {output_dir}/results.json")
print("═"*60)
