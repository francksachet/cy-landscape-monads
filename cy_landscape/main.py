#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════╗
║          CY LANDSCAPE EXPLORER — ÉTAGE 1                       ║
║          Chargement CICY + Filtrage topologique                 ║
╚══════════════════════════════════════════════════════════════════╝

Pipeline :
  1. Charger la base CICY (7890 variétés, ~266 paires de Hodge)
  2. Calculer les invariants dérivés (χ, N_gen, etc.)
  3. Filtrer les candidats Modèle Standard (3 générations)
  4. Produire les visualisations du paysage
  5. Exporter les résultats

Usage :
  python main.py
"""

import os
import sys

# Ajouter le répertoire parent au path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cy_landscape.core.database import CYDatabase
from cy_landscape.core.visualization import (
    plot_hodge_shield,
    plot_generation_distribution,
    plot_sm_candidates_detail,
)


def separator(title: str):
    print(f"\n{'─' * 60}")
    print(f"  {title}")
    print(f"{'─' * 60}")


def main():
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    # ══════════════════════════════════════════════════════════════
    # ÉTAPE 1 : Chargement de la base CICY
    # ══════════════════════════════════════════════════════════════
    separator("ÉTAPE 1 — Chargement de la base CICY")

    db = CYDatabase()
    db.load_cicy(include_mirrors=True)

    summary = db.summary()
    print(f"\n  Base chargée avec succès.")
    print(f"  Sources         : {summary['sources']}")
    print(f"  Topologies      : {summary['topologies_distinctes']}")
    print(f"  Variétés (×mult): {summary['varietes_totales']}")
    print(f"  h¹¹ ∈ {list(summary['h11_range'])}")
    print(f"  h²¹ ∈ {list(summary['h21_range'])}")
    print(f"  χ  ∈ {list(summary['euler_range'])}")
    print(f"  N_gen ∈ {list(summary['n_gen_range'])}")
    print(f"  Auto-miroir     : {summary['auto_miroir']} topologies")

    # ══════════════════════════════════════════════════════════════
    # ÉTAPE 2 : Distribution des générations
    # ══════════════════════════════════════════════════════════════
    separator("ÉTAPE 2 — Distribution du nombre de générations")

    dist = db.generation_distribution()
    print(f"\n  {'N_gen':>5}  {'Topologies':>10}  {'Variétés':>10}  {'% topo':>8}")
    print(f"  {'─'*5}  {'─'*10}  {'─'*10}  {'─'*8}")

    total_topo = db.total_topologies
    for ngen in sorted(dist.keys()):
        d = dist[ngen]
        pct = 100.0 * d["topologies"] / total_topo
        marker = " ◀ SM" if ngen == 3 else ""
        print(f"  {ngen:>5}  {d['topologies']:>10}  {d['varietes']:>10}  {pct:>7.1f}%{marker}")

    # ══════════════════════════════════════════════════════════════
    # ÉTAPE 3 : Filtrage des candidats Modèle Standard
    # ══════════════════════════════════════════════════════════════
    separator("ÉTAPE 3 — Candidats Modèle Standard (3 générations)")

    sm_candidates = db.filter_standard_model_candidates()
    sm_total_var = sum(m.count for m in sm_candidates)

    print(f"\n  Filtre appliqué : |h¹¹ − h²¹| = 3")
    print(f"  Résultat : {len(sm_candidates)} topologies distinctes, "
          f"{sm_total_var} variétés")
    print(f"\n  Taux de sélection : "
          f"{100.0 * len(sm_candidates) / total_topo:.1f}% des topologies")

    # Lister les candidats (groupés par h11)
    print(f"\n  {'h¹¹':>4}  {'h²¹':>4}  {'χ':>5}  {'×count':>7}  {'Source':>12}")
    print(f"  {'─'*4}  {'─'*4}  {'─'*5}  {'─'*7}  {'─'*12}")

    for m in sorted(sm_candidates, key=lambda x: (x.h11, x.h21)):
        print(f"  {m.h11:>4}  {m.h21:>4}  {m.euler:>5}  {m.count:>7}  {m.source:>12}")

    # ══════════════════════════════════════════════════════════════
    # ÉTAPE 4 : Visualisations
    # ══════════════════════════════════════════════════════════════
    separator("ÉTAPE 4 — Génération des visualisations")

    # 4a. Hodge shield complet
    path1 = plot_hodge_shield(
        db, highlight_ngen=3,
        output_path=os.path.join(output_dir, "hodge_shield.png"),
    )
    print(f"  ✓ Hodge shield         → {path1}")

    # 4b. Distribution des générations
    path2 = plot_generation_distribution(
        db,
        output_path=os.path.join(output_dir, "generation_distribution.png"),
    )
    print(f"  ✓ Distribution N_gen   → {path2}")

    # 4c. Détail des candidats SM
    path3 = plot_sm_candidates_detail(
        sm_candidates,
        output_path=os.path.join(output_dir, "sm_candidates.png"),
    )
    print(f"  ✓ Candidats SM détail  → {path3}")

    # ══════════════════════════════════════════════════════════════
    # ÉTAPE 5 : Export des données
    # ══════════════════════════════════════════════════════════════
    separator("ÉTAPE 5 — Export des données")

    csv_path = os.path.join(output_dir, "cicy_database.csv")
    json_path = os.path.join(output_dir, "cicy_database.json")
    db.to_csv(csv_path)
    db.to_json(json_path)
    print(f"  ✓ CSV  → {csv_path}")
    print(f"  ✓ JSON → {json_path}")

    # Export des candidats SM séparément
    sm_csv = os.path.join(output_dir, "sm_candidates.csv")
    import csv
    with open(sm_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["h11", "h21", "euler", "n_gen", "count", "source"])
        for m in sorted(sm_candidates, key=lambda x: (x.h11, x.h21)):
            writer.writerow([m.h11, m.h21, m.euler, m.n_gen, m.count, m.source])
    print(f"  ✓ Candidats SM CSV → {sm_csv}")

    # ══════════════════════════════════════════════════════════════
    # RÉSUMÉ
    # ══════════════════════════════════════════════════════════════
    separator("RÉSUMÉ")
    print(f"""
  ┌─────────────────────────────────────────────────┐
  │  Base CICY chargée : {db.total_topologies:>4} topologies              │
  │  Variétés totales  : {db.total_manifolds:>5}                       │
  │  Candidats SM      : {len(sm_candidates):>4} topologies ({sm_total_var} var.) │
  │                                                 │
  │  Prochaine étape : Étage 2                      │
  │  → Construction de fibrés vectoriels sur        │
  │    les {len(sm_candidates)} topologies candidates               │
  │  → Détermination du groupe de jauge             │
  └─────────────────────────────────────────────────┘
""")



# ======================================================================
# POINT D'ENTREE OBSOLETE  --  refuse de tourner
# ======================================================================
#
# Le seul scan maintenu est `cy_landscape.main_optimized` (§3 du document
# de reference). Ce fichier-ci est anterieur et n'a recu AUCUNE des
# corrections du pipeline :
#
#   - annulation d'anomalie, c2(TX) - c2(V) effective (§5.21) : absente.
#     C'est une condition PHYSIQUE. 70 entrees sur 115 du catalogue
#     `scan_wilson2`, produit sans elle, ne sont pas des modeles ;
#   - phase des twists du critere de Hoppe (§5.15) : absente. Elle a
#     demontre un faux positif du catalogue ;
#   - non-degenerescence des monades (§4.6), c2 croise (§4.1), d_1 de la
#     suite de Koszul (§4.2), bornes rigoureuses sur H^i (§4.3, §4.4),
#     phase 0 de Hoppe (§4.5) : selon les fichiers, absentes ou partielles.
#
# Le laisser executable sans le dire reviendrait a offrir un chemin qui
# produit des resultats d'apparence normale et faux -- exactement le motif
# que le §8 proscrit. Il est conserve pour l'historique, pas pour l'usage.

def _refuser_point_entree_obsolete():
    import sys
    print("Ce point d'entree est OBSOLETE et ne doit plus etre utilise.\n"
          "Il n'a ni l'annulation d'anomalie (condition physique, §5.21),\n"
          "ni la phase des twists de Hoppe (§5.15), ni plusieurs correctifs\n"
          "du §4. Utiliser :\n\n"
          "    python -m cy_landscape.main_optimized cicylist.txt ...\n",
          file=sys.stderr)
    sys.exit(2)


if __name__ == "__main__":
    _refuser_point_entree_obsolete()
    main()
