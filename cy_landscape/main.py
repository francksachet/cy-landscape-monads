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


if __name__ == "__main__":
    main()
