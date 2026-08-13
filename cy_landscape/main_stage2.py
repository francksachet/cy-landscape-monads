#!/usr/bin/env python3
"""
╔══════════════════════════════════════════════════════════════════╗
║          CY LANDSCAPE EXPLORER — ÉTAGE 2                       ║
║          Construction de fibrés vectoriels                      ║
║          Détermination du groupe de jauge                       ║
╚══════════════════════════════════════════════════════════════════╝

Pipeline :
  1. Charger les candidats SM de l'étage 1
  2. Générer les géométries détaillées
  3. Scanner l'espace des fibrés (sommes de line bundles)
  4. Analyser : groupe de jauge, indice chiral, anomalies
  5. Classer et visualiser les résultats
"""

import os
import sys
import json
import csv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cy_landscape.core.database import CYDatabase
from cy_landscape.core.geometry import generate_model_geometry
from cy_landscape.core.bundles import BundleAnalysis, GAUGE_GROUP_TABLE
from cy_landscape.core.scanner import BundleScanner, ScanConfig
from cy_landscape.core.viz_stage2 import (
    plot_gauge_group_distribution,
    plot_score_landscape,
    plot_bundle_summary_table,
    plot_physics_chain,
)


def separator(title: str):
    print(f"\n{'═' * 64}")
    print(f"  {title}")
    print(f"{'═' * 64}")


def main():
    output_dir = "output_stage2"
    os.makedirs(output_dir, exist_ok=True)

    # ══════════════════════════════════════════════════════════════
    # ÉTAPE 1 : Recharger les candidats SM de l'étage 1
    # ══════════════════════════════════════════════════════════════
    separator("ÉTAPE 1 — Chargement des candidats SM (étage 1)")

    db = CYDatabase()
    db.load_cicy(include_mirrors=True)
    sm_candidates = db.filter_standard_model_candidates()

    # Ne garder que les variétés originales (pas les miroirs, qui sont redondantes)
    unique_candidates = []
    seen = set()
    for m in sm_candidates:
        pair = tuple(sorted([m.h11, m.h21]))
        if pair not in seen:
            unique_candidates.append(m)
            seen.add(pair)

    print(f"\n  Candidats SM uniques (modulo miroir) : {len(unique_candidates)}")
    for m in unique_candidates:
        print(f"    CY({m.h11}, {m.h21})  χ={m.euler:+d}  ×{m.count}")

    # ══════════════════════════════════════════════════════════════
    # ÉTAPE 2 : Générer les géométries détaillées
    # ══════════════════════════════════════════════════════════════
    separator("ÉTAPE 2 — Génération des géométries CY")

    geometries = []
    for i, cy in enumerate(unique_candidates):
        geom = generate_model_geometry(cy.h11, cy.h21, seed=42 + i)
        geometries.append(geom)
        print(f"  CY({geom.h11},{geom.h21}) : "
              f"espace ambiant P^{' × P^'.join(str(d) for d in geom.ambient_dims)}, "
              f"h11={geom.h11}, d_ijk shape={geom.intersection_numbers.shape}")

    # ══════════════════════════════════════════════════════════════
    # ÉTAPE 3 : Scan des fibrés vectoriels
    # ══════════════════════════════════════════════════════════════
    separator("ÉTAPE 3 — Scan de l'espace des fibrés")

    print(f"\n  Stratégie :")
    print(f"    • Rangs explorés : SU(3) → E₆, SU(4) → SO(10), SU(5) → SU(5)")
    print(f"    • Type : sommes de fibrés en droites (line bundle sums)")
    print(f"    • Charges : entiers dans [-3, +3]")
    print(f"    • Contrainte automatique : c₁(V) = 0")
    print()

    config = ScanConfig(
        bundle_ranks=[3, 4, 5],
        charge_max=3,
        max_bundles_per_rank=500,
        min_score=30.0,
        target_generations=3,
    )

    scanner = BundleScanner(config)

    for i, geom in enumerate(geometries):
        results = scanner.scan_geometry(geom)
        n_good = sum(1 for r in results if r.sm_score >= 50)
        print(f"  CY({geom.h11},{geom.h21}) : "
              f"{len(results)} fibrés retenus, {n_good} avec score ≥ 50")

    stats = scanner.stats
    print(f"\n  ── Statistiques du scan ──")
    print(f"  Géométries scannées  : {stats['geometries_scanned']}")
    print(f"  Fibrés testés        : {stats['bundles_tested']}")
    print(f"  c₁ = 0 satisfait     : {stats['bundles_valid_c1']}")
    print(f"  3 générations        : {stats['bundles_3gen']}")
    print(f"  Anomalie OK          : {stats['bundles_anomaly_ok']}")
    print(f"  Retenus (score≥30)   : {stats['bundles_retained']}")

    # ══════════════════════════════════════════════════════════════
    # ÉTAPE 4 : Analyse des résultats
    # ══════════════════════════════════════════════════════════════
    separator("ÉTAPE 4 — Analyse des résultats")

    all_results = scanner.results
    if not all_results:
        print("\n  ⚠ Aucun fibré viable trouvé. Élargir les paramètres de scan.")
        return

    # Par groupe de jauge
    by_group = scanner.results_by_gauge_group()
    print(f"\n  Distribution par groupe de jauge :")
    print(f"  {'Groupe':>8}  {'Fibrés':>8}  {'Score moyen':>12}  {'Score max':>10}")
    print(f"  {'─'*8}  {'─'*8}  {'─'*12}  {'─'*10}")
    for group in sorted(by_group.keys()):
        rs = by_group[group]
        scores = [r.sm_score for r in rs]
        print(f"  {group:>8}  {len(rs):>8}  {sum(scores)/len(scores):>11.1f}  "
              f"{max(scores):>9.1f}")

    # Top 15
    top = scanner.top_results(15)
    print(f"\n  ── Top 15 fibrés ──")
    print(f"  {'#':>3}  {'CY':>12}  {'Rang':>5}  {'Groupe':>8}  "
          f"{'c₁=0':>5}  {'N_gen':>6}  {'Anom':>5}  {'Score':>6}")
    print(f"  {'─'*3}  {'─'*12}  {'─'*5}  {'─'*8}  {'─'*5}  {'─'*6}  {'─'*5}  {'─'*6}")

    for i, r in enumerate(top):
        ck = lambda b: "  ✓" if b else "  ✗"
        print(f"  {i+1:>3}  ({r.geometry.h11:>3},{r.geometry.h21:>3})  "
              f"SU({r.structure_group_rank})  {r.gauge_group:>8}  "
              f"{ck(r.c1_vanishes)}  {r.chiral_index_value:>5.1f}  "
              f"{ck(r.anomaly_cancelled)}  {r.sm_score:>5.1f}")

    # ══════════════════════════════════════════════════════════════
    # ÉTAPE 5 : Visualisations
    # ══════════════════════════════════════════════════════════════
    separator("ÉTAPE 5 — Visualisations")

    # 5a. Distribution des groupes de jauge
    p1 = plot_gauge_group_distribution(
        all_results,
        output_path=os.path.join(output_dir, "gauge_groups.png"),
    )
    if p1:
        print(f"  ✓ Groupes de jauge    → {p1}")

    # 5b. Paysage des scores
    p2 = plot_score_landscape(
        all_results,
        output_path=os.path.join(output_dir, "score_landscape.png"),
    )
    if p2:
        print(f"  ✓ Paysage des scores  → {p2}")

    # 5c. Tableau des meilleurs
    p3 = plot_bundle_summary_table(
        top,
        output_path=os.path.join(output_dir, "top_bundles.png"),
    )
    if p3:
        print(f"  ✓ Top fibrés (table)  → {p3}")

    # 5d. Chaîne physique du meilleur candidat
    if top:
        p4 = plot_physics_chain(
            top[0],
            output_path=os.path.join(output_dir, "physics_chain.png"),
        )
        if p4:
            print(f"  ✓ Chaîne physique     → {p4}")

    # ══════════════════════════════════════════════════════════════
    # ÉTAPE 6 : Export
    # ══════════════════════════════════════════════════════════════
    separator("ÉTAPE 6 — Export des données")

    # JSON complet
    export_data = {
        "scan_config": {
            "bundle_ranks": config.bundle_ranks,
            "charge_max": config.charge_max,
            "min_score": config.min_score,
        },
        "statistics": stats,
        "results": [r.to_dict() for r in all_results],
    }
    json_path = os.path.join(output_dir, "stage2_results.json")
    with open(json_path, "w") as f:
        json.dump(export_data, f, indent=2)
    print(f"  ✓ JSON complet → {json_path}")

    # CSV des top résultats
    csv_path = os.path.join(output_dir, "top_bundles.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["rank", "h11", "h21", "euler", "bundle_rank",
                         "gauge_group", "c1_zero", "chiral_index",
                         "anomaly_ok", "sm_score"])
        for i, r in enumerate(top):
            writer.writerow([
                i + 1, r.geometry.h11, r.geometry.h21, r.geometry.euler,
                r.structure_group_rank, r.gauge_group,
                r.c1_vanishes, round(r.chiral_index_value, 2),
                r.anomaly_cancelled, round(r.sm_score, 1),
            ])
    print(f"  ✓ Top CSV     → {csv_path}")

    # ══════════════════════════════════════════════════════════════
    # RÉSUMÉ
    # ══════════════════════════════════════════════════════════════
    separator("RÉSUMÉ ÉTAGE 2")

    best = top[0] if top else None
    n_viable = sum(1 for r in all_results if r.sm_score >= 70)

    print(f"""
  ┌──────────────────────────────────────────────────────────┐
  │  Géométries analysées : {stats['geometries_scanned']:>4}                              │
  │  Fibrés testés        : {stats['bundles_tested']:>6}                            │
  │  Fibrés retenus       : {stats['bundles_retained']:>6} (score ≥ 30)               │
  │  Hautement viables    : {n_viable:>6} (score ≥ 70)               │""")

    if best:
        print(f"  │                                                          │")
        print(f"  │  ★ Meilleur candidat :                                   │")
        print(f"  │    CY({best.geometry.h11},{best.geometry.h21}) + "
              f"SU({best.structure_group_rank}) → {best.gauge_group:<6}                     │")
        print(f"  │    Score SM = {best.sm_score:.0f}/100                "
              f"                       │")

    print(f"""  │                                                          │
  │  Prochaine étape : Étage 3                               │
  │  → Calcul de la cohomologie des fibrés                   │
  │  → Extraction du spectre de particules                   │
  │  → Comparaison au contenu du Modèle Standard             │
  └──────────────────────────────────────────────────────────┘
""")


if __name__ == "__main__":
    main()
