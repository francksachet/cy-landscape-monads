#!/usr/bin/env python3
"""
ETAGE 3 — Cohomologie des fibres et spectre de particules.

Pipeline :
  1. Recharger les meilleurs fibres de l'etage 2
  2. Calculer la cohomologie de chaque fibre
  3. Extraire le spectre de matiere (representations du groupe de jauge)
  4. Comparer au contenu du Modele Standard
  5. Score final combine et classement
"""

import os, sys, json, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cy_landscape.core.database import CYDatabase
from cy_landscape.core.geometry import generate_model_geometry
from cy_landscape.core.bundles import BundleAnalysis
from cy_landscape.core.scanner import BundleScanner, ScanConfig
from cy_landscape.core.cohomology import (
    extract_spectrum, FullAnalysis, bundle_cohomology,
)
from cy_landscape.core.viz_stage3 import (
    plot_particle_spectrum,
    plot_sm_comparison,
    plot_final_ranking,
    plot_full_pipeline,
)


def sep(title):
    print(f"\n{'=' * 64}")
    print(f"  {title}")
    print(f"{'=' * 64}")


def main():
    out = "output_stage3"
    os.makedirs(out, exist_ok=True)

    # =====================================================
    # ETAPE 1 : Recharger et re-scanner (etages 1+2)
    # =====================================================
    sep("ETAPE 1 - Reconstruction du pipeline (etages 1+2)")

    db = CYDatabase()
    db.load_cicy(include_mirrors=True)
    sm_candidates = db.filter_standard_model_candidates()

    # Uniques modulo miroir
    unique = []
    seen = set()
    for m in sm_candidates:
        pair = tuple(sorted([m.h11, m.h21]))
        if pair not in seen:
            unique.append(m)
            seen.add(pair)

    print(f"  {len(unique)} geometries CY candidates (modulo miroir)")

    # Generer geometries
    geometries = []
    for i, cy in enumerate(unique):
        geom = generate_model_geometry(cy.h11, cy.h21, seed=42 + i)
        geometries.append(geom)

    # Scanner les fibres
    config = ScanConfig(
        bundle_ranks=[3, 4, 5],
        charge_max=3,
        max_bundles_per_rank=500,
        min_score=30.0,
    )
    scanner = BundleScanner(config)
    for geom in geometries:
        scanner.scan_geometry(geom)

    # Garder les meilleurs de l'etage 2
    top_bundles = scanner.top_results(50)
    print(f"  {len(top_bundles)} meilleurs fibres retenus pour l'analyse spectrale")

    # =====================================================
    # ETAPE 2 : Cohomologie et spectre de matiere
    # =====================================================
    sep("ETAPE 2 - Calcul de la cohomologie et extraction du spectre")

    full_analyses = []
    for i, ba in enumerate(top_bundles):
        # Calculer le spectre
        spectrum = extract_spectrum(ba)

        # Creer l'analyse complete
        fa = FullAnalysis(
            bundle_analysis=ba,
            spectrum=spectrum,
        )
        fa.compute_final_score()
        full_analyses.append(fa)

    # Trier par score final
    full_analyses.sort(key=lambda a: a.final_score, reverse=True)

    print(f"\n  Analyses completes : {len(full_analyses)}")
    print(f"  Score final max    : {full_analyses[0].final_score:.1f}")
    print(f"  Score final min    : {full_analyses[-1].final_score:.1f}")

    # =====================================================
    # ETAPE 3 : Resultats detailles
    # =====================================================
    sep("ETAPE 3 - Resultats detailles")

    # Top 15
    print(f"\n  {'#':>3}  {'CY':>10}  {'Jauge':>8}  {'N_gen':>5}  "
          f"{'Higgs':>5}  {'Exot':>4}  {'S.fibre':>7}  {'S.spec':>6}  {'FINAL':>6}")
    print(f"  {'='*3}  {'='*10}  {'='*8}  {'='*5}  {'='*5}  {'='*4}  {'='*7}  {'='*6}  {'='*6}")

    for i, fa in enumerate(full_analyses[:15]):
        ba = fa.bundle_analysis
        sp = fa.spectrum
        gen_ok = "ok" if sp.generation_match else "  "
        print(f"  {i+1:>3}  ({ba.geometry.h11:>3},{ba.geometry.h21:>3})  "
              f"{sp.gauge_group:>8}  {sp.n_generations:>3}{gen_ok}  "
              f"{sp.n_higgs_candidates:>5}  {sp.n_exotics:>4}  "
              f"{ba.sm_score:>6.0f}  {sp.sm_compatibility:>5.0f}  "
              f"{fa.final_score:>5.0f}")

    # Detail du meilleur
    sep("DETAIL DU MEILLEUR CANDIDAT")
    best = full_analyses[0]
    print(f"\n{best.summary()}")

    # Cohomologie detaillee du meilleur
    cohom = bundle_cohomology(best.bundle_analysis.bundle,
                               best.bundle_analysis.geometry)
    print(f"\n  Cohomologie detaillee :")
    for key in ["V", "V_dual", "wedge2V", "end_V"]:
        c = cohom[key]
        print(f"    H^i(Y, {key:>7}) = "
              f"[h0={c[0]:>3}, h1={c[1]:>3}, h2={c[2]:>3}, h3={c[3]:>3}]")

    # Statistiques par groupe de jauge
    sep("STATISTIQUES PAR GROUPE DE JAUGE")
    groups = {}
    for fa in full_analyses:
        g = fa.spectrum.gauge_group
        if g not in groups:
            groups[g] = {"count": 0, "scores": [], "gen3": 0, "higgs": 0}
        groups[g]["count"] += 1
        groups[g]["scores"].append(fa.final_score)
        if fa.spectrum.generation_match:
            groups[g]["gen3"] += 1
        if fa.spectrum.higgs_present:
            groups[g]["higgs"] += 1

    print(f"\n  {'Groupe':>8}  {'Total':>5}  {'3-gen':>5}  {'Higgs':>5}  "
          f"{'Score moy':>9}  {'Score max':>9}")
    print(f"  {'='*8}  {'='*5}  {'='*5}  {'='*5}  {'='*9}  {'='*9}")
    for g in sorted(groups.keys()):
        d = groups[g]
        print(f"  {g:>8}  {d['count']:>5}  {d['gen3']:>5}  {d['higgs']:>5}  "
              f"{np.mean(d['scores']):>8.1f}  {max(d['scores']):>8.1f}")

    # =====================================================
    # ETAPE 4 : Visualisations
    # =====================================================
    sep("ETAPE 4 - Visualisations")

    p1 = plot_particle_spectrum(full_analyses,
                                 output_path=os.path.join(out, "particle_spectra.png"))
    if p1: print(f"  ok Spectres de particules  -> {p1}")

    p2 = plot_sm_comparison(full_analyses[0],
                             output_path=os.path.join(out, "sm_comparison.png"))
    if p2: print(f"  ok Comparaison SM          -> {p2}")

    p3 = plot_final_ranking(full_analyses,
                             output_path=os.path.join(out, "final_ranking.png"))
    if p3: print(f"  ok Classement final        -> {p3}")

    p4 = plot_full_pipeline(full_analyses[0],
                             output_path=os.path.join(out, "full_pipeline.png"))
    if p4: print(f"  ok Pipeline complet        -> {p4}")

    # =====================================================
    # ETAPE 5 : Export
    # =====================================================
    sep("ETAPE 5 - Export")

    # JSON
    export = {
        "pipeline": "CY Landscape Explorer - Etage 3",
        "total_analyses": len(full_analyses),
        "top_results": [fa.to_dict() for fa in full_analyses[:20]],
    }
    jp = os.path.join(out, "stage3_results.json")
    with open(jp, "w") as f:
        json.dump(export, f, indent=2)
    print(f"  ok JSON -> {jp}")

    # CSV
    cp = os.path.join(out, "final_ranking.csv")
    with open(cp, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["rank", "h11", "h21", "euler", "bundle_rank",
                     "gauge_group", "n_gen", "gen_match", "n_higgs",
                     "n_exotics", "score_bundle", "score_spectrum",
                     "score_final"])
        for i, fa in enumerate(full_analyses[:20]):
            ba = fa.bundle_analysis
            sp = fa.spectrum
            w.writerow([i+1, ba.geometry.h11, ba.geometry.h21,
                        ba.geometry.euler, ba.structure_group_rank,
                        sp.gauge_group, sp.n_generations,
                        sp.generation_match, sp.n_higgs_candidates,
                        sp.n_exotics, round(ba.sm_score, 1),
                        round(sp.sm_compatibility, 1),
                        round(fa.final_score, 1)])
    print(f"  ok CSV  -> {cp}")

    # =====================================================
    # RESUME FINAL
    # =====================================================
    sep("RESUME FINAL — PIPELINE COMPLET")

    n_gen3 = sum(1 for fa in full_analyses if fa.spectrum.generation_match)
    n_higgs = sum(1 for fa in full_analyses if fa.spectrum.higgs_present)
    n_clean = sum(1 for fa in full_analyses if fa.spectrum.exotic_free)
    n_excellent = sum(1 for fa in full_analyses if fa.final_score >= 70)

    best = full_analyses[0]
    ba = best.bundle_analysis
    sp = best.spectrum

    print(f"""
  +----------------------------------------------------------+
  |  CY LANDSCAPE EXPLORER — RESULTATS COMPLETS              |
  +----------------------------------------------------------+
  |                                                          |
  |  Etage 1 : {len(unique):>4} topologies CY a 3 generations          |
  |  Etage 2 : {len(top_bundles):>4} fibres vectoriels viables              |
  |  Etage 3 : {len(full_analyses):>4} spectres de matiere calcules          |
  |                                                          |
  |  Avec 3 generations exactes : {n_gen3:>4}                      |
  |  Avec Higgs               : {n_higgs:>4}                      |
  |  Sans exotiques            : {n_clean:>4}                      |
  |  Score final >= 70         : {n_excellent:>4}                      |
  |                                                          |
  |  MEILLEUR CANDIDAT :                                     |
  |    CY({ba.geometry.h11},{ba.geometry.h21}) + SU({ba.structure_group_rank}) -> {sp.gauge_group:>6}               |
  |    {sp.n_generations} generations, {sp.n_higgs_candidates} Higgs, {sp.n_exotics} exotiques            |
  |    Score final : {best.final_score:.0f}/100                              |
  +----------------------------------------------------------+
""")


if __name__ == "__main__":
    import numpy as np
    main()
