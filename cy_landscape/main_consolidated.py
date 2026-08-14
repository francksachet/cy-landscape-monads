#!/usr/bin/env python3
"""
PIPELINE CONSOLIDE -- Fondations mathematiques exactes.

Corrections appliquees :
  P1. Nombres d'intersection exacts depuis les matrices de configuration
  P2. Cohomologie par Bott-Borel-Weil + suite de Koszul
  P3. Test de semi-stabilite de Mumford-Takemoto
"""

import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cy_landscape.data.cicy_matrices import (
    get_all_configs, validate_config,
)
from cy_landscape.core.intersection import (
    compute_intersection_numbers,
    compute_euler_from_intersection,
    compute_c2_tangent,
    compute_hodge_favorable,
)
from cy_landscape.core.exact_cohomology import (
    koszul_cohomology,
    bundle_cohomology_exact,
)
from cy_landscape.core.stability import (
    check_semi_stability,
    find_stability_region,
)
from cy_landscape.core.cohomology import (
    extract_spectrum_su5, extract_spectrum_so10, extract_spectrum_e6,
    ParticleSpectrum, FullAnalysis,
)
from cy_landscape.core.bundles import (
    CICYGeometry, LineBundleSum, BundleAnalysis, GAUGE_GROUP_TABLE,
)


def sep(title):
    w = 66
    print(f"\n{'=' * w}")
    print(f"  {title}")
    print(f"{'=' * w}")


def main():
    out = "output_consolidated"
    os.makedirs(out, exist_ok=True)

    # =====================================================
    # PHASE 1 : Validation et calculs exacts de geometrie
    # =====================================================
    sep("PHASE 1 -- Geometrie exacte depuis les matrices de config")

    configs = get_all_configs()
    valid_geoms = []

    print(f"\n  {'Nom':<25} {'Amb':>5} {'K':>3} {'dim':>4} {'CY':>4} "
          f"{'h11':>4} {'h21':>4} {'chi':>6} {'|dh|':>4}")
    print(f"  {'='*25} {'='*5} {'='*3} {'='*4} {'='*4} "
          f"{'='*4} {'='*4} {'='*6} {'='*4}")

    for cfg in configs:
        v = validate_config(cfg)
        if not v["valid"]:
            print(f"  {cfg['name']:<25} INVALIDE: dim={v['dim_CY']}, CY={v['cy_condition']}")
            continue

        ambient = cfg["ambient"]
        config = cfg["config"]
        m = len(ambient)

        # Calcul EXACT des nombres d'intersection
        try:
            d_ijk = compute_intersection_numbers(ambient, config)
        except Exception as e:
            print(f"  {cfg['name']:<25} ERREUR intersection: {e}")
            continue

        # Verification: tous les d_ijk doivent etre >= 0 (effectivite)
        # (pas toujours vrai, certains peuvent etre negatifs pour les non-favorables)

        # Calcul de chi et Hodge
        try:
            chi = compute_euler_from_intersection(ambient, config, d_ijk)
            h11 = m  # Favorable
            h21 = h11 - chi // 2
        except Exception as e:
            print(f"  {cfg['name']:<25} ERREUR Euler: {e}")
            continue

        # c2(TY)
        try:
            c2 = compute_c2_tangent(ambient, config, d_ijk)
        except Exception as e:
            c2 = np.ones(m) * 24  # Fallback
            print(f"  {cfg['name']:<25} WARN c2: fallback")

        # Verification des Hodge connus
        status = ""
        if cfg.get("h_known"):
            h_exp = cfg["h_known"]
            if h11 == h_exp[0] and h21 == h_exp[1]:
                status = " OK"
            else:
                status = f" !!(attendu {h_exp})"

        n_gen = abs(h11 - h21)
        amb_str = "x".join(f"P{n}" for n in ambient)
        print(f"  {cfg['name']:<25} {amb_str:>5} {config.shape[0]:>3} "
              f"{3:>4} {'oui':>4} {h11:>4} {h21:>4} {chi:>6} {n_gen:>4}{status}")

        # Construire l'objet geometrie
        geom = CICYGeometry(
            ambient_dims=ambient,
            config_matrix=config,
            h11=h11, h21=h21,
            intersection_numbers=d_ijk,
            c2_tangent=c2,
        )
        valid_geoms.append((cfg, geom))

    print(f"\n  Total : {len(valid_geoms)} geometries valides sur {len(configs)}")

    # Filtrer les candidats SM (|h11 - h21| <= 5 pour inclure les proches)
    sm_geoms = [(c, g) for c, g in valid_geoms if abs(g.h11 - g.h21) <= 5]
    exact3 = [(c, g) for c, g in valid_geoms if abs(g.h11 - g.h21) == 3]
    print(f"  Candidats |dh| <= 5 : {len(sm_geoms)}")
    print(f"  Candidats |dh| = 3  : {len(exact3)}")

    if not sm_geoms:
        sm_geoms = valid_geoms[:8]
        print(f"  (aucun a 3 gen -> analyse des {len(sm_geoms)} premiers)")

    # =====================================================
    # PHASE 2 : Scan de fibres avec cohomologie exacte
    # =====================================================
    sep("PHASE 2 -- Scan de fibres + cohomologie Koszul exacte")

    all_analyses = []

    for cfg, geom in sm_geoms:
        m = geom.h11
        print(f"\n  --- {cfg['name']} : CY({geom.h11},{geom.h21}), "
              f"chi={geom.euler}, |dh|={geom.n_gen} ---")

        # Generer des fibres de rang 3, 4, 5
        n_tested = 0
        n_valid = 0

        for rank in [3, 4, 5]:
            gauge_info = GAUGE_GROUP_TABLE.get(rank, {})
            gauge_group = gauge_info.get("group", f"?")

            bundles = _generate_bundles(m, rank, charge_max=2)

            for charges in bundles:
                n_tested += 1
                bundle = LineBundleSum(charges=charges)

                if not bundle.c1_vanishes:
                    continue

                # Cohomologie EXACTE (Koszul)
                try:
                    cohom = bundle_cohomology_exact(
                        geom.ambient_dims,
                        geom.config_matrix,
                        charges,
                    )
                except Exception:
                    continue

                # Extraire le spectre
                if gauge_group == "SU(5)":
                    spectrum = extract_spectrum_su5(cohom)
                elif gauge_group == "SO(10)":
                    spectrum = extract_spectrum_so10(cohom)
                elif gauge_group in ("E₆", "E6"):
                    spectrum = extract_spectrum_e6(cohom)
                else:
                    continue

                # Stabilite (P3)
                stab = check_semi_stability(charges, geom.intersection_numbers)

                # Construire l'analyse du fibre
                ba = BundleAnalysis(geometry=geom, bundle=bundle).compute()

                # Analyse complete
                fa = FullAnalysis(bundle_analysis=ba, spectrum=spectrum)
                fa.compute_final_score()

                # Bonus/malus stabilite
                if stab["semi_stable"]:
                    fa.final_score = min(100, fa.final_score + 5)
                    stab_flag = "S"
                else:
                    fa.final_score = max(0, fa.final_score - 10)
                    stab_flag = "U"

                fa._stability = stab
                fa._config_name = cfg["name"]

                all_analyses.append(fa)
                n_valid += 1

        print(f"    Testes: {n_tested}, retenus: {n_valid}")

    # =====================================================
    # PHASE 3 : Classement et resultats
    # =====================================================
    sep("PHASE 3 -- Classement final avec fondations exactes")

    all_analyses.sort(key=lambda a: a.final_score, reverse=True)

    if not all_analyses:
        print("\n  Aucun fibre viable trouve.")
        print("  Les geometries incluses n'ont peut-etre pas |h11-h21|=3.")
        print("  Etendre la base de config matrices pour plus de resultats.")
        _save_geometry_report(valid_geoms, out)
        return

    # Top resultats
    top = all_analyses[:min(20, len(all_analyses))]

    print(f"\n  {'#':>3}  {'Config':<18} {'CY':>10} {'Rk':>3} {'Jauge':>7} "
          f"{'Ngen':>4} {'Hig':>3} {'Exo':>3} {'Stab':>4} {'Score':>5}")
    print(f"  {'='*3}  {'='*18} {'='*10} {'='*3} {'='*7} "
          f"{'='*4} {'='*3} {'='*3} {'='*4} {'='*5}")

    for i, fa in enumerate(top):
        ba = fa.bundle_analysis
        sp = fa.spectrum
        stab = "S" if getattr(fa, '_stability', {}).get('semi_stable', False) else "U"
        name = getattr(fa, '_config_name', '?')
        print(f"  {i+1:>3}  {name:<18} ({ba.geometry.h11:>3},{ba.geometry.h21:>3}) "
              f"SU{ba.structure_group_rank} {sp.gauge_group:>7} "
              f"{sp.n_generations:>4} {sp.n_higgs_candidates:>3} "
              f"{sp.n_exotics:>3} {stab:>4} {fa.final_score:>5.0f}")

    # Detail du meilleur
    if top:
        sep("MEILLEUR CANDIDAT -- DETAIL COMPLET")
        best = top[0]
        ba = best.bundle_analysis
        sp = best.spectrum
        stab = getattr(best, '_stability', {})
        name = getattr(best, '_config_name', '?')

        print(f"\n  Configuration : {name}")
        print(f"  Espace ambiant: {'x'.join(f'P^{n}' for n in ba.geometry.ambient_dims)}")
        print(f"  Matrice config:")
        for row in ba.geometry.config_matrix:
            print(f"    {row.tolist()}")

        print(f"\n  Nombres d'intersection d_ijk (non nuls):")
        d = ba.geometry.intersection_numbers
        for i in range(d.shape[0]):
            for j in range(i, d.shape[1]):
                for k in range(j, d.shape[2]):
                    if d[i,j,k] != 0:
                        print(f"    d_{{{i}{j}{k}}} = {int(d[i,j,k])}")

        print(f"\n  c2(TY) . J_i = {ba.geometry.c2_tangent.tolist()}")
        print(f"  h11 = {ba.geometry.h11}, h21 = {ba.geometry.h21}, "
              f"chi = {ba.geometry.euler}")

        print(f"\n  Fibre : rang {ba.structure_group_rank}")
        print(f"  Charges :")
        for q in ba.bundle.charges:
            print(f"    {q}")
        print(f"  c1(V) = {ba.bundle.c1.tolist()} "
              f"({'= 0 ok' if ba.c1_vanishes else '!= 0 ERREUR'})")

        print(f"\n  Groupe de jauge : {sp.gauge_group}")

        print(f"\n  Spectre (cohomologie Koszul exacte) :")
        for rep, count in sp.representations.items():
            if count > 0:
                print(f"    {count:>3} x {rep}")

        gen_ok = "oui" if sp.generation_match else "NON"
        hig_ok = "oui" if sp.higgs_present else "NON"
        exo_ok = "oui" if sp.exotic_free else "NON"
        stab_ok = "oui" if stab.get('semi_stable', False) else "NON"

        print(f"\n  Validation :")
        print(f"    3 generations  : {gen_ok} (N_gen = {sp.n_generations})")
        print(f"    Higgs present  : {hig_ok} ({sp.n_higgs_candidates} candidat(s))")
        print(f"    Sans exotiques : {exo_ok} ({sp.n_exotics} exotiques)")
        print(f"    Semi-stable    : {stab_ok}")
        if not stab.get('semi_stable', False):
            print(f"      Violation  : {stab.get('worst_violation', '?')}")
            print(f"      Sous-ens.  : {stab.get('destabilizing_subset', '?')}")

        print(f"\n  Score final : {best.final_score:.0f}/100")

    # Export
    _save_results(all_analyses, valid_geoms, out)

    # Resume
    sep("RESUME -- CONSOLIDATION")
    n_stable = sum(1 for a in all_analyses
                   if getattr(a, '_stability', {}).get('semi_stable', False))
    n_3gen = sum(1 for a in all_analyses if a.spectrum.generation_match)

    print(f"""
  +------------------------------------------------------------+
  |  PIPELINE CONSOLIDE -- FONDATIONS EXACTES                  |
  +------------------------------------------------------------+
  |                                                            |
  |  P1. Nombres d'intersection : EXACTS (depuis config)       |
  |  P2. Cohomologie : EXACTE (Bott-Borel-Weil + Koszul)      |
  |  P3. Stabilite : VERIFIEE (Mumford-Takemoto)              |
  |                                                            |
  |  Geometries analysees  : {len(valid_geoms):>4}                            |
  |  Fibres testes         : {len(all_analyses):>5}                           |
  |  Semi-stables          : {n_stable:>5}                           |
  |  3 generations exactes : {n_3gen:>5}                           |""")
    if top:
        print(f"  |  Meilleur score        : {top[0].final_score:>5.0f}/100"
              f"                       |")
    print(f"  +------------------------------------------------------------+\n")


def _generate_bundles(m, rank, charge_max=2):
    """Genere des sommes de line bundles avec c1=0 automatique."""
    from itertools import product as iprod
    bundles = []
    crange = list(range(-charge_max, charge_max + 1))
    m_eff = min(m, 3)  # Limiter pour tractabilite

    # Fibres structures
    for shift in range(min(m, 3)):
        charges = []
        for r in range(rank - 1):
            q = [0] * m
            idx = (r + shift) % m
            q[idx] = 1
            charges.append(q)
        last = [-sum(q[i] for q in charges) for i in range(m)]
        if all(-charge_max-1 <= l <= charge_max+1 for l in last):
            charges.append(last)
            bundles.append(charges)

    # Fibres (1,-1) dans differentes directions
    for i in range(min(m, 3)):
        for j in range(i+1, min(m, 4)):
            charges = []
            for r in range(rank - 1):
                q = [0] * m
                q[i] = 1 if r % 2 == 0 else -1
                q[j] = -1 if r % 2 == 0 else 1
                charges.append(q)
            last = [-sum(q[k] for q in charges) for k in range(m)]
            charges.append(last)
            bundles.append(charges)

    # Echantillonnage aleatoire
    rng = np.random.RandomState(42)
    for _ in range(300):
        charges = []
        for r in range(rank - 1):
            q = [0] * m
            for d in range(m_eff):
                q[d] = int(rng.randint(-charge_max, charge_max + 1))
            charges.append(q)
        last = [-sum(q[k] for q in charges) for k in range(m)]
        if all(-charge_max-1 <= l <= charge_max+1 for l in last):
            charges.append(last)
            bundles.append(charges)

    return bundles


def _save_results(analyses, geoms, out_dir):
    """Sauvegarde JSON."""
    data = {
        "pipeline": "CY Landscape - Pipeline consolide",
        "corrections": [
            "P1: Nombres d'intersection exacts",
            "P2: Cohomologie Bott-Borel-Weil + Koszul",
            "P3: Semi-stabilite Mumford-Takemoto",
        ],
        "geometries": [
            {"name": c["name"], "ambient": c["ambient"],
             "h11": g.h11, "h21": g.h21, "euler": g.euler}
            for c, g in geoms
        ],
        "top_results": [a.to_dict() for a in analyses[:20]],
    }
    path = os.path.join(out_dir, "consolidated_results.json")
    with open(path, "w") as f:
        json.dump(data, f, indent=2, default=str)
    print(f"\n  Export: {path}")


def _save_geometry_report(geoms, out_dir):
    """Sauvegarde un rapport des geometries meme sans fibres."""
    data = {
        "geometries": [
            {"name": c["name"], "ambient": c["ambient"],
             "h11": g.h11, "h21": g.h21, "euler": g.euler,
             "n_gen": g.n_gen,
             "d_ijk_nonzero": int(np.count_nonzero(g.intersection_numbers)),
             "c2": g.c2_tangent.tolist()}
            for c, g in geoms
        ]
    }
    path = os.path.join(out_dir, "geometry_report.json")
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"\n  Rapport geometrique: {path}")



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
