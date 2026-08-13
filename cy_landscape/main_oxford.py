#!/usr/bin/env python3
"""
Pipeline CY Landscape avec donnees Oxford.

Usage:
  python -m cy_landscape.main_oxford                           # donnees embarquees (14 CICYs)
  python -m cy_landscape.main_oxford cicylist.txt              # fichier Oxford complet (7890)
  python -m cy_landscape.main_oxford cicylist.txt --max-ps 6   # limiter la taille des espaces ambiants
  python -m cy_landscape.main_oxford cicylist.txt --chi -4     # filtrer par chi
"""
import os, sys, json, argparse, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cy_landscape.core.intersection import (
    compute_intersection_numbers, compute_euler_from_intersection, compute_c2_tangent)
from cy_landscape.core.exact_cohomology import bundle_cohomology_exact
from cy_landscape.core.stability import check_semi_stability
from cy_landscape.core.cohomology import extract_spectrum_su5, extract_spectrum_so10, extract_spectrum_e6
from cy_landscape.core.bundles import CICYGeometry, LineBundleSum, BundleAnalysis, GAUGE_GROUP_TABLE


def load_cicys(args):
    """Charge les CICYs selon la source choisie."""
    if args.file:
        from cy_landscape.data.parse_oxford import load_oxford_file
        print(f"  Chargement de {args.file}...")
        entries = load_oxford_file(args.file)
        print(f"  {len(entries)} CICYs parsees depuis Oxford")
    else:
        from cy_landscape.data.oxford_cicys import get_all_oxford
        entries = get_all_oxford()
        print(f"  {len(entries)} CICYs embarquees (echantillon)")

    # Filtres
    if args.max_ps:
        entries = [e for e in entries if len(e['ambient']) <= args.max_ps]
        print(f"  Filtre max_ps<={args.max_ps}: {len(entries)} restantes")
    if args.chi is not None:
        entries = [e for e in entries if e['chi'] == args.chi]
        print(f"  Filtre chi={args.chi}: {len(entries)} restantes")

    return entries


def validate_cicys(entries):
    """Valide les CICYs: intersection numbers, chi, c2."""
    valid = []
    n_fail = 0
    for c in entries:
        try:
            d = compute_intersection_numbers(c['ambient'], c['config'])
            chi = compute_euler_from_intersection(c['ambient'], c['config'], d)
            c2 = compute_c2_tangent(c['ambient'], c['config'], d)
            if chi != c['chi']:
                n_fail += 1
                continue
            geom = CICYGeometry(
                ambient_dims=c['ambient'], config_matrix=c['config'],
                h11=c['h11'], h21=c['h21'],
                intersection_numbers=d, c2_tangent=c2)
            valid.append((c, geom))
        except Exception:
            n_fail += 1
    return valid, n_fail


def scan_bundles(geom, c, rng, n_random=80):
    """Scanne les fibres sur une geometrie CY donnee."""
    m = len(geom.ambient_dims)
    results = []

    for rank in [3, 4, 5]:
        gauge = GAUGE_GROUP_TABLE.get(rank, {}).get("group", "?")

        # Fibres structures
        bundles = []
        for s in range(min(m, 3)):
            ch = []
            for r in range(rank - 1):
                q = [0] * m; q[(r + s) % m] = 1; ch.append(q)
            last = [-sum(q[i] for q in ch) for i in range(m)]
            if all(abs(l) <= 3 for l in last):
                ch.append(last); bundles.append(ch)

        for i in range(min(m, 3)):
            for j in range(i + 1, min(m, 4)):
                ch = []
                for r in range(rank - 1):
                    q = [0] * m
                    q[i] = 1 if r % 2 == 0 else -1
                    q[j] = -1 if r % 2 == 0 else 1
                    ch.append(q)
                last = [-sum(q[k] for q in ch) for k in range(m)]
                ch.append(last); bundles.append(ch)

        # Fibres aleatoires
        for _ in range(n_random):
            ch = []
            for r in range(rank - 1):
                q = [int(rng.randint(-2, 3)) for _ in range(m)]
                ch.append(q)
            last = [-sum(q[k] for q in ch) for k in range(m)]
            if all(abs(l) <= 4 for l in last):
                ch.append(last); bundles.append(ch)

        for charges in bundles:
            b = LineBundleSum(charges=charges)
            if not b.c1_vanishes:
                continue
            try:
                cohom = bundle_cohomology_exact(geom.ambient_dims, geom.config_matrix, charges)
            except Exception:
                continue

            if gauge == "SU(5)": sp = extract_spectrum_su5(cohom)
            elif gauge == "SO(10)": sp = extract_spectrum_so10(cohom)
            else: sp = extract_spectrum_e6(cohom)

            if not sp.generation_match:
                continue

            stab = check_semi_stability(charges, geom.intersection_numbers)
            score = sp.sm_compatibility * 0.5 + 30
            if stab["semi_stable"]: score += 20

            results.append({
                'cicy': c['num'], 'h11': c['h11'], 'h21': c['h21'], 'chi': c['chi'],
                'ambient': "x".join(f"P{n}" for n in c['ambient']),
                'rank': rank, 'gauge': gauge,
                'n_gen': sp.n_generations, 'higgs': sp.n_higgs_candidates,
                'exotics': sp.n_exotics, 'stable': stab['semi_stable'],
                'score': round(score, 1),
                'reps': {k: v for k, v in sp.representations.items() if v > 0},
                'charges': [list(q) for q in charges],
            })

    return results


def main():
    parser = argparse.ArgumentParser(description="CY Landscape Explorer - Pipeline Oxford")
    parser.add_argument('file', nargs='?', default=None,
                       help="Chemin vers cicylist.txt (optionnel, sinon donnees embarquees)")
    parser.add_argument('--max-ps', type=int, default=None,
                       help="Nombre max d'espaces projectifs (ex: 6 pour accelerer)")
    parser.add_argument('--chi', type=int, default=None,
                       help="Filtrer par chi (ex: -4, -8)")
    parser.add_argument('--n-random', type=int, default=80,
                       help="Nombre de fibres aleatoires par rang (defaut: 80)")
    parser.add_argument('--output', type=str, default='output_oxford',
                       help="Dossier de sortie")
    args = parser.parse_args()

    out = args.output
    os.makedirs(out, exist_ok=True)
    rng = np.random.RandomState(42)

    print(f"\n{'='*66}")
    print(f"  CY LANDSCAPE EXPLORER — Donnees Oxford (Candelas et al.)")
    print(f"{'='*66}")

    # Phase 1
    print(f"\n[Phase 1] Chargement des CICYs")
    entries = load_cicys(args)
    if not entries:
        print("  Aucune CICY a traiter."); return

    # Phase 2
    print(f"\n[Phase 2] Validation (intersection numbers, chi)")
    valid, n_fail = validate_cicys(entries)
    print(f"  {len(valid)} valides, {n_fail} echouees")

    # Phase 3
    print(f"\n[Phase 3] Scan de fibres (cohomologie Koszul exacte)")
    all_results = []
    t0 = time.time()
    for i, (c, geom) in enumerate(valid):
        res = scan_bundles(geom, c, rng, args.n_random)
        if res:
            print(f"  CICY #{c['num']:>4} ({c['h11']},{c['h21']}) chi={c['chi']:>4}: "
                  f"{len(res)} fibres a 3 gen")
            all_results.extend(res)
        if (i + 1) % 50 == 0:
            dt = time.time() - t0
            print(f"    ... {i+1}/{len(valid)} CICYs traitees ({dt:.0f}s)")

    dt = time.time() - t0
    all_results.sort(key=lambda r: r['score'], reverse=True)

    # Phase 4
    print(f"\n{'='*66}")
    print(f"  RESULTATS")
    print(f"{'='*66}")
    print(f"  CICYs analysees  : {len(valid)}")
    print(f"  Fibres a 3 gen   : {len(all_results)}")
    n_stable = sum(1 for r in all_results if r['stable'])
    print(f"  Semi-stables     : {n_stable}")
    print(f"  Temps total      : {dt:.1f}s")

    if all_results:
        top = all_results[:20]
        print(f"\n  {'#':>2} {'CICY':>5} {'Ambient':<20} {'CY':>8} {'Rk':>3} {'Jauge':>7} "
              f"{'Hig':>3} {'Exo':>3} {'Stab':>4} {'Score':>5}")
        print(f"  {'='*2} {'='*5} {'='*20} {'='*8} {'='*3} {'='*7} {'='*3} {'='*3} {'='*4} {'='*5}")
        for i, r in enumerate(top):
            s = "S" if r['stable'] else "U"
            print(f"  {i+1:>2} #{r['cicy']:>4} {r['ambient']:<20} ({r['h11']:>2},{r['h21']:>2}) "
                  f"SU{r['rank']} {r['gauge']:>7} {r['higgs']:>3} {r['exotics']:>3} "
                  f"{s:>4} {r['score']:>5}")

        # Export JSON
        class NpEncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, (np.integer,)): return int(obj)
                if isinstance(obj, (np.floating,)): return float(obj)
                if isinstance(obj, np.ndarray): return obj.tolist()
                return super().default(obj)

        export = {'source': 'Oxford CICY list (Candelas et al.)',
                  'n_cicys': len(valid), 'n_3gen': len(all_results),
                  'n_stable': n_stable, 'results': all_results[:100]}
        with open(os.path.join(out, 'oxford_results.json'), 'w') as f:
            json.dump(export, f, indent=2, cls=NpEncoder)
        print(f"\n  Export: {out}/oxford_results.json")


if __name__ == "__main__":
    main()
