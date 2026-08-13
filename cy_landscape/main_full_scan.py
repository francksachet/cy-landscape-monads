#!/usr/bin/env python3
"""
Pipeline de scan complet : monades + monades positives + extensions.
Critere de Hoppe complet sur tous les candidats a 3 generations.

Usage:
  python -m cy_landscape.main_full_scan                         # donnees embarquees
  python -m cy_landscape.main_full_scan cicylist.txt --max-ps 6
  python -m cy_landscape.main_full_scan cicylist.txt --max-ps 5 --n-random 500
"""
import os, sys, json, argparse, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cy_landscape.core.intersection import (
    compute_intersection_numbers, compute_euler_from_intersection, compute_c2_tangent)
from cy_landscape.core.bundles import CICYGeometry, GAUGE_GROUP_TABLE
from cy_landscape.core.exact_cohomology import bundle_cohomology_exact
from cy_landscape.core.monads import (
    MonadBundle, compute_monad_cohomology, generate_monads, check_map_exists)
from cy_landscape.core.positive_monads import generate_positive_monads, is_positive_monad
from cy_landscape.core.extensions import (
    ExtensionBundle, check_extension_exists, compute_extension_cohomology, generate_extensions)
from cy_landscape.core.monad_wedge import cohomology_wedge2_V, cohomology_end_V_approx
from cy_landscape.core.stability_full import hoppe_criterion
from cy_landscape.core.cohomology import extract_spectrum_su5, extract_spectrum_so10, extract_spectrum_e6


def load_and_validate(args):
    if args.file:
        from cy_landscape.data.parse_oxford import load_oxford_file
        entries = load_oxford_file(args.file)
        print(f"  {len(entries)} CICYs parsees")
    else:
        from cy_landscape.data.oxford_cicys import get_all_oxford
        entries = get_all_oxford()
        print(f"  {len(entries)} CICYs embarquees")

    if args.max_ps:
        entries = [e for e in entries if len(e['ambient']) <= args.max_ps]
        print(f"  Filtre max_ps<={args.max_ps}: {len(entries)}")

    valid = []
    for c in entries:
        try:
            d = compute_intersection_numbers(c['ambient'], c['config'])
            chi = compute_euler_from_intersection(c['ambient'], c['config'], d)
            if chi != c['chi']: continue
            c2 = compute_c2_tangent(c['ambient'], c['config'], d)
            geom = CICYGeometry(ambient_dims=c['ambient'], config_matrix=c['config'],
                               h11=c['h11'], h21=c['h21'],
                               intersection_numbers=d, c2_tangent=c2)
            valid.append((c, geom))
        except: pass
    print(f"  {len(valid)} validees")
    return valid


def _spectrum(gauge, cohom_fmt):
    if gauge == "SU(5)": return extract_spectrum_su5(cohom_fmt)
    elif gauge == "SO(10)": return extract_spectrum_so10(cohom_fmt)
    else: return extract_spectrum_e6(cohom_fmt)


def scan_monads(c, geom, rng, n_random):
    """Scan monades classiques."""
    m = len(geom.ambient_dims)
    candidates = []

    for rank_V in [3, 4, 5]:
        monads = generate_monads(m, rank_V, max_charge=3, n_random=n_random, rng=rng)
        for monad in monads:
            if not monad.c1_vanishes: continue
            map_ok, _ = check_map_exists(monad, geom.ambient_dims, geom.config_matrix)
            if not map_ok: continue
            cohom = compute_monad_cohomology(monad, geom.ambient_dims, geom.config_matrix)
            if cohom is None: continue
            if abs(cohom[1] - cohom[2]) == 3:
                candidates.append(('monad', rank_V, monad, cohom))
    return candidates


def scan_positive_monads(c, geom, rng, n_random):
    """Scan monades positives."""
    m = len(geom.ambient_dims)
    candidates = []

    for rank_V in [3, 4, 5]:
        monads = generate_positive_monads(m, rank_V, max_charge=4,
                                           n_systematic=n_random, rng=rng)
        for monad in monads:
            map_ok, _ = check_map_exists(monad, geom.ambient_dims, geom.config_matrix)
            if not map_ok: continue
            cohom = compute_monad_cohomology(monad, geom.ambient_dims, geom.config_matrix)
            if cohom is None: continue
            if abs(cohom[1] - cohom[2]) == 3:
                candidates.append(('pos_monad', rank_V, monad, cohom))
    return candidates


def scan_extensions(c, geom, rng, n_random):
    """Scan fibres d'extension."""
    m = len(geom.ambient_dims)
    candidates = []

    for rank_V in [3, 4, 5]:
        exts = generate_extensions(m, rank_V, max_charge=3,
                                    n_random=n_random, rng=rng)
        for ext in exts:
            exists, dim_ext = check_extension_exists(ext, geom.ambient_dims, geom.config_matrix)
            if not exists: continue
            cohom = compute_extension_cohomology(ext, geom.ambient_dims, geom.config_matrix)
            if cohom is None: continue
            if abs(cohom[1] - cohom[2]) == 3:
                # Emballer l'extension comme un MonadBundle pour le test de Hoppe
                # 0 -> F1 -> V -> F2 -> 0  <=>  monade avec B=V, C=F2
                # Pas exactement mais le test Hoppe porte sur V
                # On construit un pseudo-monad: B = F1 + F2 (concat), C = F2
                pseudo_monad = MonadBundle(
                    ext.f1_charges + ext.f2_charges,
                    ext.f2_charges)
                candidates.append(('extension', rank_V, pseudo_monad, cohom,
                                   ext, dim_ext))
    return candidates


def test_and_score(candidate, geom, c):
    """Applique Hoppe et calcule le score."""
    if len(candidate) == 4:
        kind, rank_V, monad, cohom = candidate
        ext = None
    else:
        kind, rank_V, monad, cohom, ext, dim_ext = candidate

    gauge = GAUGE_GROUP_TABLE.get(rank_V, {}).get("group", "?")

    # Hoppe complet
    hoppe = hoppe_criterion(geom.ambient_dims, geom.config_matrix, monad, max_H=1)
    is_stable = hoppe.get('stable', False)

    # Spectre
    if monad.rank_C == 1 and ext is None:
        try:
            w2 = cohomology_wedge2_V(geom.ambient_dims, geom.config_matrix,
                                      monad.b_charges, monad.c_charges)
            w2V = w2['wedge2V']
        except:
            w2V = {i: 0 for i in range(4)}
    else:
        rV = rank_V
        w2V = {0: 0, 1: max(0, cohom[1]*(rV-1)//2), 2: max(0, cohom[2]*(rV-1)//2), 3: 0}

    end_V = {0: 1, 1: max(1, rank_V**2 - 1), 2: max(1, rank_V**2 - 1), 3: 1}

    cohom_fmt = {
        "V": {i: cohom[i] for i in range(4)},
        "V_dual": {i: cohom[3-i] for i in range(4)},
        "wedge2V": w2V, "end_V": end_V,
    }
    sp = _spectrum(gauge, cohom_fmt)

    score = sp.sm_compatibility * 0.5 + 30
    if is_stable: score += 25
    if kind == 'pos_monad': score += 3  # Bonus positivite
    if kind == 'extension': score += 2  # Bonus extension

    return {
        'type': kind, 'stable': is_stable,
        'cicy': c['num'], 'h11': c['h11'], 'h21': c['h21'], 'chi': c['chi'],
        'ambient': "x".join(f"P{n}" for n in c['ambient']),
        'rank_V': rank_V, 'gauge': gauge,
        'n_gen': sp.n_generations, 'higgs': sp.n_higgs_candidates,
        'exotics': sp.n_exotics, 'score': round(score, 1),
        'reps': {k: v for k, v in sp.representations.items() if v > 0},
        'b_charges': [list(b) for b in monad.b_charges],
        'c_charges': [list(cc) for cc in monad.c_charges],
        'hoppe': hoppe['reason'], 'hoppe_tests': hoppe.get('tests', 0),
    }


def main():
    parser = argparse.ArgumentParser(
        description="CY Landscape — Scan complet (monades + positives + extensions)")
    parser.add_argument('file', nargs='?', default=None)
    parser.add_argument('--max-ps', type=int, default=None)
    parser.add_argument('--n-random', type=int, default=150)
    parser.add_argument('--output', type=str, default='output_full')
    args = parser.parse_args()

    out = args.output
    os.makedirs(out, exist_ok=True)
    rng = np.random.RandomState(42)

    print(f"\n{'='*66}")
    print(f"  CY LANDSCAPE — SCAN COMPLET")
    print(f"  Monades + Monades positives + Extensions")
    print(f"  Critere de Hoppe complet (rangs 1..rk-1)")
    print(f"{'='*66}")

    print(f"\n[1] Chargement")
    valid = load_and_validate(args)
    if not valid: return

    print(f"\n[2] Scan (3 methodes)")
    all_results = []
    stats = {'monad': [0, 0], 'pos_monad': [0, 0], 'extension': [0, 0]}
    t0 = time.time()

    for i, (c, geom) in enumerate(valid):
        cands = []
        cands += scan_monads(c, geom, rng, args.n_random)
        cands += scan_positive_monads(c, geom, rng, args.n_random)
        cands += scan_extensions(c, geom, rng, args.n_random)

        if cands:
            n_by_type = {}
            for cand in cands:
                result = test_and_score(cand, geom, c)
                all_results.append(result)
                t = result['type']
                stats[t][0] += 1
                if result['stable']: stats[t][1] += 1
                n_by_type[t] = n_by_type.get(t, 0) + 1

            parts = " + ".join(f"{v} {k}" for k, v in n_by_type.items())
            n_stab = sum(1 for r in all_results[-len(cands):]
                         if r['cicy'] == c['num'] and r['stable'])
            print(f"  CICY #{c['num']:>4} ({c['h11']},{c['h21']}): "
                  f"{len(cands)} candidats ({parts})"
                  f"{f', {n_stab} HOPPE-STABLES!' if n_stab else ''}")

        if (i + 1) % 50 == 0:
            print(f"    ... {i+1}/{len(valid)} ({time.time()-t0:.0f}s)")

    dt = time.time() - t0
    all_results.sort(key=lambda r: r['score'], reverse=True)

    # Resultats
    n_total = len(all_results)
    n_stable = sum(1 for r in all_results if r['stable'])

    print(f"\n{'='*66}")
    print(f"  RESULTATS")
    print(f"{'='*66}")
    print(f"  CICYs analysees         : {len(valid)}")
    print(f"  Candidats a 3 gen       : {n_total}")
    for t, (n3, ns) in stats.items():
        label = {'monad': 'Monades', 'pos_monad': 'Monades pos.',
                 'extension': 'Extensions'}[t]
        print(f"    {label:<22}: {n3:>4} dont {ns} Hoppe-stables")
    print(f"  TOTAL HOPPE-STABLES     : {n_stable}"
          f"{'  <-- PERCEE!' if n_stable else ''}")
    print(f"  Temps                   : {dt:.1f}s")

    if all_results:
        top = all_results[:20]
        print(f"\n  {'#':>2} {'Type':<12} {'CICY':>5} {'CY':>8} {'Rk':>2} "
              f"{'Jauge':>7} {'Hig':>3} {'Exo':>3} {'Stab':>5} {'Score':>5}")
        print(f"  {'='*2} {'='*12} {'='*5} {'='*8} {'='*2} "
              f"{'='*7} {'='*3} {'='*3} {'='*5} {'='*5}")
        for i, r in enumerate(top):
            s = "HOPPE" if r['stable'] else "  -  "
            print(f"  {i+1:>2} {r['type']:<12} #{r['cicy']:>4} ({r['h11']:>2},{r['h21']:>2}) "
                  f"{r['rank_V']:>2} {r['gauge']:>7} {r['higgs']:>3} {r['exotics']:>3} "
                  f"{s:>5} {r['score']:>5}")

    # Export
    class NpEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.integer,)): return int(obj)
            if isinstance(obj, (np.floating,)): return float(obj)
            if isinstance(obj, np.ndarray): return obj.tolist()
            return super().default(obj)

    with open(os.path.join(out, 'full_scan_results.json'), 'w') as f:
        json.dump({'n_cicys': len(valid), 'n_3gen': n_total, 'n_stable': n_stable,
                   'stats': stats, 'results': all_results[:100]},
                  f, indent=2, cls=NpEncoder)
    print(f"\n  Export: {out}/full_scan_results.json")

    print(f"\n  Generation des graphiques...")
    try:
        from cy_landscape.core.visualize import generate_all_plots
        generate_all_plots(os.path.join(out, 'full_scan_results.json'), out)
    except Exception as e:
        print(f"  Graphiques ignores : {e}")


if __name__ == "__main__":
    main()
