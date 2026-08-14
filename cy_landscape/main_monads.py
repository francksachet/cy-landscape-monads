#!/usr/bin/env python3
"""
Pipeline monade : scan de fibres monades sur les CICYs d'Oxford.

Usage:
  python -m cy_landscape.main_monads                         # donnees embarquees
  python -m cy_landscape.main_monads cicylist.txt            # fichier Oxford
  python -m cy_landscape.main_monads cicylist.txt --max-ps 5 --n-random 200
"""
import os, sys, json, argparse, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cy_landscape.core.intersection import (
    compute_intersection_numbers, compute_euler_from_intersection, compute_c2_tangent)
from cy_landscape.core.bundles import CICYGeometry, GAUGE_GROUP_TABLE
from cy_landscape.core.cohomology import extract_spectrum_su5, extract_spectrum_so10, extract_spectrum_e6
from cy_landscape.core.monad_wedge import cohomology_wedge2_V, cohomology_end_V_approx
from cy_landscape.core.stability_full import hoppe_criterion
from cy_landscape.core.monads import (
    MonadBundle, check_map_exists, compute_monad_cohomology,
    check_monad_stability, generate_monads)


def load_and_validate(args):
    """Charge et valide les CICYs."""
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
            c2 = compute_c2_tangent(c['ambient'], c['config'], d)
            if chi != c['chi']: continue
            geom = CICYGeometry(
                ambient_dims=c['ambient'], config_matrix=c['config'],
                h11=c['h11'], h21=c['h21'],
                intersection_numbers=d, c2_tangent=c2)
            valid.append((c, geom))
        except: pass
    print(f"  {len(valid)} validees")
    return valid


def scan_monads_on_cicy(c, geom, rng, n_random=100, verbose=False):
    """Scanne les fibres monades sur une CICY donnee."""
    m = len(geom.ambient_dims)
    results = []

    for rank_V in [3, 4, 5]:
        gauge = GAUGE_GROUP_TABLE.get(rank_V, {}).get("group", "?")
        monads = generate_monads(m, rank_V, max_charge=3, n_random=n_random, rng=rng)

        for monad in monads:
            if not monad.c1_vanishes:
                continue

            # Verifier que la map existe
            map_ok, _ = check_map_exists(monad, geom.ambient_dims, geom.config_matrix)
            if not map_ok:
                continue

            # Cohomologie
            cohom = compute_monad_cohomology(monad, geom.ambient_dims, geom.config_matrix)
            if cohom is None:
                continue

            # Cohomologie EXACTE de wedge2V (si rk(C) = 1)
            if monad.rank_C == 1:
                try:
                    w2_res = cohomology_wedge2_V(
                        geom.ambient_dims, geom.config_matrix,
                        monad.b_charges, monad.c_charges)
                    w2V = w2_res['wedge2V']
                except Exception:
                    rV = monad.rank_V
                    w2V = {0:0, 1: max(0, cohom[1]*(rV-1)//2), 2: max(0, cohom[2]*(rV-1)//2), 3:0}
            else:
                rV = monad.rank_V
                w2V = {0:0, 1: max(0, cohom[1]*(rV-1)//2), 2: max(0, cohom[2]*(rV-1)//2), 3:0}
            
            hV = {i: cohom[i] for i in range(4)}
            end_V = cohomology_end_V_approx(
                geom.ambient_dims, geom.config_matrix,
                monad.b_charges, monad.c_charges, hV)
            
            cohom_fmt = {
                "V": {i: cohom[i] for i in range(4)},
                "V_dual": {i: cohom[3-i] for i in range(4)},
                "wedge2V": w2V,
                "end_V": end_V,
            }
            if gauge == "SU(5)": sp = extract_spectrum_su5(cohom_fmt)
            elif gauge == "SO(10)": sp = extract_spectrum_so10(cohom_fmt)
            else: sp = extract_spectrum_e6(cohom_fmt)

            if not sp.generation_match:
                continue

            # Stabilite (test monade)
            # Test de stabilite: Hoppe complet (rangs 1..rk-1)
            stab = hoppe_criterion(geom.ambient_dims, geom.config_matrix,
                                   monad, max_H=1)

            score = sp.sm_compatibility * 0.5 + 30
            if stab.get('stable', stab.get('semi_stable', False)): score += 20
            label = "MONADE"

            results.append({
                'type': 'monad',
                'cicy': c['num'], 'h11': c['h11'], 'h21': c['h21'], 'chi': c['chi'],
                'ambient': "x".join(f"P{n}" for n in c['ambient']),
                'rank_V': rank_V, 'rank_B': monad.rank_B, 'rank_C': monad.rank_C,
                'gauge': gauge,
                'n_gen': sp.n_generations, 'higgs': sp.n_higgs_candidates,
                'exotics': sp.n_exotics, 'stable': stab.get('stable', stab.get('semi_stable', False)),
                'stability_reason': stab['reason'],
                'score': round(score, 1),
                'reps': {k: v for k, v in sp.representations.items() if v > 0},
                'b_charges': [list(b) for b in monad.b_charges],
                'c_charges': [list(c_) for c_ in monad.c_charges],
                'cohomology': {k: v for k, v in cohom.items() if isinstance(k, int)},
            })

    return results


def main():
    parser = argparse.ArgumentParser(description="CY Landscape — Scan de fibres monades")
    parser.add_argument('file', nargs='?', default=None)
    parser.add_argument('--max-ps', type=int, default=None)
    parser.add_argument('--n-random', type=int, default=100)
    parser.add_argument('--output', type=str, default='output_monads')
    args = parser.parse_args()

    out = args.output
    os.makedirs(out, exist_ok=True)
    rng = np.random.RandomState(42)

    print(f"\n{'='*66}")
    print(f"  CY LANDSCAPE — FIBRES MONADES")
    print(f"  0 -> V -> B -> C -> 0")
    print(f"{'='*66}")

    print(f"\n[1] Chargement")
    valid = load_and_validate(args)
    if not valid: return

    print(f"\n[2] Scan de monades")
    all_results = []
    t0 = time.time()

    for i, (c, geom) in enumerate(valid):
        res = scan_monads_on_cicy(c, geom, rng, args.n_random)
        if res:
            n_stab = sum(1 for r in res if r['stable'])
            print(f"  CICY #{c['num']:>4} ({c['h11']},{c['h21']}): "
                  f"{len(res)} monades 3-gen, {n_stab} stables")
            all_results.extend(res)
        if (i + 1) % 20 == 0:
            print(f"    ... {i+1}/{len(valid)} ({time.time()-t0:.0f}s)")

    dt = time.time() - t0
    all_results.sort(key=lambda r: r['score'], reverse=True)

    # Resultats
    n_total = len(all_results)
    n_stable = sum(1 for r in all_results if r['stable'])

    print(f"\n{'='*66}")
    print(f"  RESULTATS — FIBRES MONADES")
    print(f"{'='*66}")
    print(f"  CICYs analysees     : {len(valid)}")
    print(f"  Monades a 3 gen     : {n_total}")
    print(f"  SEMI-STABLES        : {n_stable}  {'<-- PERCEE!' if n_stable > 0 else ''}")
    print(f"  Temps               : {dt:.1f}s")

    if all_results:
        top = all_results[:20]
        print(f"\n  {'#':>2} {'CICY':>5} {'Ambient':<18} {'CY':>8} {'V':>2}={'B':>2}-{'C':>2} "
              f"{'Jauge':>7} {'Hig':>3} {'Exo':>3} {'Stab':>4} {'Score':>5}")
        print(f"  {'='*2} {'='*5} {'='*18} {'='*8} {'='*2} {'='*2} {'='*2} "
              f"{'='*7} {'='*3} {'='*3} {'='*4} {'='*5}")
        for i, r in enumerate(top):
            s = "S" if r['stable'] else "U"
            print(f"  {i+1:>2} #{r['cicy']:>4} {r['ambient']:<18} ({r['h11']:>2},{r['h21']:>2}) "
                  f"{r['rank_V']:>2}={r['rank_B']:>2}-{r['rank_C']:>2} "
                  f"{r['gauge']:>7} {r['higgs']:>3} {r['exotics']:>3} {s:>4} {r['score']:>5}")

        # Detail du meilleur stable
        best_stable = next((r for r in all_results if r['stable']), None)
        if best_stable:
            r = best_stable
            print(f"\n  === MEILLEUR FIBRE MONADE SEMI-STABLE ===")
            print(f"  CICY #{r['cicy']} ({r['ambient']}), h=({r['h11']},{r['h21']})")
            print(f"  Monade : 0 -> V(rk {r['rank_V']}) -> B(rk {r['rank_B']}) -> C(rk {r['rank_C']}) -> 0")
            print(f"  B = {r['b_charges']}")
            print(f"  C = {r['c_charges']}")
            print(f"  Groupe de jauge : {r['gauge']}")
            print(f"  Spectre : {r['reps']}")
            print(f"  H^i(V) = {r['cohomology']}")
            print(f"  Semi-stable : oui ({r['stability_reason']})")
            print(f"  Score : {r['score']}/100")

    # Export
    class NpEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, (np.integer,)): return int(obj)
            if isinstance(obj, (np.floating,)): return float(obj)
            if isinstance(obj, np.ndarray): return obj.tolist()
            return super().default(obj)

    with open(os.path.join(out, 'monad_results.json'), 'w') as f:
        json.dump({'n_cicys': len(valid), 'n_3gen': n_total, 'n_stable': n_stable,
                   'results': all_results[:100]}, f, indent=2, cls=NpEncoder)
    print(f"\n  Export: {out}/monad_results.json")



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
