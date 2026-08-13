"""Saturation du generateur : combien de monades DISTINCTES par graine ?

On ne mesure que la generation et le prefiltre chi -- pas la stabilite, qui
coute l'essentiel du scan mais ne change pas la question : si le generateur
sature, tout ce qui en decoule sature aussi.
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from cy_landscape.core.positive_monads import generate_positive_monads
from cy_landscape.core.intersection import (compute_intersection_numbers,
    compute_euler_from_intersection, compute_c2_tangent)
from cy_landscape.core.chi_exact import ChiCalculator
from cy_landscape.data.parse_oxford import load_oxford_file

GRAINE = int(sys.argv[1]); MAXCH = 5; NRAND = 8000
DEB = int(sys.argv[2]) if len(sys.argv)>2 else 0
FIN = int(sys.argv[3]) if len(sys.argv)>3 else 10**6
BUDGET = float(sys.argv[4]) if len(sys.argv)>4 else 140.0
DST = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'monades.jsonl')

w = json.load(open('wilson_gros_gamma.json'))
entries = {e['num']: e for e in load_oxford_file('cicylist.txt')}
cibles = {int(k): {3*o for o in v['ordres']} for k, v in w.items()}

t0 = time.time(); n_tot = 0; n_chi = 0
with open(DST, 'a', encoding='utf-8') as fh:
    for idx, num in enumerate(sorted(cibles)):
        if idx < DEB or idx >= FIN: continue
        if time.time()-t0 > BUDGET: print(f'  arret budget a idx={idx}'); break
        e = entries.get(num)
        if e is None: continue
        amb, cfg = e['ambient'], e['config']
        m = len(amb)
        try:
            d = compute_intersection_numbers(amb, cfg)
            c2 = compute_c2_tangent(amb, cfg, d)
            chical = ChiCalculator(amb, d, c2)
        except Exception:
            continue
        for rank_V in (3, 4, 5):
            rng = np.random.RandomState(GRAINE)
            for mo in generate_positive_monads(m, rank_V, max_charge=MAXCH,
                                               n_systematic=NRAND, rng=rng,
                                               sampling_threshold=None,
                                               sampling_budget=100, seed=GRAINE):
                if not mo.c1_vanishes: continue
                n_tot += 1
                chi = chical.bundle(mo.b_charges) - chical.bundle(mo.c_charges)
                if abs(chi) not in cibles[num]: continue
                n_chi += 1
                fh.write(json.dumps({'g': GRAINE, 'cicy': num,
                                     'b': [list(x) for x in mo.b_charges],
                                     'c': [list(x) for x in mo.c_charges]}) + '\n')
print(f"graine {GRAINE} [{DEB}:{FIN}] : {n_tot} monades engendrees, {n_chi} passent le prefiltre chi"
      f"  ({time.time()-t0:.0f}s)")
