#!/usr/bin/env python3
"""
verify_hoppe.py -- Test de Hoppe manquant (H = 0), sur les seuls survivants.

--------------------------------------------------------------------------
Ce qui est teste ici et ne l'a pas ete pendant le scan
--------------------------------------------------------------------------
Le critere de Hoppe, pour V de premiere classe de Chern nulle, s'enonce :

    V est stable  <=>  h0(wedge^p V) = 0  pour tout p = 1..rk(V)-1

`hoppe_fast` teste h0(wedge^p V(-H)) pour H = e_i seulement, jamais pour
H = 0. Pour H ample, H0(wedge^p V(-H)) est un SOUS-espace de H0(wedge^p V) :
le test effectue pendant le scan est donc strictement plus faible que le
critere. Il ne peut pas rejeter davantage, seulement moins.

Aucune CICY n'est rescannee : chaque fibre est reconstruit a partir des
charges B et C stockees dans le JSONL.

--------------------------------------------------------------------------
Reduction des tests par isomorphismes (c1(V) = 0, donc det V = O)
--------------------------------------------------------------------------
wedge^(rk-1) V est isomorphe a V*, et wedge^(rk-2) V a wedge^2 V*. Avec la
dualite de Serre sur un CY3 (K trivial), h0(F*) = h3(F). D'ou :

  rk = 3 : p=1 -> h0(V)          p=2 -> h3(V)
  rk = 4 : p=1 -> h0(V)          p=2 -> h0(w2V)     p=3 -> h3(V)
  rk = 5 : p=1 -> h0(V)          p=2 -> h0(w2V)
           p=3 -> h3(w2V)        p=4 -> h3(V)

Consequence notable : en rang 3, le critere se reduit entierement a
h0(V) = h3(V) = 0. Les candidats E6 deja filtres sur ces deux conditions
passent donc automatiquement -- il n'y a rien de neuf a apprendre sur eux
par ce chemin.

ATTENTION : `_wedge3_h0_twisted` de `hoppe_fast.py` implemente h0(w3V) =
h3(V), identite valable en rang 4 seulement (sa docstring le dit), mais
`hoppe_fast` l'appelle des que rk >= 4, donc aussi en rang 5 ou elle est
fausse. Ce script ne l'utilise pas et applique le tableau ci-dessus.

--------------------------------------------------------------------------
Verdicts
--------------------------------------------------------------------------
  STABLE      tous les tests passent
  NON STABLE  au moins un h0 non nul
  PARTIEL     rank_C >= 2 : wedge^2 V n'est pas calculable par ce chemin
              (monad_wedge n'exp1oite la suite exacte que pour rank_C = 1).
              Ne pas compter ces fibres comme stables.

Un controle croise gratuit est affiche en fin de rapport : en rang 3,
h0(w2V) calcule par monad_wedge DOIT valoir h3(V). Tout ecart signale un
defaut de `monad_wedge.py` lui-meme.

Usage:
    python verify_hoppe.py output_optimized cicylist.txt
    python verify_hoppe.py output_optimized cicylist.txt --input results_ranked.jsonl
"""
import os, sys, json, argparse
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def _sortie_tolerante():
    """
    Empeche un plantage d'encodage sur une console Windows.

    Les etiquettes de groupe de jauge contiennent des indices Unicode
    ("E\u2086"), que la console cp1252 ne sait pas encoder. `print` levait
    alors UnicodeEncodeError, le script mourait sans rien ecrire, et la
    commande suivante de l'enchainement tournait sur un fichier absent.
    Cas reel : scan_exh2, ou `equivariance.py` est mort ainsi et a fait
    echouer `equivariance_f.py` derriere lui.
    """
    import sys as _sys
    for flux in (_sys.stdout, _sys.stderr):
        try:
            flux.reconfigure(errors='replace')
        except (AttributeError, ValueError):
            pass


def main():
    _sortie_tolerante()
    ap = argparse.ArgumentParser()
    ap.add_argument('output_dir')
    ap.add_argument('cicyfile', nargs='?', default=None,
                    help="cicylist.txt (sinon : CICYs embarquees)")
    ap.add_argument('--input', default='results_clean.jsonl')
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--verbose', action='store_true',
                    help="Afficher aussi les fibres rejetes")
    args = ap.parse_args()

    from cy_landscape.core.monads import MonadBundle, compute_monad_cohomology
    from cy_landscape.core.cache import set_geometry
    from cy_landscape.core.hoppe_fast import _monad_h0_twisted
    from cy_landscape.core.monad_wedge import cohomology_wedge2_V

    if args.cicyfile:
        from cy_landscape.data.parse_oxford import load_oxford_file
        entries = load_oxford_file(args.cicyfile)
    else:
        from cy_landscape.data.oxford_cicys import get_all_oxford
        entries = get_all_oxford()
    by_num = {e['num']: e for e in entries}

    src = os.path.join(args.output_dir, args.input)
    if not os.path.exists(src):
        print(f"Introuvable : {src}")
        return 1
    rs = []
    with open(src, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line:
                rs.append(json.loads(line))
    if args.limit:
        rs = rs[:args.limit]

    print(f"\n{'='*80}")
    print(f"  TEST DE HOPPE MANQUANT (H = 0) -- {len(rs)} entrees")
    print(f"{'='*80}")

    verdicts = Counter()
    out_rows = []
    coherence = []          # (cicy, h0w2, h3V) pour les rang 3
    shown = 0

    for r in rs:
        c = by_num.get(r['cicy'])
        if c is None:
            verdicts['cicy_absente_du_fichier'] += 1
            continue
        if r.get('type') == 'extension':
            verdicts['extension_ignoree'] += 1
            continue

        b = [list(x) for x in r['b_charges']]
        cc = [list(x) for x in r['c_charges']]
        monad = MonadBundle(b, cc)
        rk = monad.rank_V
        m = len(c['ambient'])
        H0 = [0] * m

        try:
            set_geometry(c['ambient'], c['config'])
            cohV = compute_monad_cohomology(monad, c['ambient'], c['config'])
            if cohV is None:
                verdicts['cohomologie_indisponible'] += 1
                continue
            h0V = _monad_h0_twisted(c['ambient'], c['config'], monad, H0)
            h3V = int(cohV[3])

            w2 = None
            if monad.rank_C == 1:
                w2 = cohomology_wedge2_V(c['ambient'], c['config'], b, cc)['wedge2V']
        except Exception as e:
            verdicts[f'erreur:{type(e).__name__}'] += 1
            continue

        tests = {'p1_h0V': h0V}
        complet = True

        if rk == 3:
            tests['p2_h3V'] = h3V
            if w2 is not None:
                coherence.append((r['cicy'], int(w2[0]), h3V))
        elif rk == 4:
            if w2 is None:
                complet = False
            else:
                tests['p2_h0w2V'] = int(w2[0])
            tests['p3_h3V'] = h3V
        elif rk == 5:
            if w2 is None:
                complet = False
            else:
                tests['p2_h0w2V'] = int(w2[0])
                tests['p3_h3w2V'] = int(w2[3])
            tests['p4_h3V'] = h3V
        else:
            complet = False

        echecs = {k: v for k, v in tests.items() if v > 0}
        if echecs:
            verdict = 'NON STABLE'
        elif complet:
            verdict = 'STABLE'
        else:
            verdict = 'PARTIEL (rank_C>=2)'
        verdicts[verdict] += 1

        r2 = dict(r)
        r2['hoppe_H0'] = {'tests': tests, 'verdict': verdict}
        if w2 is not None:
            r2['wedge2V_recalcule'] = {str(k): int(v) for k, v in w2.items()}
        out_rows.append(r2)

        if verdict != 'NON STABLE' or args.verbose:
            if shown == 0:
                print(f"\n    {'#CICY':>6} {'type':<10} {'jauge':>7} {'rk':>2}  "
                      f"tests h0(wedge^p V)                      verdict")
            shown += 1
            detail = " ".join(f"{k}={v}" for k, v in tests.items())
            print(f"    {r['cicy']:>6} {r['type']:<10} {r['gauge']:>7} {rk:>2}  "
                  f"{detail:<40} {verdict}")

    print(f"\n{'='*80}")
    print(f"  BILAN")
    print(f"{'='*80}")
    for k, v in verdicts.most_common():
        print(f"    {k:<28} {v:>6}")

    stables = [r for r in out_rows if r['hoppe_H0']['verdict'] == 'STABLE']
    if stables:
        print(f"\n  Fibres reellement stables : {len(stables)}")
        print(f"    {'#CICY':>6} {'jauge':>7} {'rk':>2} {'H stocke':>9} "
              f"{'h1(w2V) recalc':>15} {'cohomologie':>16}")
        for r in sorted(stables, key=lambda x: x.get('higgs', 0)):
            w2 = r.get('wedge2V_recalcule') or {}
            print(f"    {r['cicy']:>6} {r['gauge']:>7} {r['rank_V']:>2} "
                  f"{r.get('higgs',0):>9} {str(w2.get('1','-')):>15} "
                  f"{str(r.get('cohomology')):>16}")

    if coherence:
        ecarts = [(n, a, b_) for n, a, b_ in coherence if a != b_]
        print(f"\n  Controle croise rang 3 : h0(w2V) doit valoir h3(V)")
        print(f"    verifie sur {len(coherence)} fibres, {len(ecarts)} ecart(s)")
        for n, a, b_ in ecarts[:10]:
            print(f"      CICY #{n} : h0(w2V)={a}  h3(V)={b_}")
        if ecarts:
            print(f"    -> ecart = defaut de monad_wedge.py, pas des candidats")

    out = os.path.join(args.output_dir, 'results_verified.jsonl')
    with open(out, 'w', encoding='utf-8') as f:
        for r in out_rows:
            f.write(json.dumps(r) + '\n')
    print(f"\n  Ecrit : {out}")
    print(f"{'='*80}\n")
    return 0


if __name__ == '__main__':
    sys.exit(main())
