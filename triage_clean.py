#!/usr/bin/env python3
"""
triage_clean.py -- Second tri, sur results_clean.jsonl.

Ne recalcule rien non plus. Recorrige a la lecture deux quantites que le
pipeline reporte de facon inutilisable, puis reclasse les candidats.

--------------------------------------------------------------------------
1. Le nombre d'exotiques stocke est structurellement nul
--------------------------------------------------------------------------
Dans `cohomology.py` :

  - SO(10) : `sp.n_exotics = 0` est ecrit en dur.
  - SU(5)  : `n_exotics = max(0, n_10 + n_10bar - n_gen - 2*n_anti)`.
             Or n_gen = |a-b| et n_anti = min(a,b), donc
             |a-b| + 2*min(a,b) = a+b identiquement : cette expression
             vaut 0 quelles que soient les valeurs. Elle ne peut pas
             etre non nulle.
  - E6     : seul cas ou n_exotics = n_anti est effectivement reporte.

Le champ "exotics = 0" des candidats SO(10) et SU(5) n'est donc pas un
resultat : c'est une constante. Il vaut aussi 25 points gratuits dans
`compute_sm_compatibility`, ce qui fausse le classement par score en
faveur de SO(10)/SU(5) contre E6.

--------------------------------------------------------------------------
2. La quantite physique pertinente : n_anti = min(h1, h2)
--------------------------------------------------------------------------
h1(V) et h2(V) = h1(V*) comptent les familles et les anti-familles. Leur
difference donne les 3 generations chirales ; leur MINIMUM compte les
paires vectorielles (16 + 16bar, ou 10 + 10bar, ou 27 + 27bar). Ces paires
ne sont pas protegees par la chiralite : elles peuvent en principe
s'apparier et devenir massives, mais rien ici ne le demontre, et un
modele avec n_anti = 0 est qualitativement meilleur qu'un modele avec
n_anti = 3.

Une cohomologie [0, 3, 6, 0] signifie donc 3 generations ET 3 paires
vectorielles -- pas un spectre propre a 3 familles.

Le seul spectre reellement propre est n_anti = 0, c'est-a-dire une
cohomologie [0, 3, 0, 0] ou [0, 0, 3, 0].

Usage:
    python triage_clean.py output_optimized
"""
import os, sys, json, argparse
from collections import Counter, defaultdict


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
    ap.add_argument('--top', type=int, default=25)
    args = ap.parse_args()

    src = os.path.join(args.output_dir, 'results_clean.jsonl')
    if not os.path.exists(src):
        print(f"Introuvable : {src}  (lancer d'abord audit_results.py)")
        return 1

    rs = []
    with open(src, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line:
                rs.append(json.loads(line))

    for r in rs:
        h = r.get('cohomology') or [0, 0, 0, 0]
        r['n_anti'] = min(h[1], h[2])          # paires vectorielles reelles
        r['higgs_fiable'] = 'wedge2_heuristique' not in (r.get('warnings') or [])

    print(f"\n{'='*70}")
    print(f"  TRIAGE DE {src}   ({len(rs)} candidats)")
    print(f"{'='*70}")

    print(f"\n  Par type de construction :")
    for k, v in Counter(r['type'] for r in rs).most_common():
        print(f"    {k:<12} {v:>6}")

    ordres = Counter(tuple(r.get('ordres_gamma') or []) for r in rs)
    if any(k for k in ordres):
        print(f"\n  Par ordre du groupe librement agissant :")
        for k, v in ordres.most_common():
            print(f"    |Gamma| in {str(list(k)) if k else 'aucun':<14} {v:>6}")
        print(f"    (l'indice affiche est celui d'AMONT ; sur le quotient il")
        print(f"     est divise par |Gamma|, et les Higgs se decomposent selon")
        print(f"     les representations de Gamma -- pas une simple division)")

    print(f"\n  Par groupe de jauge :")
    for k, v in Counter(r['gauge'] for r in rs).most_common():
        print(f"    {k:<12} {v:>6}")

    print(f"\n  Paires vectorielles n_anti = min(h1,h2) "
          f"-- 0 = spectre reellement propre :")
    dist = Counter(r['n_anti'] for r in rs)
    for k in sorted(dist):
        marque = "   <<< propre" if k == 0 else ""
        print(f"    n_anti = {k:<3} {dist[k]:>6}{marque}")

    print(f"\n  Signatures de cohomologie les plus frequentes :")
    for k, v in Counter(tuple(r['cohomology']) for r in rs).most_common(8):
        print(f"    {str(list(k)):<18} {v:>6}")

    # Familles : meme fibre (memes charges) reapparaissant sur plusieurs CICYs
    fam = defaultdict(set)
    for r in rs:
        sig = (tuple(tuple(b) for b in r['b_charges']),
               tuple(tuple(c) for c in r['c_charges']))
        fam[sig].add(r['cicy'])
    multi = {k: v for k, v in fam.items() if len(v) > 1}
    print(f"\n  Fibres distincts (par charges)      : {len(fam)}")
    print(f"  Dont partages par plusieurs CICYs   : {len(multi)}")
    if multi:
        top = max(multi.items(), key=lambda kv: len(kv[1]))
        print(f"  Le plus repandu present sur        : {len(top[1])} CICYs")
        print(f"    B = {[list(b) for b in top[0][0]]}")
        print(f"    C = {[list(c) for c in top[0][1]]}")
        print(f"    -> un meme jeu de charges reapparaissant sur des dizaines")
        print(f"       de CICYs distinctes merite une verification manuelle :")
        print(f"       soit une vraie universalite, soit un artefact du calcul.")

    propres = [r for r in rs if r['n_anti'] == 0]
    print(f"\n{'='*70}")
    print(f"  CANDIDATS SANS PAIRE VECTORIELLE : {len(propres)}")
    print(f"{'='*70}")

    if not propres:
        print(f"\n  Aucun. Tous les candidats retenus portent au moins une paire")
        print(f"  16+16bar / 10+10bar / 27+27bar en plus des 3 generations.")
        print(f"  Ce n'est pas redhibitoire -- ces paires peuvent devenir")
        print(f"  massives -- mais aucun modele du lot n'a le spectre propre")
        print(f"  vise au depart.")
        pool = sorted(rs, key=lambda r: (r['n_anti'],
                                         0 if r.get('higgs', 0) >= 1 else 1,
                                         r.get('higgs', 0)))
        entete = "Meilleurs disponibles (n_anti croissant)"
    else:
        pool = sorted(propres, key=lambda r: (0 if r.get('higgs', 0) >= 1 else 1,
                                              r.get('higgs', 0)))
        entete = "Candidats propres"

    print(f"\n  {entete} :")
    print(f"    {'#CICY':>6} {'type':<10} {'jauge':>7} {'rk':>2} {'n_anti':>6} "
          f"{'H':>4} {'H fiable':>9} {'cohomologie':>16}")
    for r in pool[:args.top]:
        print(f"    {r['cicy']:>6} {r['type']:<10} {r['gauge']:>7} {r['rank_V']:>2} "
              f"{r['n_anti']:>6} {r.get('higgs',0):>4} "
              f"{'oui' if r['higgs_fiable'] else 'non':>9} "
              f"{str(r['cohomology']):>16}")

    out = os.path.join(args.output_dir, 'results_ranked.jsonl')
    with open(out, 'w', encoding='utf-8') as f:
        for r in sorted(rs, key=lambda r: (r['n_anti'],
                                           0 if r.get('higgs', 0) >= 1 else 1)):
            f.write(json.dumps(r) + '\n')
    print(f"\n  Ecrit : {out}")
    print(f"{'='*70}\n")
    return 0


if __name__ == '__main__':
    sys.exit(main())
