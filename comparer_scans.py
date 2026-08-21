#!/usr/bin/env python3
"""
comparer_scans.py -- Deux balayages d'equivariance, cote a cote, dans les
DEUX SENS.

A QUOI CA REPOND
----------------
`scan_wilson4` melangeait au moins trois versions du code : le sigma du 5.34,
et des lignes « hors domaine » heritees de l'epoque ou `rank_c_max` valait 1.
On savait mesurer les secondes seulement la ou elles se contredisaient
elles-memes -- 41 identites portant a la fois « hors domaine » et « ok ».
Pour les 5 689 identites « hors domaine » SEULES, il aurait fallu rejouer
`domaine_valide` sur chacune.

`scan_wilson5`, recalcule d'un bloc dans un seul etat du code, rend cette
mesure gratuite : il suffit de comparer. Et il la rend MEILLEURE, parce
qu'elle ne porte plus seulement sur le domaine mais sur tous les verdicts.

CE QU'IL DECLARE
----------------
Les deux sens, toujours :

  - identites presentes d'un seul cote (la couverture n'est pas la meme :
    l'ancien balayage etait plafonne a 16 realisations par tache) ;
  - RELIQUES : « hors domaine » d'un cote, evaluees de l'autre -- ce sont
    les candidats ecartes a tort, et le 5.23 dit ce que cela coute : un
    candidat jamais evalue ne se distingue pas d'une absence ;
  - le SENS INVERSE, evaluees puis ecartees. Un rapport qui ne compterait
    que les gains serait un plaidoyer, pas une mesure ;
  - les SURVIT gagnes et les SURVIT perdus, separement. Le solde seul
    masquerait deux erreurs qui se compensent.

Usage :
    python -u comparer_scans.py scan_wilson4 scan_wilson5 \
           --sortie comparaison_w4_w5.json
"""
import os
import sys
import json
import argparse
from collections import defaultdict, Counter

HORS = 'hors domaine (modele S/I non valide)'


def lire(chemin):
    """{(cicy, b, c): {n, survit, etats}} -- l'identite d'un candidat."""
    par = defaultdict(lambda: {'n': 0, 'survit': 0, 'etats': Counter()})
    with open(chemin, encoding='utf-8') as f:
        for ligne in f:
            if not ligne.strip():
                continue
            d = json.loads(ligne)
            cle = (d.get('cicy'), json.dumps(d.get('b_charges')),
                   json.dumps(d.get('c_charges')))
            e = par[cle]
            e['n'] += 1
            e['survit'] += bool(d.get('survit'))
            e['etats'][d.get('etat')] += 1
    return par


def principal():
    ap = argparse.ArgumentParser()
    ap.add_argument('avant')
    ap.add_argument('apres')
    ap.add_argument('--sortie', default=None)
    args = ap.parse_args()

    ch = lambda d: os.path.join(d, 'results_equivariance_f.jsonl')
    for d in (args.avant, args.apres):
        if not os.path.exists(ch(d)):
            sys.exit(f"  ABSENT : {ch(d)}")

    print(f"  lecture de {args.avant} ...", flush=True)
    a = lire(ch(args.avant))
    print(f"    {len(a)} identites", flush=True)
    print(f"  lecture de {args.apres} ...", flush=True)
    b = lire(ch(args.apres))
    print(f"    {len(b)} identites", flush=True)

    ka, kb = set(a), set(b)
    communes = ka & kb
    print(f"\n{'=' * 70}")
    print(f"  COUVERTURE")
    print(f"{'=' * 70}")
    print(f"    communes            : {len(communes)}")
    print(f"    seulement {args.avant:<12}: {len(ka - kb)}")
    print(f"    seulement {args.apres:<12}: {len(kb - ka)}")

    seul_hors = lambda p, k: set(p[k]['etats']) == {HORS}
    reliques = [k for k in communes if seul_hors(a, k) and not seul_hors(b, k)]
    inverse = [k for k in communes if seul_hors(b, k) and not seul_hors(a, k)]

    print(f"\n{'=' * 70}")
    print(f"  DOMAINE -- les deux sens")
    print(f"{'=' * 70}")
    print(f"    « hors domaine » dans {args.avant}, EVALUEES dans {args.apres}")
    print(f"      -> {len(reliques)} identites, sur "
          f"{len({k[0] for k in reliques})} CICYs")
    print(f"         elles apportent {sum(b[k]['survit'] for k in reliques)} "
          f"SURVIT qui n'existaient pas")
    print(f"    sens inverse (evaluees avant, ecartees apres)")
    print(f"      -> {len(inverse)} identites, sur "
          f"{len({k[0] for k in inverse})} CICYs")

    dif = [k for k in communes if a[k]['survit'] != b[k]['survit']]
    gagne = sum(max(0, b[k]['survit'] - a[k]['survit']) for k in dif)
    perdu = sum(max(0, a[k]['survit'] - b[k]['survit']) for k in dif)
    ta = sum(v['survit'] for v in a.values())
    tb = sum(v['survit'] for v in b.values())
    print(f"\n{'=' * 70}")
    print(f"  SURVIT -- les deux sens")
    print(f"{'=' * 70}")
    print(f"    identites dont le compte change : {len(dif)}")
    print(f"      gagnes {gagne}   perdus {perdu}   "
          f"(solde {gagne - perdu:+d})")
    print(f"    total : {ta} -> {tb}")
    if gagne and perdu and ta == tb:
        print(f"    ATTENTION : le total est identique alors que {gagne} "
              f"SURVIT apparaissent et {perdu} disparaissent.")
        print(f"    Un solde nul n'est pas une absence de changement.")
    print(f"    CICYs les plus concernees : "
          f"{Counter(k[0] for k in dif).most_common(8)}")

    if args.sortie:
        with open(args.sortie, 'w', encoding='utf-8') as f:
            json.dump({'avant': args.avant, 'apres': args.apres,
                       'couverture': {'communes': len(communes),
                                      'seulement_avant': len(ka - kb),
                                      'seulement_apres': len(kb - ka)},
                       'reliques': [list(k) for k in reliques],
                       'inverse': [list(k) for k in inverse],
                       'survit': {'change': [list(k) for k in dif],
                                  'gagnes': gagne, 'perdus': perdu,
                                  'total_avant': ta, 'total_apres': tb}},
                      f)
        print(f"\n  rapport : {args.sortie}")


if __name__ == '__main__':
    principal()
