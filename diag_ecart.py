#!/usr/bin/env python3
"""
diag_ecart.py -- Pourquoi une identite ne reproduit-elle pas la reference ?

A QUOI CA REPOND
----------------
`ancres_port.py` compare a la maille de l'IDENTITE (cicy, b, c, groupe) et
suppose, pour deduire son facteur r, que les r realisations de Gamma rendent
le meme verdict. Quand elles ne le rendent PAS, il signale un ecart -- et il a
raison de le faire, parce que les deux lectures possibles sont opposees :

    (a) les realisations DIVERGENT, chacune sur son propre quotient X/Gamma.
        C'est un fait, mesure au §5.35 : sigma est une propriete de la
        REALISATION, pas du nom du groupe. Rien n'est faux ;
    (b) le balayage PERD des survivants -- verdict usurpe dans l'autre sens,
        le danger nomme au §5.34.

Les deux donnent le meme comptage agrege. Seule la ventilation REALISATION
PAR REALISATION les distingue, et c'est ce que ce script imprime.

CE QU'IL REGARDE
----------------
Chaque ligne du JSONL porte son `_lot` = ('T', tache, tranche). Une tranche
vaut `--taille-lot` realisations consecutives (16 par defaut), donc la
tranche est un PROXY de la realisation, pas la realisation elle-meme : le
balayage n'ecrit pas l'indice de realisation. On ventile donc par tranche, et
par (lambda, dim_equivariant) -- deux realisations de sigma differents
donnent des dimensions d'espace equivariant differentes, et c'est la
signature la plus lisible qu'on ait sans rejouer le calcul.

Lecture seule. N'ecrit rien.

Usage :
    python -u diag_ecart.py scan_wilson6 --cicy 480 --groupe "Z2 x Z2"
    python -u diag_ecart.py scan_wilson6 --cicy 480          # tous les groupes
"""
import os
import sys
import json
import argparse
from collections import Counter, defaultdict


def verdict(d):
    """La meme classification que `ancres_port.lire_scan`."""
    if d.get('survit'):
        return 'SURVIT'
    if d.get('indetermine'):
        return 'indetermine'
    if d.get('fibre') is False:
        return 'elimine (lieu de base)'
    if d.get('hoppe_complet') is False:
        return 'elimine (Hoppe)'
    if d.get('h0_equivariant'):
        return 'elimine (h0 equivariant)'
    if d.get('h0_generique'):
        return 'elimine (h0 generique)'
    return 'elimine (autre)'


def principal():
    ap = argparse.ArgumentParser()
    ap.add_argument('cible', help='dossier de scan, ou chemin du JSONL')
    ap.add_argument('--cicy', type=int, required=True)
    ap.add_argument('--groupe', default=None)
    ap.add_argument('--strate', default=None,
                    help='ne garder qu une strate, ex. "1,4"')
    ap.add_argument('--max-identites', type=int, default=6)
    args = ap.parse_args()
    strate_voulue = None
    if args.strate:
        strate_voulue = tuple(int(x) for x in args.strate.split(','))

    chemin = args.cible
    if os.path.isdir(chemin):
        chemin = os.path.join(chemin, 'results_equivariance_f.jsonl')
    if not os.path.exists(chemin):
        sys.exit(f"  ABSENT : {chemin}")

    # {identite: [lignes]}
    par = defaultdict(list)
    n = 0
    with open(chemin, encoding='utf-8') as f:
        for ligne in f:
            if not ligne.strip():
                continue
            if f'"cicy": {args.cicy}' not in ligne:
                continue
            d = json.loads(ligne)
            if d.get('cicy') != args.cicy or d.get('etat') != 'ok':
                continue
            if args.groupe and d.get('groupe') != args.groupe:
                continue
            if strate_voulue is not None:
                rc_ = len(d.get('c_charges') or [])
                if (rc_, d.get('rang_V')) != strate_voulue:
                    continue
            n += 1
            cle = (json.dumps(d.get('b_charges')),
                   json.dumps(d.get('c_charges')), d.get('groupe'))
            par[cle].append(d)

    print(f"\n  {chemin}")
    print(f"  #{args.cicy}"
          + (f" / {args.groupe}" if args.groupe else "")
          + f" : {n} lignes `ok`, {len(par)} identites\n")

    # RESUME : une ligne par identite, pour que le detail qui suit ne soit
    # pas la seule vue -- six identites affichees sur cent dix-huit, c'est un
    # echantillon, et un echantillon qu'on ne declare pas se lit comme un
    # tout (§8).
    print(f"  RESUME -- {len(par)} identites, toutes, par nombre de lignes :")
    entete_t = "lots 'T'"
    print(f"    {'lignes':>7} {entete_t:>9}  {'strate':<12} verdicts")
    print(f"    (les lignes des lots 'C' -- controle d'orbite, doublons "
          f"deliberes du §5.25 -- sont comptees a part)")
    for cle, lignes in sorted(par.items(), key=lambda x: -len(x[1])):
        rc = len(json.loads(cle[1]))
        rv = len(json.loads(cle[0])) - rc
        lots_t = {tuple(d['_lot']) for d in lignes
                  if d.get('_lot') and d['_lot'][0] == 'T'}
        lots_c = {tuple(d['_lot']) for d in lignes
                  if d.get('_lot') and d['_lot'][0] == 'C'}
        n_c = sum(1 for d in lignes
                  if d.get('_lot') and d['_lot'][0] == 'C')
        vs = Counter(verdict(d) for d in lignes
                     if not (d.get('_lot') and d['_lot'][0] == 'C'))
        detail = ", ".join(f"{v} {k}" for k, v in vs.most_common())
        marque = f"  [+{n_c} lignes de CONTROLE, {len(lots_c)} lots 'C']" if n_c else ""
        print(f"    {len(lignes) - n_c:>7} {len(lots_t):>9}  ({rc},{rv})".ljust(36)
              + f"  {detail}{marque}")
    print()

    for i, (cle, lignes) in enumerate(
            sorted(par.items(), key=lambda x: -len(x[1]))[:args.max_identites], 1):
        b, c, g = cle
        rc = len(json.loads(c))
        rv = len(json.loads(b)) - rc
        print("=" * 78)
        print(f"  [{i}] {g}   strate (rank_C={rc}, rang_V={rv})   "
              f"{len(lignes)} lignes")
        print(f"      b = {b}")
        print(f"      c = {c}")

        # --- ventilation globale ------------------------------------
        glob = Counter(verdict(d) for d in lignes)
        print(f"\n      verdicts, tous confondus :")
        for k, v in glob.most_common():
            print(f"        {v:>6}  {k}")

        # --- ventilation par TRANCHE (proxy de realisation) ----------
        # C'est ici que se decide la lecture (a) contre la lecture (b).
        par_tr = defaultdict(Counter)
        for d in lignes:
            lot = d.get('_lot') or []
            tr = lot[2] if len(lot) > 2 else None
            par_tr[tr][verdict(d)] += 1
        profils = Counter()
        for tr, c_ in sorted(par_tr.items(), key=lambda x: (x[0] is None, x[0])):
            profils[tuple(sorted(c_.items()))] += 1
        print(f"\n      {len(par_tr)} tranches, {len(profils)} PROFILS distincts :")
        for prof, combien in profils.most_common():
            detail = ", ".join(f"{v} {k}" for k, v in prof)
            print(f"        {combien:>4} tranches  ->  {detail}")
        if len(profils) == 1:
            print(f"        -> toutes les tranches rendent le MEME profil : "
                  f"les realisations s'accordent.")
            print(f"           Un ecart avec la reference ne vient donc PAS "
                  f"de sigma.")
        else:
            print(f"        -> les tranches DIVERGENT. C'est la signature de "
                  f"sigma (§5.35) :")
            print(f"           deux realisations, deux quotients X/Gamma, "
                  f"deux verdicts legitimes.")

        # --- signature de sigma : la dimension de l'espace equivariant
        dims = Counter((d.get('lambda') and tuple(d['lambda']),
                        d.get('dim_equivariant'), d.get('dim_totale'))
                       for d in lignes)
        print(f"\n      (lambda, dim equivariant / dim totale) :")
        for (lam, de, dt), v in dims.most_common(12):
            print(f"        {v:>6}  lambda={lam}  {de}/{dt}")
        if len({d.get('dim_equivariant') for d in lignes}) > 1:
            print(f"        -> plusieurs dimensions d'espace equivariant pour "
                  f"la MEME identite :")
            print(f"           les realisations ne portent pas le meme sigma.")
        print()

    if not par:
        print("  aucune ligne trouvee -- le lot ne contient pas encore cette "
              "identite ?\n")
    return 0


if __name__ == '__main__':
    sys.exit(principal())
