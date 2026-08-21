#!/usr/bin/env python3
"""
verdict_z4.py -- Stabilite des candidats Z4 du 2.3, maintenant que le chemin
wedge et le certificat de surjectivite traitent rank_C = 2.

Ces candidats ont trois generations DEMONTREES (H^1(V) = representation
reguliere, §5.34) et h^0(V) equivariant nul, mais aucun verdict de stabilite :
ce sont des monades a rank_C = 2, et les deux fonctions qui decident --
`hoppe_sur_espace` et `f_sans_point_base` -- y etaient restreintes a
rank_C = 1.

DECLARE LES DEUX COTES. Pour chaque (candidat, realisation, lambda) :
h^0(V) generique et equivariant, le detail p par p du critere de Hoppe, et le
certificat de surjectivite avec son degre -- ou la raison de son echec. Un
`None` reste un `None` : non calculable n'est pas negatif, et surtout pas
positif.

Usage :
    python -u verdict_z4.py cicyquotients.m cicylist.txt \
           [--input scan_wilson5/results_equivariant.jsonl]
"""
import os
import sys
import json
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from cy_landscape.core.braun_symmetry import (parse_symmetries, ordres_rt,
                                              matrice_mod_p)
from cy_landscape.core.gamma_action import choisir_premier, racine_primitive
from cy_landscape.core.covariant_ring import (resoudre_covariants,
                                              tirer_covariants, CovariantRing,
                                              verifier_covariance,
                                              permutation_facteurs_numerique)
from cy_landscape.core.equivariant_monad import (espace_f_equivariant,
                                                 hoppe_sur_espace,
                                                 f_sans_point_base,
                                                 h0_V_sur_espace,
                                                 decomposition_h1_V)
from cy_landscape.core.sections import domaine_valide
from cy_landscape.data.parse_oxford import load_oxford_file
from wilson_match import parse_braun, parse_cicylist, apparier
from equivariance_f import h0_V_generique

GRAINE = 0


def charges_non_certifiees(amb, cfg, b, c):
    """
    Les charges que `domaine_valide` exige et que le critere de degre ne
    certifie pas. On les NOMME au lieu de rendre un booleen : « hors
    domaine » sur 1 charge de 27 et sur 27 de 27 sont deux situations
    differentes, et le booleen les confond.
    """
    from cy_landscape.core.exact_cohomology import koszul_cohomology_ex
    m = len(amb)
    ch = [(f'b{i}', list(x)) for i, x in enumerate(b)]
    ch += [(f'c{j}', list(y)) for j, y in enumerate(c)]
    ch += [(f'b{i}+b{j}', [b[i][k] + b[j][k] for k in range(m)])
           for i in range(len(b)) for j in range(i + 1, len(b))]
    for j, y in enumerate(c):
        ch += [(f'c{j}+b{i}', [y[k] + x[k] for k in range(m)])
               for i, x in enumerate(b)]
    ko = []
    for nom, x in ch:
        r = koszul_cohomology_ex(amb, cfg, x)
        if not r['certified_by_degree'][0]:
            ko.append((nom, tuple(x)))
    return ko, len(ch)


def traiter(cicy, b, c, e, syms, nom_groupe='Z4', forcer=False):
    amb, cfg = e['ambient'], np.asarray(e['config'])
    rang_V = len(b) - len(c)
    print(f"\n{'=' * 78}")
    print(f"  #{cicy}   ambient {amb}   rank_B = {len(b)}, rank_C = {len(c)}, "
          f"rank_V = {rang_V}")
    print(f"    b = {b}")
    print(f"    c = {c}")
    ko, n_ch = charges_non_certifiees(amb, cfg, b, c)
    if ko:
        print(f"    HORS DOMAINE sur {len(ko)} charge(s) de {n_ch} : "
              + ", ".join(f"{nom} = {x}" for nom, x in ko))
        if not forcer:
            print("    -> pas de verdict. --forcer pour calculer quand meme, "
                  "sous reserve declaree.")
            return []
        print("    -> --forcer : on calcule, ET ON L'ESTAMPILLE. Le modele "
              "S/I n'est pas certifie sur")
        print("       cette charge, donc AUCUN de ces chiffres n'est un "
              "resultat tant que le bloc")
        print("       `corr` du §5.32 ne l'a pas etablie. Ils disent seulement "
              "ce que le calcul rendrait.")
    elif not domaine_valide(amb, cfg, b, c, rank_c_max=None):
        print("    HORS DOMAINE pour une autre raison que la certification")
        if not forcer:
            return []
    dets = [sum(x[k] for x in b) - sum(y[k] for y in c)
            for k in range(len(amb))]
    print(f"    det V = O({dets})"
          + ("  -- c1(V) = 0, le critere de Hoppe s'applique"
             if not any(dets) else "  /!\\ c1(V) != 0"))

    lignes = []
    for n_sym, sym in enumerate(syms):
        if sym['nom'] != nom_groupe:
            continue
        ordres = ordres_rt(sym['coord']) | ordres_rt(sym['poly']) | {2, 4}
        p, _ = choisir_premier(sorted(ordres), minimum=30011)
        rac = {n: racine_primitive(p, n) for n in ordres}
        Mc = [matrice_mod_p(x, p, rac) for x in sym['coord']]
        Np = [matrice_mod_p(x, p, rac) for x in sym['poly']]
        sig = permutation_facteurs_numerique(Mc[0], amb, p)
        res = resoudre_covariants(amb, cfg, Mc, Np, p)
        if res is None:
            print(f"  realisation {n_sym} : sigma non extractible")
            continue
        v = res['par_convention'].get('N')
        if v is None or not v['non_degenere']:
            print(f"  realisation {n_sym} : aucun ideal covariant non degenere")
            continue
        co = tirer_covariants(v['base'], res['offsets'], res['dims'], p,
                              np.random.RandomState(GRAINE))
        ok, ec = verifier_covariance(amb, cfg, Mc, Np, 'N', co, p)
        if not ok:
            print(f"  realisation {n_sym} : covariance non reverifiee ({ec})")
            continue
        anneau = CovariantRing(amb, cfg, co, p)
        out = espace_f_equivariant(anneau, amb, b, c, Mc, p)
        if out['etat'] != 'ok' or not out['solutions']:
            print(f"  realisation {n_sym} : {out['etat']}, "
                  f"{len(out.get('solutions') or [])} solution(s)")
            continue

        h_gen = h0_V_generique(anneau, b, c, p,
                               np.random.RandomState(GRAINE + 5))
        print(f"\n  realisation {n_sym}   sigma = {sig}   "
              f"dim totale de l'espace des f = {out['dim_totale']}   "
              f"h0(V) generique = {h_gen}")
        print(f"  {'lambda':>14} {'dim eq':>7} {'h0 eq':>6} "
              f"{'Hoppe':>7}  {'valeurs h0(wedge^p V)':<28} surjectivite")

        for s in out['solutions']:
            lam = tuple(int(x) if x < p // 2 else int(x) - p
                        for x in s['lambda'])
            h_eq, _ = h0_V_sur_espace(anneau, amb, b, c, s['base'],
                                      out['cases'], out['offsets'],
                                      out['dims'], out['degres'], p,
                                      np.random.RandomState(GRAINE + 5))
            hop = hoppe_sur_espace(anneau, b, c, s['base'], out['offsets'],
                                   out['dims'], out['degres'], p,
                                   np.random.RandomState(GRAINE + 5))
            surj = f_sans_point_base(anneau, b, c, s['base'], out['offsets'],
                                     out['dims'], out['degres'], p,
                                     np.random.RandomState(GRAINE + 7),
                                     n_essais=2, n_degres=8)
            etat_h = {True: 'stable', False: 'NON', None: '?'}[hop['stable']]
            sj = (f"certifie en {surj['degre']}" if surj['certifie']
                  else f"NON ({surj.get('motif', 'aucun degre concluant')})")
            print(f"  {str(lam):>14} {s['dim']:>7} {h_eq:>6} {etat_h:>7}  "
                  f"{str(hop['valeurs']):<28} {sj}")
            if hop['stable'] is not True and hop['motif']:
                print(f"                 motif Hoppe : {hop['motif']}")
            lignes.append({'cicy': cicy, 'realisation': n_sym,
                           'sigma': None if sig is None else list(sig),
                           'lambda': list(lam), 'h0_generique': h_gen,
                           'h0_equivariant': h_eq,
                           'hoppe': hop['stable'],
                           'hoppe_valeurs': {str(k): v
                                             for k, v in hop['valeurs'].items()},
                           'hoppe_motif': hop['motif'],
                           'surjectif': bool(surj['certifie']),
                           'surjectif_degre': surj['degre'],
                           'stable_demontre': bool(hop['stable'] is True
                                                   and surj['certifie']
                                                   and h_eq == 0
                                                   and h_gen == 0)})
    return lignes


def principal():
    ap = argparse.ArgumentParser()
    ap.add_argument('braun_m')
    ap.add_argument('cicylist')
    ap.add_argument('--input',
                    default='scan_wilson5/results_equivariant.jsonl')
    ap.add_argument('--sortie', default='verdict_z4.json')
    ap.add_argument('--forcer', action='store_true',
                    help="calculer meme si le modele S/I n'est pas certifie "
                         "sur toutes les charges. Chaque ligne porte alors "
                         "`charges_non_certifiees` : ce sont des chiffres "
                         "sous reserve, pas des resultats.")
    args = ap.parse_args()

    entries = {e['num']: e for e in load_oxford_file(args.cicylist)}
    braun = parse_braun(args.braun_m)
    cl = parse_cicylist(args.cicylist)
    corr, _, _ = apparier(braun, cl)
    inv = {v: k for k, v in corr.items()}
    SYM = parse_symmetries(args.braun_m)

    # Les candidats du 2.3 : on les reprend du fichier d'entree, pour ne pas
    # recopier des charges a la main.
    vus = {}
    with open(args.input, encoding='utf-8') as f:
        for L in f:
            if not L.strip():
                continue
            d = json.loads(L)
            if d.get('cicy') not in (7745, 6947):
                continue
            if 'Z4' not in (d.get('groupes_utiles') or []):
                continue
            if len(d.get('c_charges') or []) != 2:
                continue
            cle = (d['cicy'], json.dumps(d['b_charges']),
                   json.dumps(d['c_charges']))
            vus.setdefault(cle, d)
    print(f"  {len(vus)} candidats Z4 a rank_C = 2 trouves dans {args.input}")

    toutes = []
    for (cicy, bs, cs), d in sorted(vus.items()):
        e = entries.get(cicy)
        nb = inv.get(cicy)
        if e is None or nb not in SYM:
            print(f"  #{cicy} : absente ou non appariee")
            continue
        lg = traiter(cicy, json.loads(bs), json.loads(cs), e,
                     SYM[nb]['symetries'], forcer=args.forcer)
        ko, n_ch = charges_non_certifiees(e['ambient'], np.asarray(e['config']),
                                          json.loads(bs), json.loads(cs))
        for x in lg:
            x['charges_non_certifiees'] = [[nom, list(v)] for nom, v in ko]
            x['charges_totales'] = n_ch
            if ko:
                # REGLE DES FILTRES : un verdict pris sur un modele non
                # certifie n'est pas un verdict. On le dit dans la donnee,
                # pas seulement a l'ecran.
                x['stable_demontre'] = False
                x['sous_reserve'] = True
        toutes.extend(lg)

    print(f"\n{'=' * 78}")
    print("  BILAN, LES DEUX COTES")
    print(f"{'=' * 78}")
    dem = [x for x in toutes if x['stable_demontre']]
    tue = [x for x in toutes if x['hoppe'] is False]
    ind = [x for x in toutes if x['hoppe'] is None
           or (x['hoppe'] is True and not x['surjectif'])]
    print(f"    lignes (candidat, realisation, lambda) evaluees : {len(toutes)}")
    print(f"      stabilite DEMONTREE                : {len(dem)}")
    print(f"      ELIMINES par Hoppe                 : {len(tue)}")
    print(f"      indetermines (un test non concluant): {len(ind)}")
    for x in ind:
        motif = (x['hoppe_motif'] if x['hoppe'] is not True
                 else 'surjectivite non certifiee')
        print(f"        #{x['cicy']} real {x['realisation']} "
              f"lambda {x['lambda']} : {motif}")
    with open(args.sortie, 'w', encoding='utf-8') as f:
        json.dump(toutes, f, indent=1)
    print(f"\n  rapport : {args.sortie}")


if __name__ == '__main__':
    principal()
