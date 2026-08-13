#!/usr/bin/env python3
"""
equivariance_f.py -- Test d'equivariance de f : B -> C, et stabilite restreinte.

--------------------------------------------------------------------------
Ce que ce script ajoute a `equivariance.py`
--------------------------------------------------------------------------
`equivariance.py` teste une condition NECESSAIRE portant sur les CHARGES :
sigma doit permuter l'ensemble des b_i et des c_j. Elle ne dit rien des
polynomes, et elle est vide sur les CICYs ou Gamma n'agit que par des phases.

Ce script descend au niveau des polynomes, en trois etapes :

  1. POLYNOMES COVARIANTS. Les coefficients des K polynomes definissants
     sont resolus pour que Gamma preserve l'ideal (`covariant_ring`), au lieu
     d'etre tires au hasard. Sans cela l'action ne descend pas au quotient et
     tout ce qui suit serait un calcul sur une autre variete.

  2. EXISTENCE D'UN f EQUIVARIANT. On resout S_g(f_{j,i}) = lambda_g
     f_{rho(j),pi(i)} pour tous les generateurs (`equivariant_monad`), en
     ENUMERANT les lambda_g admissibles a partir de l'ordre projectif de
     l'operateur -- ce qui traite le relevement projectif au lieu de le
     supposer trivial.

  3. STABILITE RESTREINTE. C'est l'etape qui mord. Sur un Gamma qui agit par
     phases, le sous-espace equivariant represente une fraction 1/|Gamma| de
     l'espace des f et n'est donc JAMAIS vide : conclure de l'etape 2 que le
     fibre descend serait une erreur. On recalcule donc h^0(V) avec un f tire
     dans le sous-espace equivariant. Si h^0(V) devient non nul, V n'est plus
     stable et le candidat tombe.

--------------------------------------------------------------------------
Domaine
--------------------------------------------------------------------------
Les etapes 2 et 3 utilisent le modele R_a = S_a / I_a, qui ne represente
H^0(Y, O(a)) que sur le domaine verifie par `sections.domaine_valide`.
Hors de ce domaine le script n'affiche pas de verdict : il marque
`hors domaine`. Un candidat hors domaine n'est ni retenu ni elimine.

Usage:
    python equivariance_f.py cicyquotients.m cicylist.txt scan_wilson2
    python equivariance_f.py cicyquotients.m cicylist.txt scan_wilson2 --cicy 6947
"""
import os
import sys
import json
import argparse

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cy_landscape.core.braun_symmetry import (parse_symmetries, ordres_rt,
                                              matrice_mod_p)
from cy_landscape.core.gamma_action import choisir_premier, racine_primitive
from cy_landscape.core.covariant_ring import (resoudre_covariants,
                                              tirer_covariants,
                                              verifier_covariance,
                                              CovariantRing)
from cy_landscape.core.equivariant_monad import (espace_f_equivariant,
                                                 h0_V_sur_espace,
                                                 h0_wedge2_V_sur_espace,
                                                 hoppe_sur_espace,
                                                 f_sans_point_base)
from cy_landscape.core.sections import (basis_multi, rref_mod, _mult_matrix,
                                        domaine_valide)
from cy_landscape.data.parse_oxford import load_oxford_file
from wilson_match import parse_braun, parse_cicylist, apparier


def h0_V_generique(anneau, b_charges, c_charges, p, rng, n_essais=5):
    """h^0(V) pour un f generique, dans le MEME anneau covariant.

    C'est la valeur de reference : la comparer a `h0_V_sur_espace` isole
    l'effet de la contrainte d'equivariance, sans melanger l'effet du
    passage aux polynomes covariants.
    """
    m = len(b_charges[0])
    dsrc = sum(anneau.dimY(list(b)) for b in b_charges)
    meilleur = None
    for _ in range(n_essais):
        lignes = []
        for cj in c_charges:
            ddst = anneau.dimY(list(cj))
            ligne = []
            for b in b_charges:
                deg = [cj[k] - b[k] for k in range(m)]
                if any(x < 0 for x in deg):
                    ligne.append(np.zeros((ddst, anneau.dimY(list(b))),
                                          dtype=np.int64))
                    continue
                fb = basis_multi(anneau.amb, deg)
                ligne.append(_mult_matrix(anneau, list(b), deg,
                                          (fb, rng.randint(1, p, size=len(fb))),
                                          list(cj)))
            lignes.append(np.hstack(ligne))
        rang, _ = rref_mod(np.vstack(lignes).T.copy(), p)
        h0 = dsrc - rang
        meilleur = h0 if meilleur is None else min(meilleur, h0)
    return meilleur


def analyser(cicy_num, amb, cfg, b, c, symetries, groupes=None, graine=0):
    """Renvoie une liste de lignes de resultat, une par (symetrie, lambda)."""
    lignes = []
    # rank_c_max=None : on autorise rank_C >= 2. Les fonctions generalisees
    # (h^0(V), decomposition de H^1(V)) s'y appliquent ; celles qui supposent
    # rank_C = 1 -- wedge^p V, surjectivite -- se declarent elles-memes non
    # calculables et le verdict tombe en `indetermine`, jamais en succes.
    if not domaine_valide(amb, cfg, b, c, rank_c_max=None):
        return [{'groupe': '-', 'etat': 'hors domaine (modele S/I non valide)'}]

    for sym in symetries:
        if groupes and sym['nom'] not in groupes:
            continue
        ordres = ordres_rt(sym['coord']) | ordres_rt(sym['poly']) | {2}
        try:
            p, _ = choisir_premier(sorted(ordres), minimum=30011)
            rac = {n: racine_primitive(p, n) for n in ordres}
            Mc = [matrice_mod_p(x, p, rac) for x in sym['coord']]
            Np = [matrice_mod_p(x, p, rac) for x in sym['poly']]
        except Exception as exc:
            lignes.append({'groupe': sym['nom'],
                           'etat': f'generateurs illisibles ({exc})'})
            continue

        res = resoudre_covariants(amb, cfg, Mc, Np, p)
        if res is None:
            lignes.append({'groupe': sym['nom'], 'etat': 'sigma non extractible'})
            continue
        bons = [k for k, v in res['par_convention'].items() if v['non_degenere']]
        if 'N' not in bons:
            lignes.append({'groupe': sym['nom'],
                           'etat': 'aucun ideal covariant non degenere'})
            continue
        v = res['par_convention']['N']
        co = tirer_covariants(v['base'], res['offsets'], res['dims'], p,
                              np.random.RandomState(graine))
        ok, ec = verifier_covariance(amb, cfg, Mc, Np, 'N', co, p)
        if not ok:
            lignes.append({'groupe': sym['nom'],
                           'etat': f'covariance non revérifiée (ecart {ec})'})
            continue

        anneau = CovariantRing(amb, cfg, co, p)
        out = espace_f_equivariant(anneau, amb, b, c, Mc, p)
        if out['etat'] != 'ok':
            lignes.append({'groupe': sym['nom'], 'etat': out['etat']})
            continue
        if not out['solutions']:
            lignes.append({'groupe': sym['nom'], 'etat': 'aucun f equivariant',
                           'dim_totale': out['dim_totale'], 'elimine': True})
            continue

        rang_V = len(b) - len(c)
        # Hoppe, c1(V) = 0 : rk 3 -> h0(V), h3(V) ; rk >= 4 -> + h0(w2V).
        # w2V n'est calculable par ce chemin que si rank_C = 1, comme dans
        # hoppe_fast. Au-dela on ne conclut pas.
        besoin_w2 = rang_V >= 4
        w2_possible = besoin_w2 and len(c) == 1

        h_gen = h0_V_generique(anneau, b, c, p, np.random.RandomState(graine + 5))
        base_tot = np.eye(out['dim_totale'], dtype=np.int64)
        w2_gen = None
        if w2_possible:
            w2_gen, _ = h0_wedge2_V_sur_espace(
                anneau, b, c, base_tot, out['offsets'], out['dims'],
                out['degres'], p, np.random.RandomState(graine + 5))

        for s in out['solutions']:
            h_eq, _ = h0_V_sur_espace(anneau, amb, b, c, s['base'],
                                      out['cases'], out['offsets'], out['dims'],
                                      out['degres'], p,
                                      np.random.RandomState(graine + 5))
            w2_eq = None
            if w2_possible:
                w2_eq, _ = h0_wedge2_V_sur_espace(
                    anneau, b, c, s['base'], out['offsets'], out['dims'],
                    out['degres'], p, np.random.RandomState(graine + 5))
            lam = tuple(int(x) if x < p // 2 else int(x) - p for x in s['lambda'])

            # « survit » exige TOUS les tests disponibles a ce rang. Un test
            # non calculable (w2 hors taille, ou rank_C >= 2) rend le verdict
            # indetermine -- jamais favorable.
            if h_gen != 0 or h_eq != 0:
                survit, indetermine = False, False
            elif not besoin_w2:
                survit, indetermine = True, False
            elif not w2_possible or w2_eq is None or w2_gen is None:
                survit, indetermine = False, True
            else:
                survit = (w2_gen == 0 and w2_eq == 0)
                indetermine = False

            # Surjectivite : V = ker(f) n'est un FIBRE que si f est surjective
            # en tout point. Rien ne le garantit sur le sous-espace
            # equivariant. Le critere est SUFFISANT (voir f_sans_point_base) :
            # un echec ne prouve pas la non-surjectivite, il empeche seulement
            # de conclure -- d'ou `indetermine` et non une elimination.
            # Teste uniquement quand tout le reste passe, pour le cout.
            # Critere de Hoppe COMPLET sur le sous-espace : h0(wedge^p V) = 0
            # pour p = 1..rk-1. Inclut h^3(V) = h0(wedge^{rk-1} V), qui n'etait
            # jusqu'ici traite nulle part sous contrainte. Un p non calculable
            # rend le verdict indetermine, jamais favorable.
            hoppe = None
            if survit and len(c) == 1:
                hoppe = hoppe_sur_espace(
                    anneau, b, c, s['base'], out['offsets'], out['dims'],
                    out['degres'], p, np.random.RandomState(graine + 5))
                if hoppe['stable'] is not True:
                    survit = False
                    indetermine = (hoppe['stable'] is None)

            surj = None
            if survit:
                surj = f_sans_point_base(
                    anneau, b, c, s['base'], out['offsets'], out['dims'],
                    out['degres'], p, np.random.RandomState(graine + 7),
                    n_essais=2, n_degres=8)
                if not surj['certifie']:
                    survit, indetermine = False, True

            lignes.append({
                'groupe': sym['nom'], 'etat': 'ok', 'lambda': lam,
                'dim_equivariant': s['dim'], 'dim_totale': out['dim_totale'],
                'rang_V': rang_V,
                'h0_generique': h_gen, 'h0_equivariant': h_eq,
                'h0w2_generique': w2_gen, 'h0w2_equivariant': w2_eq,
                'hoppe_complet': None if hoppe is None else hoppe['stable'],
                'hoppe_valeurs': None if hoppe is None else
                                 {str(k): v for k, v in hoppe['valeurs'].items()},
                'surjectif_certifie': None if surj is None else bool(surj['certifie']),
                'surjectif_degre': None if surj is None else surj['degre'],
                'survit': bool(survit), 'indetermine': bool(indetermine),
            })
    return lignes


def _sortie_tolerante():
    """
    Empeche un plantage d'encodage sur une console Windows.

    Les etiquettes de groupe de jauge contiennent des indices Unicode
    ("E₆"), que la console cp1252 ne sait pas encoder : `print` levait
    UnicodeEncodeError des la premiere ligne de resultat, apres plusieurs
    minutes de calcul et sans rien ecrire. On passe donc stdout en
    errors='replace' -- le caractere devient '?', le calcul continue. Le
    fichier JSONL de sortie, lui, est ecrit en UTF-8 explicite et garde
    l'etiquette exacte.
    """
    for flux in (sys.stdout, sys.stderr):
        try:
            flux.reconfigure(errors='replace')
        except (AttributeError, ValueError):
            pass


def main():
    _sortie_tolerante()
    ap = argparse.ArgumentParser()
    ap.add_argument('braun_m')
    ap.add_argument('cicylist')
    ap.add_argument('output_dir')
    ap.add_argument('--input', default='results_equivariant.jsonl')
    ap.add_argument('--cicy', type=int, default=None)
    ap.add_argument('--tous-groupes', action='store_true',
                    help="Ne pas se limiter aux groupes d'ordre compatible "
                         "avec l'indice (champ groupes_utiles).")
    args = ap.parse_args()

    entries = {e['num']: e for e in load_oxford_file(args.cicylist)}
    braun = parse_braun(args.braun_m)
    cl = parse_cicylist(args.cicylist)
    corr, _, _ = apparier(braun, cl)
    inv = {v: k for k, v in corr.items()}
    SYM = parse_symmetries(args.braun_m)

    src = os.path.join(args.output_dir, args.input)
    rs = [json.loads(l) for l in open(src, encoding='utf-8') if l.strip()]
    if args.cicy:
        rs = [r for r in rs if r['cicy'] == args.cicy]

    print(f"\n{'=' * 96}")
    print("  EQUIVARIANCE DE f  --  polynomes covariants, puis stabilite restreinte")
    print(f"{'=' * 96}")
    print(f"  {'CICY':>5} {'jauge':>7} {'rk':>2} {'groupe':<11} {'lambda':>14} "
          f"{'dim eq':>6} {'N':>5} {'h0 gen':>6} {'h0 eq':>5} "
          f"{'w2 gen':>6} {'w2 eq':>5}  verdict")

    sortie = []
    survivants = 0
    indetermines = 0
    ecartes = 0
    for r in rs:
        e = entries.get(r['cicy'])
        num_b = inv.get(r['cicy'])
        if e is None or num_b not in SYM:
            continue
        amb, cfg = e['ambient'], np.asarray(e['config'])
        b = [list(x) for x in r['b_charges']]
        c = [list(x) for x in r['c_charges']]
        groupes = None if args.tous_groupes else set(r.get('groupes_utiles') or [])
        if groupes is not None and not groupes:
            groupes = set(r.get('equivariant_possible') or [])
        lignes = analyser(r['cicy'], amb, cfg, b, c, SYM[num_b]['symetries'],
                          groupes=groupes)
        # Identite du candidat, recopiee sur CHAQUE ligne de sortie.
        ident = {k: r.get(k) for k in
                 ('cicy', 'gauge', 'rank_V', 'cohomology',
                  'b_charges', 'c_charges', 'groupes_utiles',
                  'equivariant_possible', 'ordres_gamma')}

        for L in lignes:
            if L['etat'] != 'ok':
                # ECRITES ELLES AUSSI. La version precedente ne persistait que
                # les lignes 'ok' : les « hors domaine », « charges non
                # permutees » et « espace trop grand » n'existaient que dans la
                # sortie console. Sur le balayage `scan_gros_gamma`, cela a
                # rendu le JSONL trompeur -- il montrait 0 couple sur un groupe
                # d'ordre compatible, alors que 26 candidats en portaient un et
                # avaient simplement ete ecartes en amont, pour une raison que
                # le fichier ne contenait pas. Un fichier de resultats doit
                # dire pourquoi un cas n'a pas ete traite, sinon son silence se
                # lit a tort comme une absence de candidats.
                print(f"  {r['cicy']:>5} {r.get('gauge', ''):>7} "
                      f"{L['groupe']:<11} {L['etat']}")
                ecartes += 1
                sortie.append({**ident, **L, 'survit': False,
                               'indetermine': True})
                continue
            if L['survit']:
                verdict = (f"SURVIT (Hoppe complet + surjectif en "
                           f"{L['surjectif_degre']})")
            elif L.get('hoppe_complet') is False:
                verdict = f"tue par Hoppe complet : {L.get('hoppe_valeurs')}"
            elif L.get('hoppe_complet') is None and L.get('indetermine') \
                    and L.get('surjectif_certifie') is None:
                verdict = "indetermine : Hoppe complet non calculable"
            elif L.get('surjectif_certifie') is False:
                verdict = "indetermine : surjectivite de f non certifiee"
            elif L.get('indetermine'):
                verdict = "indetermine (w2V non calculable)"
            elif L['h0_generique'] != 0:
                verdict = "deja non stable avec f generique"
            elif L['h0_equivariant'] != 0:
                verdict = "tue par h0(V) equivariant"
            elif L.get('h0w2_generique'):
                verdict = "deja non stable : h0(w2V) generique != 0"
            else:
                verdict = "tue par h0(w2V) equivariant"
            survivants += bool(L['survit'])
            indetermines += bool(L.get('indetermine'))
            fmt = lambda x: '-' if x is None else str(x)
            print(f"  {r['cicy']:>5} {r.get('gauge', ''):>7} {L['rang_V']:>2} "
                  f"{L['groupe']:<11} {str(L['lambda']):>14} "
                  f"{L['dim_equivariant']:>6} {L['dim_totale']:>5} "
                  f"{L['h0_generique']:>6} {L['h0_equivariant']:>5} "
                  f"{fmt(L['h0w2_generique']):>6} "
                  f"{fmt(L['h0w2_equivariant']):>5}  {verdict}")
            sortie.append({**ident, **L})

    dst = os.path.join(args.output_dir, 'results_equivariance_f.jsonl')
    with open(dst, 'w', encoding='utf-8') as fh:
        for x in sortie:
            fh.write(json.dumps(x, default=int) + '\n')
    print(f"\n  Couples (candidat, lambda) qui survivent   : {survivants}")
    print(f"  Indetermines (un test non calculable)     : {indetermines}"
          f"   <- ni retenus ni elimines")
    print(f"  Ecartes avant evaluation                  : {ecartes}"
          f"   <- hors domaine, charges non permutees, etc.")
    print(f"  Toutes ces lignes sont dans le JSONL, champ 'etat'.")
    print(f"  Ecrit : {dst}")
    print(f"{'=' * 96}\n")


if __name__ == '__main__':
    sys.exit(main())
