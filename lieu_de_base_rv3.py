#!/usr/bin/env python3
"""
lieu_de_base_rv3.py -- f a-t-elle un lieu de base, sur la strate
rank_C = 1 / rang_V = 3 ? Reponse EXACTE, pas un certificat.

CE QUE CE SCRIPT DECIDE
-----------------------
Les 472 candidats de cette strate passent tous le critere de Hoppe complet
(h^0(V) = h^0(wedge^2 V) = 0, mesure non vide : sources 10 et 33) et aucun ne
voit sa surjectivite certifiee. Or `f_sans_point_base` est un critere
SUFFISANT : son echec ne demontre rien. Deux lectures s'opposaient — la liste
de degres s'arrete trop tot, ou f a reellement un lieu de base — et elles
donnent des conclusions inverses.

La strate se laisse trancher exactement, parce qu'elle est UNE SEULE
configuration repetee. Mesure sur les 472 : trois formes, identiques a
permutation pres des facteurs,

    b = 3 x O(e_k) + O(e_i + e_j),   c = O(3 e_k + e_i + e_j),
    tous les facteurs porteurs sont des P^1, et l'ideal est INERTE
    (dim (S/I)_d = dim S_d sur toutes les charges en jeu).

Les CICYs hebergeantes ne participent donc pas au calcul : le probleme vit sur
P^1 x P^1 x P^1.

L'ARGUMENT
----------
Avec z la coordonnee du facteur k et (x, y) celles des deux autres :

    f_4 est de degre (0,0,3) : un CUBIQUE BINAIRE en z, donc 3 racines.
    f_1, f_2, f_3 sont de degre (1,1,2).

Un zero commun exige f_4 = 0, donc z = z_r pour une des trois racines. En
z_r, chaque f_i devient une forme (1,1) sur P^1 x P^1, c'est-a-dire une
matrice 2x2 M_i. Les trois formes ont un zero commun si et seulement si

    il existe x != 0 tel que rang [ x^T M_1 ; x^T M_2 ; x^T M_3 ] <= 1

(car il faut alors un y != 0 dans le noyau commun). Les trois mineurs 2x2 de
cette matrice sont des QUADRATIQUES BINAIRES en x : leur pgcd est non trivial
si et seulement si elles ont une racine commune.

Tout est exact et instantane. Aucune recherche de degre, aucun plafond.

CE QUE LA REPONSE VAUT, ET DANS QUEL SENS
-----------------------------------------
TROUVER un point de base sur GF(p) le DEMONTRE sur GF(p) : c'est une
demonstration positive, dans la direction que le certificat J_d = R_d ne peut
pas fournir. N'en trouver aucun sur GF(p) ne demontre rien en caracteristique
nulle -- meme reserve que partout ailleurs dans ce depot.

Les racines de f_4 peuvent vivre dans une extension de GF(p) : le script le
DIT au lieu de conclure. Un candidat dont le cubique n'a aucune racine dans
GF(p) est rendu `indetermine`, pas `sans lieu de base`.

LA GARDE F ∩ Y (obligatoire, ajoutee apres coup)
------------------------------------------------
Le temoin fixe les coordonnees des trois facteurs PORTEURS et laisse les autres
libres : les f_i s'annulent donc sur toute la sous-variete F = {p_0} x (facteurs
non porteurs), et non en un point. Or un point de base doit etre SUR Y. Tant que
F ∩ Y n'est pas demontree non vide, le temoin ne temoigne de rien -- et la
substitution elle-meme n'est pas bien definie, les f_i etant des elements de S/I
dont un representant n'a de valeur intrinseque qu'en un point de Y.

`rencontre_F_Y.nombre_intersection` en fournit la preuve quand elle existe :
prod_i [d_i] . H^{dim F - K} > 0 dans l'anneau de Chow de F entraine F ∩ Y ≠ vide.
Aucun verdict `lieu_de_base = True` n'est rendu sans ce nombre strictement
positif ; un nombre nul rend le candidat `indetermine`, pas `elimine`.

Usage :
    python -u lieu_de_base_rv3.py cicyquotients.m cicylist.txt \
           [--source tous_indetermines.jsonl] [-n 40]
"""
import os
import sys
import json
import argparse
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from rencontre_F_Y import rencontre


def _pgcd_poly(a, b, p):
    """pgcd de deux polynomes en une variable sur GF(p), coefficients croissants."""
    def net(u):
        u = [int(t) % p for t in u]
        while u and u[-1] == 0:
            u.pop()
        return u
    a, b = net(a), net(b)
    while b:
        inv = pow(b[-1], p - 2, p)
        r = list(a)
        while len(r) >= len(b) and r:
            k = len(r) - len(b)
            f = (r[-1] * inv) % p
            for i, cb in enumerate(b):
                r[i + k] = (r[i + k] - f * cb) % p
            r = net(r)
        a, b = b, r
    return a


def _eval_mono(S, coeffs, valeurs, p):
    """
    Evalue le polynome (S, coeffs) au point donne.

    `basis_multi` indexe un monome PAR FACTEUR : c'est un tuple de tuples,
    un par facteur de l'ambiant, donnant les exposants des coordonnees de ce
    facteur. `valeurs` suit la meme forme -- une liste par facteur.
    """
    tot = 0
    for mono, cf in zip(S, coeffs):
        cf = int(cf) % p
        if not cf:
            continue
        t = cf
        for k, exps in enumerate(mono):
            vk = valeurs[k]
            for i, ex in enumerate(exps):
                if ex:
                    t = (t * pow(int(vk[i]) % p, int(ex), p)) % p
        tot = (tot + t) % p
    return tot % p


def analyser(anneau, amb, cfg, b, c, base, offsets, dims, degres, p, rng):
    """Rend un verdict sur le lieu de base, ou une raison de ne pas conclure."""
    from cy_landscape.core.equivariant_monad import _f_depuis_vecteur
    m = len(amb)
    port = [k for k in range(m) if any(x[k] for x in list(b) + list(c))]
    if len(port) != 3 or any(amb[k] != 1 for k in port):
        return {'etat': 'forme inattendue', 'porteurs': port}
    # l'indice k : celui ou une charge b vaut 1 trois fois
    compte = Counter()
    for x in b:
        for k in port:
            if x[k]:
                compte[k] += 1
    kz = [k for k, n in compte.items() if n == 3]
    if len(kz) != 1:
        return {'etat': 'forme inattendue (pas de facteur triple)'}
    kz = kz[0]
    kx, ky = [k for k in port if k != kz]

    # --- GARDE : le futur temoin sera-t-il seulement SUR Y ? --------------
    # Elle est placee AVANT tout calcul : un temoin qu'on ne pourra pas
    # situer sur Y ne vaut pas la peine d'etre cherche, et surtout ne doit
    # pas pouvoir etre rendu.
    rc = rencontre(amb, cfg, port)
    if rc['nombre'] is None or rc['nombre'] <= 0:
        return {'etat': 'F ne rencontre pas Y de facon demontree',
                'F_Y': rc['nombre'], 'dim_F': rc['dim_F'], 'K': rc['K'],
                'equations_inertes': rc['equations_inertes'],
                'porteurs': port}

    v = (rng.randint(0, p, size=base.shape[0]) @ base) % p
    fs = {}
    for i in range(len(b)):
        if (0, i) not in offsets:
            continue
        fs[i] = _f_depuis_vecteur(anneau, v, offsets, dims, degres, 0, i)
    if len(fs) != len(b):
        return {'etat': 'une case de f est absente'}
    i4 = [i for i in range(len(b)) if not b[i][kz]][0]   # celui de degre (0,0,3)
    autres = [i for i in range(len(b)) if i != i4]

    def point(vx, vy, vz):
        """Un point de l'ambiant, facteur par facteur. Les facteurs non
        porteurs n'apparaissent dans aucun monome (exposant nul partout) :
        la valeur qu'on leur donne est sans effet."""
        val = []
        for k in range(m):
            if k == kx:
                val.append(list(vx))
            elif k == ky:
                val.append(list(vy))
            elif k == kz:
                val.append(list(vz))
            else:
                val.append([1] + [0] * amb[k])
        return val

    # --- racines du cubique binaire f_4, sur P^1(GF(p)) -----------------
    deg4, S4, co4 = fs[i4]
    racines = []
    for t in range(p):
        if _eval_mono(S4, co4, point((1, 0), (1, 0), (1, t)), p) == 0:
            racines.append((1, t))
    if _eval_mono(S4, co4, point((1, 0), (1, 0), (0, 1)), p) == 0:
        racines.append((0, 1))
    if not racines:
        return {'etat': 'aucune racine de f_4 dans GF(p)',
                'n_racines': 0}

    # --- en chaque racine : les trois formes (1,1) ont-elles un zero ? ---
    def verifier(vx, vy, vz):
        """
        RECALCULE les quatre f_i au point exhibe. Un temoin qu'on affirme
        sans le substituer n'est pas un temoin : c'est la seule ligne de ce
        script qui ne depende d'aucun raisonnement en amont.
        """
        return {i: _eval_mono(fs[i][1], fs[i][2], point(vx, vy, vz), p)
                for i in fs}

    for vz in racines:
        M = []
        for i in autres:
            deg, S, co = fs[i]
            M.append([[_eval_mono(S, co, point((1, 0), (1, 0), vz), p),
                       _eval_mono(S, co, point((1, 0), (0, 1), vz), p)],
                      [_eval_mono(S, co, point((0, 1), (1, 0), vz), p),
                       _eval_mono(S, co, point((0, 1), (0, 1), vz), p)]])
        # x = (1, u) : a_i(u) = M_i[0][0] + u M_i[1][0], b_i(u) = M_i[0][1] + u M_i[1][1]
        A = [[Mi[0][0], Mi[1][0]] for Mi in M]
        B = [[Mi[0][1], Mi[1][1]] for Mi in M]

        def mul(u, w):
            r = [0] * (len(u) + len(w) - 1)
            for i, a in enumerate(u):
                for j, bb in enumerate(w):
                    r[i + j] = (r[i + j] + a * bb) % p
            return r

        quad = []
        for i in range(3):
            for j in range(i + 1, 3):
                q = [(x - y) % p for x, y in zip(mul(A[i], B[j]), mul(A[j], B[i]))]
                quad.append(q)
        g = quad[0]
        for q in quad[1:]:
            g = _pgcd_poly(g, q, p)

        # Candidats x : les racines du pgcd des trois mineurs, plus le point
        # a l'infini x = (0,1) qu'aucune parametrisation x = (1,u) n'atteint.
        cand_x = []
        if not g:                       # les trois mineurs sont identiquement nuls
            cand_x = [(1, u) for u in range(min(p, 4))] + [(0, 1)]
        elif len(g) > 1:
            cand_x = [(1, u) for u in range(p)
                      if sum(cf * pow(u, i, p)
                             for i, cf in enumerate(g)) % p == 0]
        cand_x.append((0, 1))

        for vx in cand_x:
            # r_i = x^T M_i. Un y non nul commun existe ssi les r_i non nuls
            # sont tous proportionnels, c'est-a-dire ssi le rang MOD P de la
            # matrice 3x2 vaut au plus 1. Un rang reel n'a rien a voir ici :
            # np.linalg.matrix_rank sur des residus est un faux ami.
            R = [[(vx[0] * Mi[0][0] + vx[1] * Mi[1][0]) % p,
                  (vx[0] * Mi[0][1] + vx[1] * Mi[1][1]) % p] for Mi in M]
            non_nuls = [r for r in R if r[0] % p or r[1] % p]
            if not non_nuls:
                vy = (1, 0)             # tout y convient
            else:
                r0 = non_nuls[0]
                if any((r0[0] * r[1] - r0[1] * r[0]) % p for r in non_nuls):
                    continue            # rang 2 : pas de y commun
                vy = ((-r0[1]) % p, r0[0] % p)   # noyau de r0
            vals = verifier(vx, vy, vz)
            return {'etat': 'ok', 'lieu_de_base': all(x == 0 for x in vals.values()),
                    'z': list(vz), 'x': list(vx), 'y': list(vy),
                    'F_Y': rc['nombre'], 'dim_F': rc['dim_F'], 'K': rc['K'],
                    'n_racines': len(racines),
                    'verification': {str(k): int(x) for k, x in vals.items()},
                    'temoin_verifie': all(x == 0 for x in vals.values())}
    return {'etat': 'ok', 'lieu_de_base': False, 'F_Y': rc['nombre'],
            'n_racines': len(racines)}


def principal():
    ap = argparse.ArgumentParser()
    ap.add_argument('braun_m')
    ap.add_argument('cicylist')
    ap.add_argument('--source', default='tous_indetermines.jsonl')
    ap.add_argument('-n', type=int, default=40)
    ap.add_argument('--sortie', default='lieu_de_base_rv3.jsonl')
    args = ap.parse_args()

    from cy_landscape.core.braun_symmetry import (parse_symmetries, ordres_rt,
                                                  matrice_mod_p)
    from cy_landscape.core.gamma_action import (choisir_premier,
                                                racine_primitive)
    from cy_landscape.core.covariant_ring import (resoudre_covariants,
                                                  tirer_covariants,
                                                  CovariantRing,
                                                  verifier_covariance)
    from cy_landscape.core.equivariant_monad import espace_f_equivariant
    from cy_landscape.data.parse_oxford import load_oxford_file
    from wilson_match import parse_braun, parse_cicylist, apparier

    E = {x['num']: x for x in load_oxford_file(args.cicylist)}
    SYM = parse_symmetries(args.braun_m)
    corr, _, _ = apparier(parse_braun(args.braun_m),
                          parse_cicylist(args.cicylist))
    inv = {v: k for k, v in corr.items()}

    cand = [d for d in map(json.loads, open(args.source, encoding='utf-8'))
            if d.get('etat') == 'ok' and d.get('rank_C') == 1
            and d.get('rang_V') == 3]
    print(f"  {len(cand)} candidats dans la strate ; on en traite {args.n}\n",
          flush=True)

    bilan = Counter()
    with open(args.sortie, 'w', encoding='utf-8') as fh:
        for n, d in enumerate(cand[:args.n], 1):
            b, c = json.loads(d['b']), json.loads(d['c'])
            e = E[d['cicy']]
            amb, cfg = e['ambient'], np.asarray(e['config'])
            syms = [s for s in SYM[inv[d['cicy']]]['symetries']
                    if s['nom'] == d['groupe']]
            if not syms:
                continue
            sym = syms[0]
            ordres = sorted(ordres_rt(sym['coord']) | ordres_rt(sym['poly'])
                            | {2})
            p, _ = choisir_premier(ordres, minimum=30011)
            rac = {k: racine_primitive(p, k) for k in ordres}
            Mc = [matrice_mod_p(x, p, rac) for x in sym['coord']]
            Np = [matrice_mod_p(x, p, rac) for x in sym['poly']]
            res = resoudre_covariants(amb, cfg, Mc, Np, p)
            if res is None:
                continue
            co = tirer_covariants(res['par_convention']['N']['base'],
                                  res['offsets'], res['dims'], p,
                                  np.random.RandomState(0))
            ok, _ = verifier_covariance(amb, cfg, Mc, Np, 'N', co, p)
            if not ok:
                continue
            A = CovariantRing(amb, cfg, co, p)
            out = espace_f_equivariant(A, amb, b, c, Mc, p)
            if out['etat'] != 'ok' or not out['solutions']:
                continue
            for s in out['solutions']:
                r = analyser(A, amb, cfg, b, c, s['base'], out['offsets'],
                             out['dims'], out['degres'], p,
                             np.random.RandomState(7))
                r.update({'cicy': d['cicy'], 'groupe': d['groupe']})
                fh.write(json.dumps(r, default=str) + '\n')
                fh.flush()
                if r.get('etat') != 'ok':
                    cle = r['etat']
                elif r['lieu_de_base']:
                    cle = f"LIEU DE BASE (temoin verifie, sur Y : F.Y = {r['F_Y']})"
                elif r.get('x') is not None:
                    # un temoin propose mais qui ne s'annule pas : c'est le
                    # raisonnement qui est faux, pas la geometrie.
                    cle = 'TEMOIN REFUSE PAR SUBSTITUTION'
                else:
                    cle = 'pas de lieu de base sur GF(p)'
                bilan[cle] += 1
                print(f"  [{n}] #{d['cicy']:<5} {d['groupe']:<9} "
                      f"{cle}"
                      + (f"   z={r.get('z')} x={r.get('x')} y={r.get('y')}"
                         if r.get('x') is not None else "")
                      + f"   ({r.get('n_racines', '?')} racines de f_4 "
                        f"dans GF(p))", flush=True)

    print(f"\n{'=' * 70}")
    print("  BILAN, LES DEUX COTES")
    print(f"{'=' * 70}")
    for k, v in bilan.most_common():
        print(f"    {v:>5}  {k}")
    print("\n  Trouver un point de base le DEMONTRE (sur GF(p)). N'en trouver")
    print("  aucun ne demontre rien : le certificat J_d = R_d ne concluait pas")
    print("  non plus, et pour la meme raison de sens.")
    print(f"\n  detail : {args.sortie}")


if __name__ == '__main__':
    principal()
