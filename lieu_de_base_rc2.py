#!/usr/bin/env python3
"""
lieu_de_base_rc2.py -- les 34 dernieres lignes lambda indeterminees.

CE QU'IL RESTE, ET POURQUOI C'EST UNE AUTRE GEOMETRIE
-----------------------------------------------------
Apres le §5.37, l'indetermine se reduit a 34 lignes lambda : celles qui passent
le critere de Hoppe complet a rank_C = 2 SANS certificat de surjectivite. Le
§5.37 annonce qu'il faudra « l'equivalent de lieu_de_base_rv3.py, ou le lieu de
base est celui des mineurs 2x2 ». C'est vrai, et insuffisant : le decompte de
dimensions n'est pas le meme, et c'est LUI qui decide.

La strate est deux formes, a permutation pres des facteurs :

    porteurs P^1 x P^n (n = 3 pour 28 lignes, n = 4 pour 6)
    b = 2 x O(0,1) + 3 x O(1,0)      c = O(1,1) (+) O(2,1)

Avec x sur le P^1 et y sur le P^n, f est 2 x 5 :

    ligne 0 :  L_0(x)  L_1(x)  |  A_2(y)  A_3(y)  A_4(y)
    ligne 1 :  Q_0(x)  Q_1(x)  |  B_2(x,y) B_3(x,y) B_4(x,y)

L_i lineaire, Q_i quadratique en x, A_j lineaire en y, B_j bilineaire.

f est surjective en un point ssi la matrice y est de rang 2. Le lieu de base
est donc {rang <= 1}, et non le lieu des zeros communs des f_i.

L'EXISTENCE NE FAIT PAS DE DOUTE -- C'EST LA DIMENSION QUI TRANCHE
-------------------------------------------------------------------
det[[L_0,L_1],[Q_0,Q_1]] est un CUBIQUE BINAIRE en x : il a trois racines. En
une racine x*, les deux premieres colonnes deviennent proportionnelles ; si
elles ne sont pas nulles, le rapport mu est determine, et rang <= 1 equivaut
aux TROIS formes LINEAIRES B_j(x*,y) - mu A_j(y) = 0 sur P^n. Trois formes
lineaires sur P^n avec n >= 3 ont toujours un zero non nul. Le lieu de base
n'est donc JAMAIS vide dans P^1 x P^n -- pour f equivariante comme pour f
generique. Ce n'est pas un discriminant.

Ce qui discrimine, c'est ceci. Le lieu de base dans l'ambiant est

    F = {x*} x Lambda x (facteurs non porteurs),    dim Lambda = n - rang,

et un point de base doit etre SUR Y, qui est de codimension K. Or

    dim F = dim Lambda + somme(libres),   K = 1 + n + somme(libres) - 3.

Si les trois formes sont INDEPENDANTES, dim Lambda = n - 3 et dim F = K - 1 :
une dimension de MOINS que la codimension de Y. F rate Y generiquement, et le
critere d'Euler ne s'applique meme pas -- on ne peut rien conclure, ni dans un
sens ni dans l'autre. Si l'equivariance fait CHUTER le rang, dim F >= K et la
rencontre redevient demontrable.

C'est le contraire du §5.37, ou les trois facteurs porteurs etaient des P^1 et
ou dim F = K tombait juste sur les 472. Ce script MESURE le rang au lieu de le
supposer, et rend, pour chaque lambda, de quel cote elle tombe.

CONTROLE, LES DEUX COTES
------------------------
Le meme calcul est fait sur f EQUIVARIANTE et sur f GENERIQUE dans le meme
anneau. Le rang generique est la reference : si l'equivariance ne le fait pas
chuter, elle n'impose rien ici, et le dire est un resultat.

Usage :
    python -u lieu_de_base_rc2.py cicyquotients.m cicylist.txt \
           [--source tous_indetermines.jsonl] [-n 40]
"""
import os
import sys
import json
import argparse
from itertools import combinations
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from rencontre_F_Y import nombre_intersection
from lieu_de_base_rv3 import _eval_mono


def _rref(M, p):
    """Reduction de Gauss mod p. Rend (lignes reduites, colonnes pivots)."""
    M = [[int(x) % p for x in ligne] for ligne in M]
    nl = len(M)
    nc = len(M[0]) if nl else 0
    pivots = []
    r = 0
    for c in range(nc):
        piv = next((i for i in range(r, nl) if M[i][c] % p), None)
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = pow(M[r][c], p - 2, p)
        M[r] = [(x * inv) % p for x in M[r]]
        for i in range(nl):
            if i != r and M[i][c] % p:
                f = M[i][c]
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[r])]
        pivots.append(c)
        r += 1
        if r == nl:
            break
    return M, pivots


def _noyau(M, p, nc):
    """Une base du noyau de M (mod p), M ayant nc colonnes."""
    if not M:
        return [[1 if i == t else 0 for i in range(nc)] for t in range(nc)]
    R, piv = _rref(M, p)
    libres = [c for c in range(nc) if c not in piv]
    base = []
    for c in libres:
        v = [0] * nc
        v[c] = 1
        for i, pc in enumerate(piv):
            v[pc] = (-R[i][c]) % p
        base.append(v)
    return base


def _cubique(valeurs, p):
    """Coefficients a_0..a_3 de D(1,t) a partir de D en t = 0,1,2,3."""
    V = [[pow(t, i, p) for i in range(4)] + [valeurs[t]] for t in range(4)]
    R, piv = _rref(V, p)
    assert piv == [0, 1, 2, 3], 'Vandermonde degeneree : p est trop petit'
    return [R[i][4] % p for i in range(4)]


def analyser_rc2(anneau, amb, cfg, b, c, base, offsets, dims, degres, p, rng,
                 etiquette=''):
    """Rend le verdict de dimension pour un f tire dans `base`."""
    from cy_landscape.core.equivariant_monad import _f_depuis_vecteur
    m = len(amb)
    port = [k for k in range(m) if any(x[k] for x in list(b) + list(c))]
    if len(port) != 2:
        return {'etat': f'forme inattendue : {len(port)} facteurs porteurs'}
    kx = [k for k in port if amb[k] == 1]
    if len(kx) != 1:
        return {'etat': 'forme inattendue : pas un seul P^1 porteur'}
    kx = kx[0]
    ky = [k for k in port if k != kx][0]
    n = amb[ky]
    libres = [k for k in range(m) if k not in port]

    # les deux lignes de C, ordonnees par degre croissant sur le P^1
    if len(c) != 2:
        return {'etat': f'rank_C = {len(c)}'}
    j0, j1 = sorted(range(2), key=lambda j: c[j][kx])
    # colonnes dont la ligne 0 ne depend que de x
    Cx = [i for i in range(len(b)) if c[j0][ky] - b[i][ky] == 0]
    Cy = [i for i in range(len(b)) if i not in Cx]
    if len(Cx) != 2 or len(Cy) != 3:
        return {'etat': f'partition inattendue des colonnes : {len(Cx)}/{len(Cy)}'}

    v = (rng.randint(0, p, size=base.shape[0]) @ base) % p
    f = {}
    for j in range(2):
        for i in range(len(b)):
            if (j, i) not in offsets:
                return {'etat': f'case ({j},{i}) de f absente'}
            f[(j, i)] = _f_depuis_vecteur(anneau, v, offsets, dims, degres, j, i)

    def pt(vx, vy):
        val = []
        for k in range(m):
            if k == kx:
                val.append(list(vx))
            elif k == ky:
                val.append(list(vy))
            else:
                val.append([1] + [0] * amb[k])
        return val

    def ev(j, i, vx, vy):
        _, S, co = f[(j, i)]
        return _eval_mono(S, co, pt(vx, vy), p)

    y0 = [1] + [0] * n

    # --- le cubique binaire, par interpolation (4 evaluations) -----------
    vals = []
    for t in range(4):
        u = [ev(j0, i, (1, t), y0) for i in Cx]
        w = [ev(j1, i, (1, t), y0) for i in Cx]
        vals.append((u[0] * w[1] - u[1] * w[0]) % p)
    a = _cubique(vals, p)
    identiquement_nul = not any(a)

    racines = [(1, t) for t in range(p)
               if sum(a[i] * pow(t, i, p) for i in range(4)) % p == 0]
    if not a[3] % p:
        racines.append((0, 1))           # x = (0,1) : coefficient dominant
    if identiquement_nul:
        racines = [(1, 0), (0, 1)] + [(1, t) for t in range(1, 3)]

    if not racines:
        return {'etat': 'aucune racine du cubique dans GF(p)',
                'etiquette': etiquette, 'n_racines': 0}

    for vx in racines:
        u = [ev(j0, i, vx, y0) for i in Cx]
        w = [ev(j1, i, vx, y0) for i in Cx]
        if any(u):
            mu = None
            for t in range(2):
                if u[t] % p:
                    mu = (w[t] * pow(u[t], p - 2, p)) % p
                    break
            if any((w[t] - mu * u[t]) % p for t in range(2)):
                continue                 # bloc de rang 2 : pas ici
            cas = 'mu fini'
            def forme(i, vy):
                return (ev(j1, i, vx, vy) - mu * ev(j0, i, vx, vy)) % p
        elif any(w):
            cas = 'ligne 0 nulle sur le bloc'
            def forme(i, vy):
                return ev(j0, i, vx, vy)
        else:
            # bloc 2x2 nul : la condition devient QUADRATIQUE en y.
            return {'etat': 'bloc 2x2 identiquement nul : regime quadratique',
                    'etiquette': etiquette, 'x': list(vx)}

        # coefficients des trois formes lineaires : forme(i, e_t)
        L = []
        for i in Cy:
            L.append([forme(i, [1 if s == t else 0 for s in range(n + 1)])
                      for t in range(n + 1)])
        _, piv = _rref(L, p)
        rang = len(piv)
        dim_L = n - rang
        noyau = _noyau(L, p, n + 1)
        if not noyau:
            continue
        vy = noyau[0]

        # --- VERIFICATION : les dix mineurs 2x2 au point exhibe ---------
        Mnum = [[ev(jj, i, vx, vy) for i in range(len(b))] for jj in (j0, j1)]
        mineurs = {f'{i},{k}': (Mnum[0][i] * Mnum[1][k]
                                - Mnum[0][k] * Mnum[1][i]) % p
                   for i, k in combinations(range(len(b)), 2)}
        temoin = all(x == 0 for x in mineurs.values())

        # --- F rencontre-t-elle Y ? -------------------------------------
        dims_F = [dim_L] + [amb[k] for k in libres]
        classes = [[int(cfg[jj][ky])] + [int(cfg[jj][k]) for k in libres]
                   for jj in range(len(cfg))]
        if identiquement_nul:            # x reste libre : le P^1 est dans F
            dims_F = [1] + dims_F
            classes = [[int(cfg[jj][kx])] + cl
                       for jj, cl in enumerate(classes)]
        dims_F = [d for d in dims_F if d > 0] or [0]
        if 0 in dims_F:
            dims_F = [d for d in dims_F if d]
            classes = [cl for cl in classes]
        nb = nombre_intersection([d for d in dims_F if d > 0],
                                 [[cl[t] for t, d in enumerate(dims_F)
                                   if d > 0] for cl in classes]) \
            if any(d > 0 for d in dims_F) else None

        return {'etat': 'ok', 'etiquette': etiquette, 'cas': cas,
                'x': list(vx), 'y': list(vy), 'n_racines': len(racines),
                'cubique_nul': identiquement_nul,
                'rang_formes': rang, 'dim_Lambda': dim_L, 'n': n,
                'dim_F': sum(d for d in dims_F if d > 0), 'K': len(cfg),
                'F_Y': nb, 'temoin_rang_1': temoin,
                'mineurs': {k: int(x) for k, x in mineurs.items()}}
    return {'etat': 'aucune racine ne donne un bloc proportionnel',
            'etiquette': etiquette, 'n_racines': len(racines)}


def reserve(cand, E, SYM, inv, args):
    """
    LA RESERVE MOD P, OPPOSEE.

    Trouver un lieu de base mod p le DEMONTRE mod p, et la reserve joue donc
    CONTRE une elimination : c'est le sens defavorable. On la mesure au lieu de
    la supposer -- plusieurs premiers, plusieurs tirages de l'ideal covariant.
    Il faudrait qu'ils degenerent tous de la meme facon, sans toucher le calcul
    generique qui, lui, ne rend pas de lieu de base exploitable.
    """
    from cy_landscape.core.braun_symmetry import (parse_symmetries, ordres_rt,
                                                  matrice_mod_p)
    from cy_landscape.core.gamma_action import (choisir_premier,
                                                racine_primitive)
    from cy_landscape.core.covariant_ring import (resoudre_covariants,
                                                  tirer_covariants,
                                                  CovariantRing)
    from cy_landscape.core.equivariant_monad import espace_f_equivariant

    print(f"\n{'=' * 74}\n  RESERVE MOD P : {args.reserve} candidats "
          f"x 3 premiers x 2 tirages\n{'=' * 74}")
    bilan = Counter()
    for d, idx in cand[:args.reserve]:
        b, c = json.loads(d['b']), json.loads(d['c'])
        e = E[d['cicy']]
        amb, cfg = e['ambient'], np.asarray(e['config'])
        sym = [x for x in SYM[inv[d['cicy']]]['symetries']
               if x['nom'] == d['groupe']][0]
        ordres = sorted(ordres_rt(sym['coord']) | ordres_rt(sym['poly']) | {2})
        ligne = []
        for mini in (30011, 50021, 70003):
            p, _ = choisir_premier(ordres, minimum=mini)
            rac = {k: racine_primitive(p, k) for k in ordres}
            Mc = [matrice_mod_p(x, p, rac) for x in sym['coord']]
            Np = [matrice_mod_p(x, p, rac) for x in sym['poly']]
            res = resoudre_covariants(amb, cfg, Mc, Np, p)
            if res is None:
                ligne.append('sigma?')
                continue
            for graine in (0, 5):
                co = tirer_covariants(res['par_convention']['N']['base'],
                                      res['offsets'], res['dims'], p,
                                      np.random.RandomState(graine))
                A = CovariantRing(amb, cfg, co, p)
                out = espace_f_equivariant(A, amb, b, c, Mc, p)
                if out['etat'] != 'ok' or not out['solutions']:
                    ligne.append('espace?')
                    continue
                t = idx[0] if idx[0] < len(out['solutions']) else 0
                r = analyser_rc2(A, amb, cfg, b, c, out['solutions'][t]['base'],
                                 out['offsets'], out['dims'], out['degres'], p,
                                 np.random.RandomState(7))
                ok = (r.get('etat') == 'ok' and r.get('temoin_rang_1')
                      and (r.get('F_Y') or 0) > 0)
                ligne.append(f"p={p} t{graine}:"
                             + ('elimine' if ok else str(r.get('etat'))[:16]))
                bilan['elimine' if ok else 'NON concluant'] += 1
        print(f"  #{d['cicy']:<5} {d['groupe']:<6} " + "  ".join(ligne))
    print(f"\n  {dict(bilan)}")
    tot = sum(bilan.values())
    print(f"  {bilan['elimine']} / {tot} evaluations concluent a l'elimination.")


def principal():
    ap = argparse.ArgumentParser()
    ap.add_argument('braun_m')
    ap.add_argument('cicylist')
    ap.add_argument('--source', default='tous_indetermines.jsonl')
    ap.add_argument('-n', type=int, default=40)
    ap.add_argument('--sortie', default='lieu_de_base_rc2.jsonl')
    ap.add_argument('--reserve', type=int, default=0,
                    help='mesurer la reserve mod p sur N candidats')
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

    cand = []
    for d in map(json.loads, open(args.source, encoding='utf-8')):
        if d.get('etat') != 'ok' or d.get('rank_C') != 2:
            continue
        idx = [t for t, mz in enumerate(d.get('mesures', []))
               if mz.get('hoppe') is True and not mz.get('surjectif')]
        if idx:
            cand.append((d, idx))
    total = sum(len(i) for _, i in cand)
    print(f"  {len(cand)} candidats, {total} lignes lambda en attente ; "
          f"on en traite {min(args.n, len(cand))} candidats\n", flush=True)

    if args.reserve:
        reserve(cand, E, SYM, inv, args)
        return

    bilan = Counter()
    with open(args.sortie, 'w', encoding='utf-8') as fh:
        for n_, (d, idx) in enumerate(cand[:args.n], 1):
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
            plein = np.eye(out['dim_totale'], dtype=np.int64)
            # Un cubique generique n'a pas toujours de racine dans GF(p) : le
            # volet de controle serait alors MUET, et un volet muet ne controle
            # rien. On insiste sur plusieurs tirages avant de le declarer tel.
            gen = None
            for graine in (11, 13, 17, 19, 23, 29, 31):
                gen = analyser_rc2(A, amb, cfg, b, c, plein, out['offsets'],
                                   out['dims'], out['degres'], p,
                                   np.random.RandomState(graine), 'generique')
                gen['graine'] = graine
                if gen.get('etat') == 'ok':
                    break
            for t in idx:
                if t >= len(out['solutions']):
                    continue
                s = out['solutions'][t]
                r = analyser_rc2(A, amb, cfg, b, c, s['base'], out['offsets'],
                                 out['dims'], out['degres'], p,
                                 np.random.RandomState(7), f'lambda[{t}]')
                r.update({'cicy': d['cicy'], 'groupe': d['groupe'],
                          'lambda': t, 'dim_espace': int(s['dim']),
                          'dim_totale': int(out['dim_totale']),
                          'generique': gen})
                fh.write(json.dumps(r, default=str) + '\n')
                fh.flush()
                if r.get('etat') != 'ok':
                    cle = r['etat']
                elif not r['temoin_rang_1']:
                    cle = 'TEMOIN REFUSE PAR SUBSTITUTION'
                elif r['F_Y'] is None:
                    cle = 'dim F < K : le critere ne s applique pas'
                elif r['F_Y'] > 0:
                    cle = 'LIEU DE BASE SUR Y (elimine)'
                else:
                    cle = 'F.Y = 0 : rien de demontre'
                bilan[cle] += 1
                rg = r.get('rang_formes')
                rgg = gen.get('rang_formes')
                print(f"  [{n_}] #{d['cicy']:<5} {d['groupe']:<6} "
                      f"l[{t}] dim {s['dim']}/{out['dim_totale']}  "
                      f"rang equiv {rg} / generique {rgg}  "
                      f"dim F {r.get('dim_F')} vs K {r.get('K')}  "
                      f"F.Y = {r.get('F_Y')}   {cle}", flush=True)

    print(f"\n{'=' * 74}\n  BILAN, LES DEUX COTES\n{'=' * 74}")
    for k, v in bilan.most_common():
        print(f"    {v:>5}  {k}")
    print("\n  Un lieu de base DEMONTRE sur Y elimine. Un `dim F < K` n'elimine")
    print("  pas et ne sauve pas : il dit que ce chemin ne peut pas trancher,")
    print("  et que le certificat de surjectivite reste la seule voie.")
    print(f"\n  detail : {args.sortie}")


if __name__ == '__main__':
    principal()
