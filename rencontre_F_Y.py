#!/usr/bin/env python3
"""
rencontre_F_Y.py -- LE TEMOIN DU §5.37 EST-IL SUR Y ?

LE TROU
-------
`lieu_de_base_rv3.py` exhibe un point ou les quatre f_i s'annulent, en fixant
les coordonnees des TROIS facteurs porteurs et en donnant aux autres une valeur
arbitraire. Cette liberte est legitime -- les f_i sont de degre 0 sur les
facteurs non porteurs, donc n'en dependent pas -- mais elle a une consequence
qui n'a jamais ete verifiee : les f_i s'annulent sur toute la sous-variete

    F = {p_0} x (produit des facteurs NON porteurs)

et un point de base doit etre SUR Y. Si F ne rencontre pas Y, le temoin ne
temoigne de rien. Le §5.37 conclut sans avoir pose cette question.

Elle a une seconde consequence, moins visible : les f_i sont des elements de
S/I, et la valeur d'un representant en un point HORS de Y depend du
representant choisi. La verification par substitution du §5.37 -- « les quatre
f_i valent zero » -- n'est bien definie que si le point est sur Y.

L'ARGUMENT, ET DANS QUEL SENS IL CONCLUT
----------------------------------------
Y est coupee dans l'ambiant par K polynomes. F ∩ Y est le lieu des zeros dans
F des K equations RESTREINTES. Si ce lieu est VIDE, la section correspondante
de L_1 (+) ... (+) L_K ne s'annule nulle part sur F, donc

    c_K( (+) L_i ) = prod c_1(L_i) = 0   dans H^{2K}(F),

et en particulier son accouplement avec H^{dim F - K} est nul (H ample sur F).
Par contraposition :

    prod_i [d_i] . H^{dim F - K} > 0   ==>   F ∩ Y ≠ vide.

C'est un critere SUFFISANT de RENCONTRE, et il est declare comme tel : un
nombre nul ne demontre pas que F rate Y. Le sens est le bon -- il faut prouver
la rencontre pour valider une elimination, et c'est ce que le nombre positif
donne.

L'anneau de Chow de F = prod P^{n_k} est Z[h_k]/(h_k^{n_k+1}) ; la classe de la
j-eme equation restreinte est somme_k cfg[j][k] h_k, ou k parcourt les facteurs
LIBRES. Le nombre cherche est le coefficient de prod h_k^{n_k}.

Usage :
    python -u rencontre_F_Y.py cicylist.txt [--source tous_indetermines.jsonl]
"""
import os
import sys
import json
import argparse
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def nombre_intersection(dims, classes):
    """
    Nombre d'intersection dans F = prod P^{dims[k]} des classes donnees,
    completees par H = somme h_k jusqu'a la dimension. Rend None si les
    classes sont trop nombreuses pour F (le critere ne s'applique pas).
    """
    r = len(dims)
    dimF = sum(dims)
    K = len(classes)
    if K > dimF:
        return None
    lineaires = [list(c) for c in classes] + [[1] * r] * (dimF - K)
    poly = {(0,) * r: 1}
    for lin in lineaires:
        suivant = {}
        for e, cf in poly.items():
            for k in range(r):
                if not lin[k] or e[k] + 1 > dims[k]:
                    continue
                e2 = list(e)
                e2[k] += 1
                e2 = tuple(e2)
                suivant[e2] = suivant.get(e2, 0) + cf * lin[k]
        poly = suivant
        if not poly:
            return 0
    return poly.get(tuple(dims), 0)


def _controles():
    """
    Les deux cotes du filtre, sur une geometrie ou la reponse est connue
    sans calcul.

    Dans P^1 x P^1 : deux diviseurs de classe (1,0) sont deux fibres de la
    premiere projection, DISJOINTS si distincts -- le nombre doit valoir 0.
    Un (1,0) et un (0,1) se coupent en UN point -- le nombre doit valoir 1.
    Un detecteur qui rendrait toujours un nombre positif echoue le premier ;
    un detecteur muet echoue le second.
    """
    assert nombre_intersection([1, 1], [[1, 0], [1, 0]]) == 0, 'controle negatif'
    assert nombre_intersection([1, 1], [[1, 0], [0, 1]]) == 1, 'controle positif'
    assert nombre_intersection([1, 1, 3],
                               [[1, 1, 0], [0, 0, 1], [0, 0, 1],
                                [0, 0, 1], [1, 1, 1]]) == 2, 'ancre #4078'
    # trop d'equations pour F : le critere ne dit rien, et le DIT
    assert nombre_intersection([1], [[1], [1], [1]]) is None, 'hors portee'


def rencontre(amb, cfg, porteurs):
    """
    Rend le verdict pour F = produit des facteurs non porteurs.
    cfg[j][k] = degre du j-eme polynome sur le k-eme facteur.
    """
    m = len(amb)
    libres = [k for k in range(m) if k not in porteurs]
    dims = [amb[k] for k in libres]
    classes = [[int(cfg[j][k]) for k in libres] for j in range(len(cfg))]
    inertes = [j for j, cl in enumerate(classes) if not any(cl)]
    n = nombre_intersection(dims, classes)
    return {'libres': libres, 'dims': dims, 'dim_F': sum(dims),
            'K': len(classes), 'classes': classes,
            'equations_inertes': inertes, 'nombre': n}


def principal():
    ap = argparse.ArgumentParser()
    ap.add_argument('cicylist')
    ap.add_argument('--source', default='tous_indetermines.jsonl')
    ap.add_argument('--sortie', default='rencontre_F_Y.jsonl')
    args = ap.parse_args()

    _controles()
    print("  controles internes : 2 cotes + ancre + hors-portee  ->  OK\n")

    from cy_landscape.data.parse_oxford import load_oxford_file
    E = {x['num']: x for x in load_oxford_file(args.cicylist)}

    # orientation de config : somme_j cfg[j][k] doit valoir n_k + 1 (c_1 = 0)
    e0 = E[4078]
    cfg0 = e0['config']
    somme = [sum(cfg0[j][k] for j in range(len(cfg0)))
             for k in range(len(e0['ambient']))]
    attendu = [n + 1 for n in e0['ambient']]
    print(f"  orientation de config : sommes {somme} vs n_k+1 {attendu}"
          f"   -> {'cfg[j][k]' if somme == attendu else 'A VERIFIER'}\n")
    if somme != attendu:
        print("  ARRET : l'orientation de config n'est pas celle supposee.")
        return

    cand = [d for d in map(json.loads, open(args.source, encoding='utf-8'))
            if d.get('etat') == 'ok' and d.get('rank_C') == 1
            and d.get('rang_V') == 3]
    print(f"  {len(cand)} candidats dans la strate rank_C=1 / rang_V=3\n")

    bilan = Counter()
    formes = Counter()
    with open(args.sortie, 'w', encoding='utf-8') as fh:
        for d in cand:
            b, c = json.loads(d['b']), json.loads(d['c'])
            e = E[d['cicy']]
            amb, cfg = e['ambient'], e['config']
            port = [k for k in range(len(amb))
                    if any(x[k] for x in list(b) + list(c))]
            r = rencontre(amb, cfg, port)
            r.update({'cicy': d['cicy'], 'groupe': d['groupe'],
                      'ambient': list(amb), 'porteurs': port})
            fh.write(json.dumps(r) + '\n')
            n = r['nombre']
            if n is None:
                cle = 'HORS PORTEE (dim F < K)'
            elif n > 0:
                cle = 'RENCONTRE DEMONTREE (nombre > 0)'
            elif r['equations_inertes']:
                cle = 'NUL, et des equations sont inertes sur F'
            else:
                cle = 'NUL : rien de demontre'
            bilan[cle] += 1
            formes[(tuple(amb), tuple(port), n)] += 1

    print(f"{'=' * 70}\n  BILAN, LES DEUX COTES\n{'=' * 70}")
    for k, v in bilan.most_common():
        print(f"    {v:>5}  {k}")
    print(f"\n  formes distinctes ({len(formes)}) :")
    for (amb, port, n), v in formes.most_common():
        print(f"    x{v:<4} ambient {list(amb)}  porteurs {list(port)}"
              f"  ->  F.Y = {n}")
    print("\n  Un nombre > 0 DEMONTRE que F rencontre Y, donc que le temoin du")
    print("  §5.37 est sur Y. Un nombre nul ne demontre pas le contraire : il")
    print("  laisse le temoin sans statut, et le candidat indetermine.")
    print(f"\n  detail : {args.sortie}")


if __name__ == '__main__':
    principal()
