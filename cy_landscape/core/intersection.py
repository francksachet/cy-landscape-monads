"""
intersection.py -- Calcul EXACT des nombres d'intersection triples
et des classes de Chern a partir de la matrice de configuration CICY.

Pour un CICY Y dans P^{n1} x ... x P^{nm} defini par K contraintes
de multi-degres q_{a,r}, le nombre d'intersection triple est :

  d_{ijk} = coeff de (J1^n1 ... Jm^nm) dans
            Ji * Jj * Jk * prod_{a=1}^K (sum_r q_{ar} Jr)

Ceci est un calcul polynomial exact -- aucune approximation.
"""

import numpy as np
from itertools import combinations
from math import comb
from typing import List, Tuple, Dict


# ================================================================
# Arithmetique polynomiale sparse sur Z[J1,...,Jm]
# Un polynome est un dict: {(e1,...,em): coeff}
# ================================================================

def _mono(m: int, idx: int) -> dict:
    """Monome J_idx dans m variables."""
    exp = [0] * m
    exp[idx] = 1
    return {tuple(exp): 1}


def _const(m: int, val: int) -> dict:
    """Constante dans m variables."""
    if val == 0:
        return {}
    return {tuple([0] * m): val}


def _poly_add(p1: dict, p2: dict) -> dict:
    result = dict(p1)
    for exp, coeff in p2.items():
        result[exp] = result.get(exp, 0) + coeff
        if result[exp] == 0:
            del result[exp]
    return result


def _poly_mul(p1: dict, p2: dict) -> dict:
    result = {}
    for exp1, c1 in p1.items():
        for exp2, c2 in p2.items():
            exp = tuple(e1 + e2 for e1, e2 in zip(exp1, exp2))
            result[exp] = result.get(exp, 0) + c1 * c2
    # Nettoyer les zeros
    return {k: v for k, v in result.items() if v != 0}


def _poly_coeff(poly: dict, target_exp: tuple) -> int:
    """Coefficient du monome target_exp dans le polynome."""
    return poly.get(target_exp, 0)


# ================================================================
# Nombres d'intersection triples
# ================================================================

def compute_intersection_numbers(
    ambient_dims: List[int],
    config_matrix: np.ndarray,
) -> np.ndarray:
    """
    Calcule les nombres d'intersection triples d_{ijk} d'un CICY.

    d_{ijk} = coeff de prod_r Jr^{nr} dans
              Ji Jj Jk prod_{a=1}^K (sum_r q_{ar} Jr)

    Args:
        ambient_dims: [n1, ..., nm] dimensions des espaces projectifs
        config_matrix: matrice K x m des multi-degres

    Returns:
        Tenseur symetrique d_{ijk} de shape (m, m, m)
    """
    m = len(ambient_dims)
    K = config_matrix.shape[0] if config_matrix.ndim > 1 else 1
    if config_matrix.ndim == 1:
        config_matrix = config_matrix.reshape(1, -1)

    # Monome cible : J1^n1 ... Jm^nm
    target = tuple(ambient_dims)

    # Precalculer les polynomes de contrainte : sum_r q_{ar} Jr
    constraint_polys = []
    for a in range(K):
        poly = {}
        for r in range(m):
            if config_matrix[a, r] != 0:
                poly = _poly_add(poly, {tuple(
                    1 if j == r else 0 for j in range(m)
                ): int(config_matrix[a, r])})
        constraint_polys.append(poly)

    # Degre max utile : monome cible a degre total = sum(ambient_dims).
    # Les monomes de degre superieur ne contribueront jamais au coefficient cible,
    # donc on peut les jeter apres chaque multiplication (garde-fou anti-explosion).
    # On ajoute 3 pour tenir compte du Ji*Jj*Jk qu'on multipliera plus tard.
    max_deg = sum(ambient_dims)

    def _truncate(poly, max_d):
        return {e: c for e, c in poly.items() if sum(e) <= max_d}

    # Produit de toutes les contraintes (ne depend pas de i,j,k)
    # Le produit final doit avoir des monomes de degre sum(ambient_dims) - 3
    # (car on multipliera par Ji*Jj*Jk plus tard, de degre 3).
    constraint_max_deg = max_deg - 3
    constraint_product = _const(m, 1)
    for poly in constraint_polys:
        constraint_product = _poly_mul(constraint_product, poly)
        constraint_product = _truncate(constraint_product, constraint_max_deg)

    # Calculer d_{ijk} pour chaque triplet
    d = np.zeros((m, m, m), dtype=float)

    for i in range(m):
        for j in range(i, m):
            for k in range(j, m):
                # Construire Ji * Jj * Jk
                ijk_exp = [0] * m
                ijk_exp[i] += 1
                ijk_exp[j] += 1
                ijk_exp[k] += 1
                ijk_poly = {tuple(ijk_exp): 1}

                # Multiplier par le produit des contraintes
                full = _poly_mul(ijk_poly, constraint_product)

                # Extraire le coefficient cible
                val = _poly_coeff(full, target)

                # Symetriser
                d[i, j, k] = val
                d[i, k, j] = val
                d[j, i, k] = val
                d[j, k, i] = val
                d[k, i, j] = val
                d[k, j, i] = val

    return d


# ================================================================
# Caracteristique d'Euler
# ================================================================

def compute_euler_from_intersection(
    ambient_dims: List[int],
    config_matrix: np.ndarray,
    d_ijk: np.ndarray,
) -> int:
    """
    Calcule la caracteristique d'Euler chi(Y) = integral c3(TY).

    Utilise la formule d'adjonction :
      c(TY) = c(T_ambient) / c(N)
    avec c(T_{P^n}) = (1+J)^{n+1} et N = oplus O(q_a).
    """
    m = len(ambient_dims)
    K = config_matrix.shape[0]

    # c3(TY) = c3(T_amb) - c3(N) - c1(T_amb)*c2(N) + c1(N)*c2(T_amb)
    #        + c1(N)*c2(N) - c1(N)^3 + ... (termes de la formule de Whitney)
    # Puisque c1(TY) = 0, on a c1(T_amb) = c1(N), et la formule se simplifie.

    # Methode: calculer directement via le polynome de Chern total
    # c(TY) = prod_r (1+Jr)^{nr+1} / prod_a (1 + sum_r qar Jr)
    # et integrer c3 sur Y (en utilisant les d_ijk)

    # Calculer les coefficients du polynome c(TY) modulo degre > 3
    # en utilisant l'expansion en series

    # c_ambient = prod (1+Jr)^{nr+1}
    # Expandre en composantes de Chern de chaque degre

    # Degre 1: c1_amb = sum_r (nr+1) Jr
    c1_amb = np.array([n + 1 for n in ambient_dims], dtype=float)

    # Degre 1 du normal: c1_N = sum_a sum_r qar Jr = sum_r (sum_a qar) Jr
    c1_N = np.sum(config_matrix, axis=0).astype(float)

    # Verifier CY: c1_amb = c1_N
    assert np.allclose(c1_amb, c1_N), f"Pas CY: c1_amb={c1_amb}, c1_N={c1_N}"

    # Degre 2: c2_amb[r,s] = coeff de Jr Js dans prod (1+Jr)^{nr+1}
    c2_amb = np.zeros((m, m))
    for r in range(m):
        c2_amb[r, r] = comb(ambient_dims[r] + 1, 2)
        for s in range(r + 1, m):
            c2_amb[r, s] = (ambient_dims[r] + 1) * (ambient_dims[s] + 1)
            c2_amb[s, r] = c2_amb[r, s]

    # c2_N[r,s] = sum_{a<b} qar qbs + qas qbr (termes croises des contraintes)
    c2_N = np.zeros((m, m))
    for a in range(K):
        for b in range(a + 1, K):
            for r in range(m):
                for s in range(m):
                    c2_N[r, s] += config_matrix[a, r] * config_matrix[b, s]

    c2_TY = c2_amb - c2_N  # Car c1(TY) = 0

    # Degre 3: c3 = c3_amb - c3_N - c1*c2_N + c1*c2_amb + ...
    # Plus simple: utiliser chi = integral c3
    # chi = sum_{r,s,t} c3_TY[r,s,t] * ... integrales

    # Methode directe: utiliser la formule de HRR pour chi
    # chi(O_Y) = integral Td3(Y) = chi/... non, c'est circulaire

    # Methode polynomiale exacte pour chi:
    # chi(Y) = coeff de prod Jr^nr dans
    #   prod_r (Jr/(1-e^{-Jr}))^{nr+1} * prod_r (1-e^{-Jr})^{nr+1}
    #   / prod_a (1 + sum qar Jr + ...)
    # C'est complique. Utilisons plutot la formule directe.

    # En fait, chi = 2(h11 - h21) et nous connaissons h11 = m (favorable)
    # Mais on veut calculer chi independamment...

    # Formule de Libgober-Wood pour les CICY:
    # chi = sum d_{ijk} c3_coeff_{ijk}
    # ou c3 est le coefficient de degre 3 de ln(c(TY))

    # Approche la plus simple et correcte:
    # chi(Y) = integral_Y c_3(TY)
    # c_3(TY) vit dans H^6(Y) ~ Z
    # On peut la calculer via le polynome:
    # c(TY) = prod_r (1+Jr)^{nr+1} * prod_a 1/(1+sum qar Jr)
    # et c3 est le terme de degre 3 de ce developpement

    # Expandre 1/(1+x) = 1 - x + x^2 - x^3 + ... modulo degre 4

    # Terme de degre 3 de c(TY):
    # = c3_amb - (c2_amb * c1_N_terms) + (c1_amb * c2_N_quadratic) - c3_N + ...

    # C'est devenu trop complique en composantes. Utilisons le calcul polynomial.

    # Construire le polynome c(TY) modulo degre > 3, en variables J1,...,Jm
    # c_amb = prod_r sum_{k=0}^{nr+1} C(nr+1,k) Jr^k
    c_amb = _const(m, 1)
    for r in range(m):
        # (1 + Jr)^{nr+1} tronque au degre 3
        factor = {}
        for k in range(min(ambient_dims[r] + 2, 4)):  # Degre max 3
            exp = [0] * m
            exp[r] = k
            factor[tuple(exp)] = comb(ambient_dims[r] + 1, k)
        c_amb = _poly_mul(c_amb, factor)
        # Tronquer au degre 3
        c_amb = {e: c for e, c in c_amb.items() if sum(e) <= 3}

    # c_N_inv = prod_a 1/(1 + sum_r qar Jr)
    # = prod_a (1 - (sum qar Jr) + (sum qar Jr)^2 - ...) tronque deg 3
    c_N_inv = _const(m, 1)
    for a in range(K):
        # x_a = sum_r qar Jr
        x_a = {}
        for r in range(m):
            if config_matrix[a, r] != 0:
                exp = [0] * m
                exp[r] = 1
                x_a[tuple(exp)] = int(config_matrix[a, r])

        # 1/(1+x) ~ 1 - x + x^2 - x^3 modulo deg 4
        inv_factor = _const(m, 1)
        x_power = _const(m, 1)
        for k in range(1, 4):
            x_power = _poly_mul(x_power, x_a)
            x_power = {e: c for e, c in x_power.items() if sum(e) <= 3}
            sign = (-1) ** k
            inv_factor = _poly_add(inv_factor,
                                    {e: sign * c for e, c in x_power.items()})

        c_N_inv = _poly_mul(c_N_inv, inv_factor)
        c_N_inv = {e: c for e, c in c_N_inv.items() if sum(e) <= 3}

    # c(TY) = c_amb * c_N_inv, termes de degre 3
    c_TY = _poly_mul(c_amb, c_N_inv)
    c3_terms = {e: c for e, c in c_TY.items() if sum(e) == 3}

    # Integrer c3 sur Y : integral = sum c3[e] * d_ijk...
    # Chaque monome J_r1^a1 ... J_rm^am avec sum ai = 3
    # correspond a l'integrale d_{indices} ou les indices sont
    # les Jr repetes selon leur exposant
    chi = 0
    for exp, coeff in c3_terms.items():
        # Convertir l'exposant en liste d'indices
        indices = []
        for r, e in enumerate(exp):
            indices.extend([r] * e)
        if len(indices) == 3:
            chi += coeff * d_ijk[indices[0], indices[1], indices[2]]

    return int(round(chi))


# ================================================================
# Seconde classe de Chern c2(TY)
# ================================================================

def compute_c2_tangent(
    ambient_dims: List[int],
    config_matrix: np.ndarray,
    d_ijk: np.ndarray,
) -> np.ndarray:
    """
    Calcule c2(TY) . J_i pour chaque direction de Kahler i.

    c2(TY) est dans H^4(Y). On retourne le vecteur
    c2_i = integral_Y c2(TY) wedge J_i = sum_{j,k} c2_{jk} d_{ijk}

    Formule: c2(TY) = c2(T_amb) - c2(N) puisque c1(TY) = 0.
    """
    m = len(ambient_dims)
    K = config_matrix.shape[0]

    # c2(T_amb)_{rs} = coeff de Jr Js dans prod (1+Ji)^{ni+1}
    #
    # CORRECTION : le terme de degre 2 de prod (1+J_r)^{n_r+1} vaut
    #     sum_r C(n_r+1, 2) J_r^2  +  sum_{r<s} (n_r+1)(n_s+1) J_r J_s
    # L'integration ci-dessous somme sur les couples ORDONNES (j,k) : le
    # coefficient hors diagonale doit donc etre la MOITIE de (n_r+1)(n_s+1),
    # sans quoi chaque terme croise est compte deux fois. c2_N, plus bas,
    # utilise deja la convention ordonnee -- les deux moities de la
    # soustraction doivent etre dans la meme convention.
    #
    # Verifications :
    #   bicubique P2xP2[3,3] : c2.J = (36, 36)   <- valeur connue
    #   quintique P4[5]      : c2.J = 50         <- inchange, m=1, aucun
    #                                               terme croise (ce qui
    #                                               explique que la
    #                                               validation initiale sur
    #                                               la quintique n'ait rien
    #                                               pu detecter)
    #   chi(L) entier sur 2800 tirages : 34 % avant, 100 % apres
    c2_amb = np.zeros((m, m))
    for r in range(m):
        c2_amb[r, r] = comb(ambient_dims[r] + 1, 2)
        for s in range(r + 1, m):
            val = (ambient_dims[r] + 1) * (ambient_dims[s] + 1) / 2.0
            c2_amb[r, s] = val
            c2_amb[s, r] = val

    # c2(N)_{rs} = sum_{a<b} q_{ar} q_{bs}  (partie symetrisee)
    c2_N = np.zeros((m, m))
    for a in range(K):
        for b in range(a + 1, K):
            for r in range(m):
                for s in range(m):
                    c2_N[r, s] += config_matrix[a, r] * config_matrix[b, s]

    c2_coeff = c2_amb - c2_N

    # Integrer: c2_i = sum_{j,k} c2_{jk} d_{ijk}
    c2_vec = np.zeros(m)
    for i in range(m):
        for j in range(m):
            for k in range(m):
                c2_vec[i] += c2_coeff[j, k] * d_ijk[i, j, k]

    return c2_vec


# ================================================================
# Hodge numbers (formule pour CICY favorables)
# ================================================================

def compute_hodge_favorable(
    ambient_dims: List[int],
    config_matrix: np.ndarray,
    d_ijk: np.ndarray,
) -> Tuple[int, int]:
    """
    Calcule (h11, h21) pour un CICY favorable (h11 = m).

    h11 = m (nombre de facteurs projectifs)
    h21 = h11 - chi/2 (depuis chi = 2(h11 - h21))
    """
    m = len(ambient_dims)
    h11 = m
    chi = compute_euler_from_intersection(ambient_dims, config_matrix, d_ijk)
    h21 = h11 - chi // 2
    return h11, h21


# ======================================================================
# ANNULATION D'ANOMALIE  --  condition (2.9) de arXiv:0911.1569
# ======================================================================
#
# Pour preserver la supersymetrie sans brane exotique, la classe duale a
# c2(TX) - c2(V) doit etre EFFECTIVE. Sur une CICY favorable, cela se lit
# composante par composante :
#
#     c2_r(V) <= c2_r(TX)   pour tout r
#
# Ce n'est pas un raffinement : c'est une condition de coherence de la
# theorie. Un fibre qui la viole n'est pas un modele, quelles que soient sa
# stabilite et sa cohomologie.
#
# Le pipeline ne la testait NULLE PART. Mesure sur le catalogue
# `scan_wilson2` : 70 entrees sur 115, soit 60,9 %, la violent. Les deux
# candidats du §2 la satisfont -- mais par chance, pas par construction.


def c2_somme_droites(d_ijk, charges):
    """
    c2 de (+) O(x_i), en composantes sur la base duale des J_r.

    c2( (+) O(x_i) ) = somme_{i<j} x_i . x_j, et

        somme_{i<j} x_i^s x_j^t = (1/2) [ (somme_i x_i^s)(somme_i x_i^t)
                                          - somme_i x_i^s x_i^t ]

    d'ou c2_r = (1/2) d_rst [ ... ]. Renvoie un vecteur de flottants : les
    c2 peuvent etre demi-entiers avant contraction, comme dans `c2_tangent`.
    """
    d = np.asarray(d_ijk, dtype=np.int64)
    m = d.shape[0]
    X = np.asarray(charges, dtype=np.int64).reshape(-1, m)
    tot = X.sum(axis=0)
    S = np.outer(tot, tot) - (X[:, :, None] * X[:, None, :]).sum(axis=0)
    return 0.5 * np.einsum('rst,st->r', d.astype(float), S.astype(float))


def c2_monade(d_ijk, b_charges, c_charges):
    """
    c2(V) pour 0 -> V -> B -> C -> 0.

    c(V) = c(B)/c(C) et c1(B) = c1(C) puisque c1(V) = 0, d'ou
    c2(V) = c2(B) - c2(C). En developpant, on retombe exactement sur la
    formule (2.9) de arXiv:0911.1569 :

        c2_r(V) = (1/2) d_rst [ somme_a c_a^s c_a^t - somme_i b_i^s b_i^t ]
    """
    return (c2_somme_droites(d_ijk, b_charges)
            - c2_somme_droites(d_ijk, c_charges))


def c2_extension(d_ijk, f1_charges, f2_charges):
    """
    c2(V) pour 0 -> F1 -> V -> F2 -> 0.

    La classe de Chern totale est multiplicative sur une suite exacte :
    c(V) = c(F1) c(F2), donc V a la meme classe que F1 (+) F2.
    """
    return c2_somme_droites(d_ijk, list(f1_charges) + list(f2_charges))


def anomalie_effective(c2_tangent, c2_V, tol=1e-9):
    """
    (ok, deficit) : c2(TX) - c2(V) est-elle effective ?

    `deficit` est le vecteur c2(TX) - c2(V) ; une composante negative
    identifie la direction fautive, ce qui vaut mieux qu'un booleen seul
    pour comprendre un rejet.
    """
    diff = [float(t) - float(v) for t, v in zip(c2_tangent, c2_V)]
    return all(x >= -tol for x in diff), diff
