"""
explicit_map.py -- Verification de l'approximation de map generique.

Pour une monade 0 -> V -> B -> C -> 0 avec rk(C) = 1 :
  f = (s_1, ..., s_{rB}) : B -> C
  s_i ∈ H⁰(O(c - b_i))  = polynomes multi-homogenes de degre (c - b_i)

La map induite sur H⁰ :
  f* : H⁰(B) = ⊕ H⁰(O(b_i)) -> H⁰(C) = H⁰(O(c))
  (σ_1, ..., σ_{rB}) |-> Σ s_i * σ_i

On verifie que rk(f*) = min(dim H⁰(B), dim H⁰(C)) (rang maximal).

Sur le CICY X ⊂ A = P^{n_1} x ... x P^{n_m} :
  H⁰(A, O(a)) = polynomes multi-homogenes de degre a = (a_1, ..., a_m)
  dim = ∏ C(a_i + n_i, n_i) si tous a_i >= 0, sinon 0

  H⁰(X, O(a)) est un quotient de H⁰(A, O(a)) par les relations de Koszul.
  Pour simplifier, on travaille sur l'ambiant et on projette.
"""
import numpy as np
from itertools import product as iprod
from math import comb
from typing import List, Tuple


def monomial_basis(n, d):
    """
    Base des monomes de degre d en (n+1) variables.
    Retourne une liste de tuples (e_0, ..., e_n) avec sum = d.
    """
    if n == 0:
        return [(d,)]
    basis = []
    for e0 in range(d + 1):
        for rest in monomial_basis(n - 1, d - e0):
            basis.append((e0,) + rest)
    return basis


def ambient_section_space(ambient_dims, degree):
    """
    Base de H⁰(A, O(degree)) = espace des polynomes multi-homogenes.
    Chaque element est un tuple de tuples : ((e_{1,0},...), (e_{2,0},...), ...)
    """
    if any(d < 0 for d in degree):
        return []

    factors = []
    for i, (n, d) in enumerate(zip(ambient_dims, degree)):
        factors.append(monomial_basis(n, d))

    # Produit tensoriel
    basis = []
    for combo in iprod(*factors):
        basis.append(combo)
    return basis


def ambient_dim(ambient_dims, degree):
    """dim H⁰(A, O(degree))."""
    if any(d < 0 for d in degree):
        return 0
    return int(np.prod([comb(d + n, n) for d, n in zip(degree, ambient_dims)]))


def multiply_sections(basis_s, basis_sigma, ambient_dims, target_degree):
    """
    Matrice de multiplication : s * sigma -> polynome de degre target.
    s ∈ H⁰(O(d_s)), sigma ∈ H⁰(O(d_sigma)), produit dans H⁰(O(target)).

    Retourne une matrice M telle que M[k, (i,j)] = coefficient du k-ieme
    monome target dans s_i * sigma_j.
    """
    target_basis = ambient_section_space(ambient_dims, target_degree)
    target_idx = {b: k for k, b in enumerate(target_basis)}

    n_target = len(target_basis)
    n_s = len(basis_s)
    n_sigma = len(basis_sigma)

    # Matrice : pour chaque (s, sigma), le produit est un monome cible
    # (car les bases sont des monomes, le produit de deux monomes est un monome)
    M = np.zeros((n_target, n_s * n_sigma), dtype=float)

    for i, s in enumerate(basis_s):
        for j, sigma in enumerate(basis_sigma):
            # Produit : additionner les exposants facteur par facteur
            prod = tuple(
                tuple(s[f][k] + sigma[f][k] for k in range(len(s[f])))
                for f in range(len(s))
            )
            if prod in target_idx:
                col = i * n_sigma + j
                M[target_idx[prod], col] = 1.0

    return M


def verify_generic_map(ambient_dims, config_matrix, b_charges, c_charges,
                        n_trials=5, rng=None):
    """
    Verifie l'approximation de map generique pour la monade.

    Construit explicitement la map f: H⁰(B) -> H⁰(C) sur l'ambiant,
    avec des sections s_i aleatoires, et verifie le rang.

    Retourne un dict avec les details.
    """
    if rng is None:
        rng = np.random.RandomState(42)

    m = len(ambient_dims)
    c = c_charges[0]  # rk(C) = 1

    # Dimensions sur l'ambiant
    dim_C_amb = ambient_dim(ambient_dims, c)
    dims_B_amb = [ambient_dim(ambient_dims, b) for b in b_charges]
    dim_B_amb = sum(dims_B_amb)

    if dim_C_amb == 0 or dim_B_amb == 0:
        return {'verified': True, 'reason': 'Espaces triviaux',
                'rank_expected': 0, 'rank_actual': 0}

    # Pour chaque b_i, l'espace des sections de O(c - b_i)
    map_degrees = [[c[k] - b[k] for k in range(m)] for b in b_charges]
    dims_map = [ambient_dim(ambient_dims, d) for d in map_degrees]

    if any(d == 0 for d in dims_map):
        return {'verified': False,
                'reason': 'Map partielle: certains H⁰(O(c-b_i)) = 0',
                'zero_entries': [i for i, d in enumerate(dims_map) if d == 0]}

    # Construire la matrice de la map f* : H⁰(B) -> H⁰(C)
    # f*(σ_1,...,σ_rB) = Σ s_i * σ_i
    # En coordonnees : matrice dim_C_amb x dim_B_amb

    rank_expected = min(dim_B_amb, dim_C_amb)
    ranks = []

    for trial in range(n_trials):
        # Matrice de f*
        F = np.zeros((dim_C_amb, dim_B_amb), dtype=float)

        col_offset = 0
        for i, b in enumerate(b_charges):
            deg_s = map_degrees[i]
            deg_b = list(b)

            basis_s = ambient_section_space(ambient_dims, deg_s)
            basis_b = ambient_section_space(ambient_dims, deg_b)

            if not basis_s or not basis_b:
                col_offset += dims_B_amb[i] if i < len(dims_B_amb) else 0
                continue

            # Section s_i aleatoire : combinaison lineaire des monomes
            s_coeffs = rng.randn(len(basis_s))

            # Matrice de multiplication s_i * . : H⁰(O(b_i)) -> H⁰(O(c))
            M = multiply_sections(basis_s, basis_b, ambient_dims, c)

            # Appliquer les coefficients de s_i
            # M_effective[k, j] = Σ_l s_coeffs[l] * M[k, l*n_b + j]
            n_b = len(basis_b)
            n_s = len(basis_s)
            M_eff = np.zeros((dim_C_amb, n_b), dtype=float)
            for l in range(n_s):
                M_eff += s_coeffs[l] * M[:, l * n_b:(l + 1) * n_b]

            F[:, col_offset:col_offset + n_b] = M_eff
            col_offset += n_b

        rank = np.linalg.matrix_rank(F, tol=1e-8)
        ranks.append(rank)

    avg_rank = np.mean(ranks)
    max_rank = max(ranks)
    min_rank = min(ranks)

    generic_ok = (max_rank == rank_expected)

    return {
        'verified': generic_ok,
        'rank_expected': rank_expected,
        'rank_actual_max': int(max_rank),
        'rank_actual_min': int(min_rank),
        'ranks': [int(r) for r in ranks],
        'dim_H0_B_amb': dim_B_amb,
        'dim_H0_C_amb': dim_C_amb,
        'n_trials': n_trials,
        'reason': ('Rang maximal confirme' if generic_ok
                   else f'Rang deficient: {max_rank} < {rank_expected}'),
    }
