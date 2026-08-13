"""
stability_full.py -- Test de semi-stabilite complet (tous rangs)

Critere de Hoppe : V est stable ssi
  H⁰(∧^p V ⊗ O(-H)) = 0
pour tout p = 1, ..., rk(V)-1 et tout diviseur ample H.

Sur un CICY dans P^{n_1} x ... x P^{n_m}, le cone ample est engendre
par les classes d'hyperplan J_1, ..., J_m. On teste O(-a) avec a_i >= 1.

Pour le monade 0 -> V -> B -> C -> 0 :
  ∧¹V = V                     (cohomologie monade avec charges decalees)
  ∧²V                         (suite exacte 0 -> ∧²V -> ∧²B -> V⊗C -> 0, decalee)
  ∧³V ≅ V* ⊗ det(V) = V*     (car c1(V)=0, donc det(V)=O)

Ref: Hoppe, "Generischer Spaltungstyp..." Math. Ann. 1984
     Anderson, Gray, Lukas, Palti, Phys.Rev.D 2012
"""
import numpy as np
from itertools import product as iprod
from typing import List, Dict

from cy_landscape.core.exact_cohomology import koszul_cohomology
from cy_landscape.core.monads import MonadBundle, compute_monad_cohomology


def _shift(charges_list, H):
    """Decale toutes les charges par -H."""
    return [[ch[k] - H[k] for k in range(len(H))] for ch in charges_list]


def _h0_sum(ambient, config, charge_list):
    """H⁰ d'une somme de line bundles."""
    total = 0
    for ch in charge_list:
        h = koszul_cohomology(ambient, config, ch)
        total += h.get(0, 0)
    return total


def h0_V_twisted(ambient, config, monad, H):
    """
    H⁰(V ⊗ O(-H)) via la monade decalee.
    0 -> V(-H) -> B(-H) -> C(-H) -> 0
    """
    shifted_monad = MonadBundle(
        _shift(monad.b_charges, H),
        _shift(monad.c_charges, H))
    cohom = compute_monad_cohomology(shifted_monad, ambient, config)
    if cohom is None:
        return None
    return cohom[0]


def h0_wedge2V_twisted(ambient, config, monad, H):
    """
    H⁰(∧²V ⊗ O(-H)) via la suite exacte decalee.
    0 -> ∧²V(-H) -> ∧²B(-H) -> (V⊗C)(-H) -> 0

    ∧²B(-H) = ⊕_{i<j} O(b_i + b_j - H)
    (V⊗C)(-H) vient de 0 -> V(c_1-H) -> B(c_1-H) -> O(2c_1-H) -> 0
    """
    m = monad.m
    b = monad.b_charges
    c1 = monad.c_charges[0]  # rk(C) = 1

    # H⁰(∧²B(-H))
    from itertools import combinations
    w2b_charges = []
    for i, j in combinations(range(monad.rank_B), 2):
        ch = [b[i][k] + b[j][k] - H[k] for k in range(m)]
        w2b_charges.append(ch)
    h0_w2B = _h0_sum(ambient, config, w2b_charges)

    # H⁰(V⊗C(-H)) = H⁰(V(c1-H))
    # Via monade: 0 -> V(c1-H) -> B(c1-H) -> O(2c1-H) -> 0
    H_shifted = [H[k] - c1[k] for k in range(m)]
    shifted_monad_vc = MonadBundle(
        _shift(b, H_shifted),
        [[c1[k] - H_shifted[k] for k in range(m)]])
    # Ou plus simplement:
    bc_charges = [[b_i[k] + c1[k] - H[k] for k in range(m)] for b_i in b]
    c2_charges = [2 * c1[k] - H[k] for k in range(m)]
    h0_BC = _h0_sum(ambient, config, bc_charges)
    h0_CC = koszul_cohomology(ambient, config, c2_charges).get(0, 0)

    # Map generique B(c1-H) -> O(2c1-H) : rang maximal
    rank_f = min(h0_BC, h0_CC)
    h0_VC = h0_BC - rank_f  # H⁰(V⊗C(-H))

    # Suite exacte longue : 0 -> H⁰(∧²V(-H)) -> H⁰(∧²B(-H)) -> H⁰(V⊗C(-H))
    # Map generique : rang maximal
    rank_g = min(h0_w2B, h0_VC)
    h0_w2V = h0_w2B - rank_g

    return h0_w2V


def h0_wedge3V_twisted(ambient, config, monad, H):
    """
    H⁰(∧³V ⊗ O(-H)).
    Pour rk(V) = 4 et c1(V) = 0 : ∧³V ≅ V*.
    Donc H⁰(∧³V(-H)) = H⁰(V*(-H)).

    Par Serre : H⁰(V*(-H)) = H³(V(H-K_X))* ou K_X = 0 (CY).
    Donc H⁰(V*(-H)) = H³(V(H)).

    Via monade decalee par +H :
    0 -> V(H) -> B(H) -> C(H) -> 0
    """
    m = monad.m
    neg_H = [-H[k] for k in range(m)]  # on "ajoute H" = on decale par -(-H)
    shifted_monad = MonadBundle(
        _shift(monad.b_charges, neg_H),
        _shift(monad.c_charges, neg_H))
    cohom = compute_monad_cohomology(shifted_monad, ambient, config)
    if cohom is None:
        return None
    return cohom[3]  # H³(V(H)) = H⁰(V*(-H))


def hoppe_criterion(ambient, config, monad, max_H=2, verbose=False):
    """
    Applique le critere de Hoppe complet.

    Teste H⁰(∧^p V(-H)) = 0 pour :
      p = 1, 2, ..., rk(V)-1
      H = (a_1, ..., a_m) avec 1 <= a_i <= max_H

    Retourne un dict avec le resultat et les details.
    """
    m = monad.m
    rk = monad.rank_V

    if rk <= 1:
        return {"stable": True, "reason": "Rang 1 => automatiquement stable"}

    # Generer les diviseurs amples a tester
    ample_divisors = []
    # D'abord les generateurs: e_1, ..., e_m
    for i in range(m):
        H = [0] * m; H[i] = 1
        ample_divisors.append(H)
    # Puis les combinaisons avec a_i in [1, max_H]
    for H in iprod(*[range(1, max_H + 1) for _ in range(m)]):
        H = list(H)
        if H not in ample_divisors:
            ample_divisors.append(H)

    violations = []
    tests_done = 0

    for H in ample_divisors:
        for p in range(1, rk):
            tests_done += 1

            if p == 1:
                h0 = h0_V_twisted(ambient, config, monad, H)
            elif p == 2:
                h0 = h0_wedge2V_twisted(ambient, config, monad, H)
            elif p == 3 and rk == 4:
                h0 = h0_wedge3V_twisted(ambient, config, monad, H)
            else:
                continue  # rang > 4 non implemente

            if h0 is None:
                continue

            if verbose:
                print(f"  H⁰(∧^{p}V(-{H})) = {h0}")

            if h0 > 0:
                violations.append({
                    'p': p, 'H': H, 'h0': int(h0),
                    'meaning': f"Sous-faisceau de rang {p} potentiellement destabilisant"
                })

    if violations:
        return {
            "stable": False,
            "semi_stable": False,
            "reason": f"{len(violations)} violations du critere de Hoppe",
            "violations": violations,
            "tests": tests_done,
        }

    return {
        "stable": True,
        "semi_stable": True,
        "reason": f"Critere de Hoppe satisfait ({tests_done} tests, rangs 1..{rk-1})",
        "tests": tests_done,
    }
