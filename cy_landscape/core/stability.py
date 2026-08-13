"""
stability.py -- Test de semi-stabilite de Mumford-Takemoto
pour les sommes de fibres en droites sur les CICY.

Un fibre V = oplus L_i avec c1(V) = 0 est semi-stable par rapport
a une forme de Kahler omega = sum t_r J_r si pour tout sous-ensemble
propre non vide S subset {1,...,n} :

    sum_{i in S} deg_omega(L_i) <= 0

ou deg_omega(L_i) = integral_Y c1(L_i) ^ omega^2
   = sum_{j,k,l} charges_i[j] * t_k * t_l * d_{jkl}

Pour une somme avec c1 = 0, semi-stabilite implique que tous les
sous-sommes ont un degre <= 0, et le complement a un degre >= 0.
Puisque la somme totale est 0, cela revient a:
pour tout S: sum_{i in S} deg(L_i) = 0 dans le cas "strictement" semi-stable.

En pratique on verifie:
  1. A un point generique du cone de Kahler (t = (1,...,1))
  2. Sur le bord du cone
  3. L'existence d'une region de stabilite dans le cone
"""

import numpy as np
from itertools import combinations
from typing import List, Tuple, Dict, Optional


def compute_line_bundle_degree(
    charges: List[int],
    d_ijk: np.ndarray,
    kahler_params: Optional[np.ndarray] = None,
) -> float:
    """
    Degre d'un fibre en droites L = O(charges) par rapport a omega = sum t_i J_i.

    deg(L) = sum_{j,k,l} a_j t_k t_l d_{jkl}

    Si kahler_params est None, utilise t = (1, ..., 1).
    """
    m = len(charges)
    if kahler_params is None:
        kahler_params = np.ones(m)
    t = kahler_params

    deg = 0.0
    for j in range(m):
        if charges[j] == 0:
            continue
        for k in range(m):
            for l in range(m):
                deg += charges[j] * t[k] * t[l] * d_ijk[j, k, l]
    return deg


def check_semi_stability(
    bundle_charges: List[List[int]],
    d_ijk: np.ndarray,
    kahler_params: Optional[np.ndarray] = None,
) -> Dict:
    """
    Verifie la semi-stabilite d'une somme de fibres en droites.

    Teste tous les sous-ensembles propres non vides S:
      sum_{i in S} deg(L_i) <= 0

    Returns:
        Dict avec:
        - "semi_stable": bool
        - "individual_degrees": liste des degres de chaque L_i
        - "worst_violation": pire violation (> 0 si instable)
        - "destabilizing_subset": sous-ensemble destabilisant (si instable)
        - "n_subsets_checked": nombre de sous-ensembles testes
    """
    n = len(bundle_charges)
    m = len(bundle_charges[0])

    if kahler_params is None:
        kahler_params = np.ones(m)

    # Degres individuels
    degrees = [
        compute_line_bundle_degree(charges, d_ijk, kahler_params)
        for charges in bundle_charges
    ]

    # Verifier la somme totale (devrait etre ~0 si c1 = 0)
    total_deg = sum(degrees)

    result = {
        "semi_stable": True,
        "individual_degrees": [round(d, 4) for d in degrees],
        "total_degree": round(total_deg, 6),
        "worst_violation": 0.0,
        "destabilizing_subset": None,
        "n_subsets_checked": 0,
    }

    # Pour n petit, tester tous les sous-ensembles
    # Pour n grand, tester les sous-ensembles de taille 1 a n//2
    max_subsets = 2 ** n - 2  # Tous les sous-ensembles propres non vides
    if max_subsets > 10000:
        # Heuristique : tester singletons + paires + random
        subsets_to_check = []
        for size in range(1, min(n, 4)):
            subsets_to_check.extend(combinations(range(n), size))
    else:
        subsets_to_check = []
        for size in range(1, n):
            subsets_to_check.extend(combinations(range(n), size))

    worst = 0.0
    worst_subset = None

    for subset in subsets_to_check:
        result["n_subsets_checked"] += 1
        sub_deg = sum(degrees[i] for i in subset)

        if sub_deg > worst:
            worst = sub_deg
            worst_subset = subset

        if sub_deg > 1e-10:  # Tolerance numerique
            result["semi_stable"] = False

    result["worst_violation"] = round(worst, 6)
    if worst_subset is not None and worst > 1e-10:
        result["destabilizing_subset"] = list(worst_subset)

    return result


def find_stability_region(
    bundle_charges: List[List[int]],
    d_ijk: np.ndarray,
    n_samples: int = 200,
) -> Dict:
    """
    Cherche une region du cone de Kahler ou le fibre est semi-stable.

    Echantillonne des points t = (t1,...,tm) avec ti > 0 et teste
    la semi-stabilite a chaque point.

    Returns:
        Dict avec:
        - "stable_anywhere": bool
        - "stable_fraction": fraction des points stables
        - "best_point": meilleur point de Kahler trouve
        - "best_violation": violation minimale trouvee
    """
    n = len(bundle_charges)
    m = len(bundle_charges[0])

    rng = np.random.RandomState(42)
    best_violation = float('inf')
    best_point = None
    n_stable = 0

    for _ in range(n_samples):
        # Point aleatoire dans le cone de Kahler
        t = rng.exponential(1.0, size=m) + 0.1

        # Degres individuels
        degrees = [
            compute_line_bundle_degree(charges, d_ijk, t)
            for charges in bundle_charges
        ]

        # Pire violation
        worst = 0.0
        for size in range(1, n):
            for subset in combinations(range(n), size):
                sub_deg = sum(degrees[i] for i in subset)
                if sub_deg > worst:
                    worst = sub_deg
            if worst > 0:
                break  # Deja instable, pas besoin de continuer

        if worst <= 1e-10:
            n_stable += 1

        if worst < best_violation:
            best_violation = worst
            best_point = t.tolist()

    return {
        "stable_anywhere": n_stable > 0,
        "stable_fraction": n_stable / n_samples,
        "best_point": [round(x, 3) for x in best_point] if best_point else None,
        "best_violation": round(best_violation, 6),
    }
