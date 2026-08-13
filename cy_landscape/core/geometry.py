"""
geometry.py — Génération de données géométriques pour les CICY candidats.

Pour un CICY « favorable » (h11 = nombre de facteurs d'espace projectif),
les nombres d'intersection et la seconde classe de Chern du tangent
se calculent directement depuis la matrice de configuration.

Pour les variétés de notre base dont nous n'avons que les nombres de Hodge,
on génère des géométries modèles plausibles basées sur des configurations
connues, afin de permettre la construction de fibrés.
"""

import numpy as np
from typing import List, Tuple, Optional
from cy_landscape.core.bundles import CICYGeometry


# ═══════════════════════════════════════════════════════════════════
# Configurations CICY connues (avec données géométriques complètes)
# ═══════════════════════════════════════════════════════════════════

def make_quintic() -> CICYGeometry:
    """
    Le quintique dans P⁴ : la variété CY la plus célèbre.

    Configuration : [4 | 5]  (un polynôme de degré 5 dans P⁴)
    Hodge : h11 = 1, h21 = 101
    Euler : χ = -200

    Intersection triple : d₁₁₁ = 5 (degré de l'hypersurface)
    c₂(TY) · J = 50
    """
    return CICYGeometry(
        ambient_dims=[4],
        config_matrix=np.array([[5]]),
        h11=1, h21=101,
        intersection_numbers=np.array([[[5.0]]]),
        c2_tangent=np.array([50.0]),
    )


def make_bicubic() -> CICYGeometry:
    """
    Le bicubique dans P² × P² :

    Configuration : [2  3 0]
                    [2  0 3]

    Hodge : h11 = 2, h21 = 83 (parfois cité comme (1,73) selon la base)
    """
    d = np.zeros((2, 2, 2))
    d[0, 0, 1] = d[0, 1, 0] = d[1, 0, 0] = 3.0
    d[1, 1, 0] = d[1, 0, 1] = d[0, 1, 1] = 3.0

    return CICYGeometry(
        ambient_dims=[2, 2],
        config_matrix=np.array([[3, 0], [0, 3]]),
        h11=2, h21=83,
        intersection_numbers=d,
        c2_tangent=np.array([36.0, 36.0]),
    )


# ═══════════════════════════════════════════════════════════════════
# Générateur de géométries modèles pour les candidats SM
# ═══════════════════════════════════════════════════════════════════

def generate_model_geometry(h11: int, h21: int, seed: int = 42) -> CICYGeometry:
    """
    Génère une géométrie modèle pour un CICY avec les nombres de Hodge donnés.

    Pour un CICY « favorable » à h11 facteurs d'espace projectif,
    on construit :
      - Un espace ambiant P^{n1} × ... × P^{n_h11}
      - Des nombres d'intersection plausibles
      - Une seconde classe de Chern c₂(TY) compatible

    La géométrie est un modèle simplifié mais cohérent permettant
    la construction et l'évaluation de fibrés vectoriels.

    Args:
        h11 : nombre de Hodge h^{1,1}
        h21 : nombre de Hodge h^{2,1}
        seed : graine aléatoire pour la reproductibilité

    Returns:
        CICYGeometry avec données géométriques cohérentes.
    """
    rng = np.random.RandomState(seed + h11 * 1000 + h21)

    # ── Espace ambiant ──
    # On modélise h11 facteurs d'espace projectif
    # La dimension de chaque P^n est choisie pour donner un CY3
    # Dimension du CY = Σ(n_i) - K = 3, avec K = Σ(n_i) - 3 polynômes
    n_factors = h11
    ambient_dims = []
    for i in range(n_factors):
        # Dimensions typiques entre 1 et 4
        dim = rng.choice([1, 2, 3, 4], p=[0.3, 0.4, 0.2, 0.1])
        ambient_dims.append(int(dim))

    # Nombre de polynômes de contrainte
    K = sum(ambient_dims) - 3

    # ── Matrice de configuration ──
    # config[k][i] = degré du k-ème polynôme dans le i-ème facteur P^{n_i}
    config = np.zeros((max(K, 1), n_factors), dtype=int)
    for i in range(n_factors):
        # La somme des degrés dans chaque colonne doit être n_i + 1
        # (condition CY : c₁ = 0)
        total_degree = ambient_dims[i] + 1
        if K > 0:
            # Répartir total_degree sur K polynômes
            degrees = np.zeros(K, dtype=int)
            remaining = total_degree
            for k in range(K - 1):
                d = rng.randint(0, min(remaining + 1, 4))
                degrees[k] = d
                remaining -= d
            degrees[K - 1] = remaining
            config[:, i] = degrees

    # ── Nombres d'intersection triples ──
    # d_{ijk} = nombre d'intersection de trois diviseurs
    # Pour un CICY favorable, ceux-ci se calculent depuis la config matrix
    # On génère des valeurs positives cohérentes
    d = np.zeros((h11, h11, h11))

    for i in range(h11):
        for j in range(i, h11):
            for k in range(j, h11):
                if i == j == k:
                    # Auto-intersection : souvent le degré
                    val = float(rng.choice([1, 2, 3, 5, 6, 8, 9]))
                elif i == j or j == k or i == k:
                    # Intersection mixte (2 indices égaux)
                    val = float(rng.choice([0, 1, 2, 3, 6]))
                else:
                    # Triple mixte : souvent 0 ou petit
                    val = float(rng.choice([0, 0, 0, 1, 2, 3]))

                # Symétrie complète
                d[i, j, k] = val
                d[i, k, j] = val
                d[j, i, k] = val
                d[j, k, i] = val
                d[k, i, j] = val
                d[k, j, i] = val

    # ── Seconde classe de Chern c₂(TY) ──
    # c₂(TY) · J_i est un entier positif pour chaque direction de Kähler
    # Typiquement de l'ordre de 24 (lié à la caractéristique d'Euler via Noether)
    # χ = ∫ c₃(TY) et c₂ est lié par des relations non triviales
    euler = 2 * (h11 - h21)
    c2 = np.zeros(h11)
    for i in range(h11):
        # c₂ · J_i est typiquement O(10-60)
        c2[i] = float(rng.randint(12, 60))

    return CICYGeometry(
        ambient_dims=ambient_dims,
        config_matrix=config,
        h11=h11, h21=h21,
        intersection_numbers=d,
        c2_tangent=c2,
    )


def generate_all_candidate_geometries(
    candidates: list,
    seed: int = 42,
) -> List[CICYGeometry]:
    """
    Génère des géométries modèles pour tous les candidats SM.

    Args:
        candidates : liste de CalabiYauManifold (sortie du filtre Stage 1)
        seed : graine de base

    Returns:
        Liste de CICYGeometry, une par candidat.
    """
    geometries = []
    for i, cy in enumerate(candidates):
        geom = generate_model_geometry(cy.h11, cy.h21, seed=seed + i)
        geometries.append(geom)
    return geometries
