"""
cicy_matrices.py -- Matrices de configuration CICY reelles.

Chaque entree contient :
  - "name": identifiant
  - "ambient": [n1, ..., nm] dimensions des P^ni
  - "config": matrice K x m (numpy array)
  - "h_known": (h11, h21) si connus (pour verification)

Sources : Candelas, Dale, Lutken, Schimmrigk (1988);
          Green, Hubsch, Lutken (1989); pyCICY database.
"""

import numpy as np


CICY_CONFIGS = [
    # ================================================================
    # CICY a 1 facteur projectif
    # ================================================================
    {
        "name": "quintic",
        "ambient": [4],
        "config": np.array([[5]]),
        "h_known": (1, 101),
        "description": "Le quintique dans P^4 - le CY le plus celebre",
    },

    # ================================================================
    # CICY a 2 facteurs projectifs
    # ================================================================
    {
        "name": "P1xP4_41_11",
        "ambient": [1, 4],
        "config": np.array([[1, 4], [1, 1]]),
        "h_known": None,
        "description": "Intersection complete dans P^1 x P^4",
    },
    {
        "name": "P2xP3_12_22",
        "ambient": [2, 3],
        "config": np.array([[1, 2], [2, 2]]),
        "h_known": None,
        "description": "Intersection complete dans P^2 x P^3",
    },
    {
        "name": "P1xP5_24_12_10",
        "ambient": [1, 5],
        "config": np.array([[2, 4], [0, 1], [0, 1]]),
        "h_known": None,
        "description": "Codimension 3 dans P^1 x P^5",
    },

    # ================================================================
    # CICY a 3 facteurs projectifs
    # ================================================================
    {
        "name": "P1xP1xP3_11_11_12",
        "ambient": [1, 1, 3],
        "config": np.array([[1, 1, 1], [1, 1, 3]]),
        "h_known": None,
        "description": "Dans P^1 x P^1 x P^3",
    },
    {
        "name": "P2xP2xP1_111_111",
        "ambient": [2, 2, 1],
        "config": np.array([[1, 1, 1], [1, 1, 1], [1, 1, 0]]),
        "h_known": None,
        "description": "Dans P^2 x P^2 x P^1",
    },
    {
        "name": "P1xP2xP3_11_12_11",
        "ambient": [1, 2, 3],
        "config": np.array([[1, 1, 2], [1, 2, 2]]),
        "h_known": None,
        "description": "Dans P^1 x P^2 x P^3",
    },

    # ================================================================
    # CICY a 4 facteurs projectifs
    # ================================================================
    {
        "name": "P1x4_tetraquadric",
        "ambient": [1, 1, 1, 1],
        "config": np.array([[2, 0, 0, 0], [0, 2, 0, 0],
                            [0, 0, 2, 0], [0, 0, 0, 2]]),
        "h_known": None,
        "description": "Tetra-quadrique dans (P^1)^4 - CY a 4 facteurs",
    },
    {
        "name": "P1x3xP2_111_11_11",
        "ambient": [1, 1, 1, 2],
        "config": np.array([[1, 1, 1, 1], [1, 1, 1, 2]]),
        "h_known": None,
        "description": "Dans (P^1)^3 x P^2",
    },
    {
        "name": "P1x2xP2x2",
        "ambient": [1, 1, 2, 2],
        "config": np.array([[1, 0, 1, 1], [0, 1, 1, 1],
                            [1, 1, 1, 1]]),
        "h_known": None,
        "description": "Dans (P^1)^2 x (P^2)^2",
    },

    # ================================================================
    # CICY a 5 facteurs projectifs (potentiellement 3 generations)
    # ================================================================
    {
        "name": "P1x5_penta",
        "ambient": [1, 1, 1, 1, 1],
        "config": np.array([[2, 0, 0, 0, 0], [0, 2, 0, 0, 0],
                            [0, 0, 2, 0, 0], [0, 0, 0, 2, 0],
                            [0, 0, 0, 0, 2]]),
        "h_known": None,
        "description": "5 quadriques dans (P^1)^5",
    },
    {
        "name": "P1x5_mixed_A",
        "ambient": [1, 1, 1, 1, 1],
        "config": np.array([[1, 1, 0, 0, 0], [0, 1, 1, 0, 0],
                            [0, 0, 1, 1, 0], [0, 0, 0, 1, 1],
                            [1, 0, 0, 0, 1]]),
        "h_known": None,
        "description": "Configuration cyclique dans (P^1)^5",
    },
    {
        "name": "P1x3xP2xP1",
        "ambient": [1, 1, 1, 2, 1],
        "config": np.array([[1, 1, 0, 1, 0], [0, 0, 1, 1, 1],
                            [1, 1, 1, 1, 1]]),
        "h_known": None,
        "description": "Config mixte a 5 facteurs",
    },

    # ================================================================
    # CICY a 6+ facteurs (candidats SM si |h11-h21|=3)
    # ================================================================
    {
        "name": "P1x6_A",
        "ambient": [1, 1, 1, 1, 1, 1],
        "config": np.array([[2, 0, 0, 0, 0, 0], [0, 2, 0, 0, 0, 0],
                            [0, 0, 2, 0, 0, 0], [0, 0, 0, 2, 0, 0],
                            [0, 0, 0, 0, 2, 0], [0, 0, 0, 0, 0, 2]]),
        "h_known": None,
        "description": "6 quadriques diagonales dans (P^1)^6",
    },
    {
        "name": "P1x6_B",
        "ambient": [1, 1, 1, 1, 1, 1],
        "config": np.array([[1, 1, 0, 0, 0, 0], [0, 1, 1, 0, 0, 0],
                            [0, 0, 1, 1, 0, 0], [0, 0, 0, 1, 1, 0],
                            [0, 0, 0, 0, 1, 1], [1, 0, 0, 0, 0, 1]]),
        "h_known": None,
        "description": "Configuration cyclique dans (P^1)^6",
    },
    {
        "name": "P1x4xP2x2_A",
        "ambient": [1, 1, 1, 1, 2, 2],
        "config": np.array([[1, 0, 0, 0, 1, 1], [0, 1, 0, 0, 1, 1],
                            [0, 0, 1, 0, 1, 0], [0, 0, 0, 1, 0, 1],
                            [1, 1, 1, 1, 0, 0]]),
        "h_known": None,
        "description": "Config mixte a 6 facteurs",
    },
]


def get_config_by_name(name: str) -> dict:
    """Retourne une configuration par son nom."""
    for c in CICY_CONFIGS:
        if c["name"] == name:
            return c
    raise ValueError(f"Configuration '{name}' non trouvee")


def get_all_configs() -> list:
    """Retourne toutes les configurations."""
    return CICY_CONFIGS


def validate_config(cfg: dict) -> dict:
    """
    Verifie la validite d'une configuration CICY :
    - dim(CY) = sum(ni) - K = 3
    - Condition CY : sum_a q_{a,i} = n_i + 1 pour chaque i
    """
    ambient = cfg["ambient"]
    config = cfg["config"]
    m = len(ambient)
    K = config.shape[0]

    # Dimension
    dim_cy = sum(ambient) - K
    dim_ok = (dim_cy == 3)

    # CY condition
    col_sums = config.sum(axis=0)
    expected = np.array([n + 1 for n in ambient])
    cy_ok = np.allclose(col_sums, expected)

    return {
        "name": cfg["name"],
        "m": m,
        "K": K,
        "dim_CY": dim_cy,
        "dim_ok": dim_ok,
        "cy_condition": cy_ok,
        "col_sums": col_sums.tolist(),
        "expected": expected.tolist(),
        "valid": dim_ok and cy_ok,
    }
