"""
bundles.py — Construction de fibrés vectoriels sur les variétés de Calabi-Yau.

En compactification hétérotique E₈ × E₈, un fibré vectoriel V sur le CY
brise un des facteurs E₈ selon :

    Structure group de V    →    Groupe de jauge résiduel (commutant dans E₈)
    ─────────────────────        ─────────────────────────────────────────────
    SU(2)                   →    E₇
    SU(3)                   →    E₆           ← GUT classique
    SU(4)                   →    SO(10)       ← GUT classique
    SU(5)                   →    SU(5)        ← GUT classique
    SU(3) × SU(2)           →    SU(3) × SU(2) × U(1)  ← Modèle Standard !

Le fibré doit satisfaire :
    1. Holomorphicité (automatique par construction)
    2. Stabilité (slope-stability au sens de Mumford-Takemoto)
    3. Annulation d'anomalie : c₂(V) + c₂(V') = c₂(TY)
       (simplification : c₂(V) = c₂(TY) si V' est trivial et pas de 5-branes)
    4. c₁(V) = 0 (condition SU(n) plutôt que U(n))

Types de fibrés implémentés :
    - Sommes de fibrés en droites (line bundle sums) : V = ⊕ L_i
    - Fibrés monades : 0 → V → B → C → 0  (exact sequence)
"""

from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict
import itertools
import numpy as np


# ═══════════════════════════════════════════════════════════════════
# Données géométriques d'un CY dans un produit d'espaces projectifs
# ═══════════════════════════════════════════════════════════════════

@dataclass
class CICYGeometry:
    """
    Géométrie complète d'un CICY défini dans un produit d'espaces projectifs.

    Un CICY dans P^{n1} × ... × P^{nm} est défini par K polynômes
    de multi-degrés donnés par la matrice de configuration.

    Attributes:
        ambient_dims : dimensions des espaces projectifs (n1, ..., nm)
        config_matrix : matrice de configuration K×m (degrés des polynômes)
        h11, h21 : nombres de Hodge
        intersection_numbers : nombres d'intersection triples d_{ijk}
                              (tenseur symétrique h11 × h11 × h11)
        c2_tangent : seconde classe de Chern du fibré tangent c₂(TY),
                     exprimée dans la base de H⁴ duale à {J_i ∧ J_j}
    """
    ambient_dims: List[int]            # ex: [4] pour P^4, [2, 2] pour P^2 × P^2
    config_matrix: np.ndarray          # matrice K × m
    h11: int
    h21: int
    intersection_numbers: np.ndarray   # tenseur d_{ijk}, shape (h11, h11, h11)
    c2_tangent: np.ndarray             # vecteur de longueur h11*(h11+1)/2

    @property
    def euler(self) -> int:
        return 2 * (self.h11 - self.h21)

    @property
    def n_gen(self) -> int:
        return abs(self.h11 - self.h21)

    @property
    def n_projective(self) -> int:
        """Nombre de facteurs d'espace projectif dans l'espace ambiant."""
        return len(self.ambient_dims)

    @property
    def n_kahler(self) -> int:
        """Nombre de paramètres de Kähler (= h11 pour les CICY favorables)."""
        return self.h11

    def triple_intersection(self, i: int, j: int, k: int) -> float:
        """Nombre d'intersection triple d_{ijk}."""
        return self.intersection_numbers[i, j, k]


# ═══════════════════════════════════════════════════════════════════
# Fibré vectoriel
# ═══════════════════════════════════════════════════════════════════

# Table des groupes de jauge résiduels (commutant du structure group dans E₈)
GAUGE_GROUP_TABLE = {
    1: {"group": "E₈",      "gut_viable": False, "sm_content": False},
    2: {"group": "E₇",      "gut_viable": False, "sm_content": False},
    3: {"group": "E₆",      "gut_viable": True,  "sm_content": True,
        "representations": {"27": "fondamentale", "27̄": "anti-fondamentale",
                            "1": "singlet"}},
    4: {"group": "SO(10)",   "gut_viable": True,  "sm_content": True,
        "representations": {"16": "spineur", "16̄": "anti-spineur",
                            "10": "vectorielle", "1": "singlet"}},
    5: {"group": "SU(5)",    "gut_viable": True,  "sm_content": True,
        "representations": {"10": "antisymétrique", "5̄": "anti-fondamentale",
                            "5": "fondamentale", "1": "singlet"}},
}


@dataclass
class LineBundleSum:
    """
    Somme de fibrés en droites : V = L₁ ⊕ L₂ ⊕ ... ⊕ Lₙ

    Chaque L_i est spécifié par un vecteur d'entiers (a_i1, ..., a_im)
    donnant les degrés dans chaque facteur de l'espace ambiant.

    La condition c₁(V) = 0 impose : Σ_i (a_i1, ..., a_im) = (0, ..., 0)
    (nécessaire pour que le structure group soit SU(n) et non U(n)).
    """
    charges: List[List[int]]   # Liste de n vecteurs de charges, chacun de dim h11

    @property
    def rank(self) -> int:
        return len(self.charges)

    @property
    def c1(self) -> np.ndarray:
        """Première classe de Chern c₁(V) = Σᵢ c₁(Lᵢ)."""
        return np.sum(self.charges, axis=0)

    @property
    def c1_vanishes(self) -> bool:
        """Vérifie c₁(V) = 0 (condition SU(n))."""
        return np.allclose(self.c1, 0)

    def c2(self, geometry: CICYGeometry) -> float:
        """
        Seconde classe de Chern c₂(V).

        Pour une somme de fibrés en droites :
            c₂(V) = -½ Σ_{i≠j} c₁(Lᵢ) · c₁(Lⱼ)
                   = ½ [c₁(V)² - Σᵢ c₁(Lᵢ)²]

        Le produit c₁(Lᵢ) · c₁(Lⱼ) utilise les nombres d'intersection.
        Pour un CY3, c₁·c₁ est un élément de H⁴(Y), qu'on évalue
        en contractant avec la forme de Kähler J :
            c₁(Lᵢ) · c₁(Lⱼ) = Σ_{a,b,c} aᵢₐ aⱼᵦ d_{abc} Jᶜ

        En pratique, on retourne le vecteur c₂ dans H⁴.
        Ici on calcule une version scalaire simplifiée.
        """
        n = len(self.charges)
        h = geometry.n_kahler
        d = geometry.intersection_numbers

        # c₂ = -½ Σ_{i<j} (aᵢ · aⱼ) via intersection numbers
        c2_val = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                # Produit c1(Li) ∧ c1(Lj) intégré sur le CY
                for a in range(h):
                    for b in range(h):
                        for c in range(h):
                            c2_val -= self.charges[i][a] * self.charges[j][b] * d[a, b, c]
        return c2_val

    def c3(self, geometry: CICYGeometry) -> float:
        """
        Troisième classe de Chern c₃(V).

        c₃(V) = ⅓ Σᵢ c₁(Lᵢ)³  (pour somme de line bundles avec c₁(V)=0)

        L'indice chiral (nombre net de générations) est :
            N_gen = ½|c₃(V)| = ½|Σᵢ d_{abc} aᵢₐ aᵢᵦ aᵢ꜀|
        """
        n = len(self.charges)
        h = geometry.n_kahler
        d = geometry.intersection_numbers

        c3_val = 0.0
        for i in range(n):
            for a in range(h):
                for b in range(h):
                    for c in range(h):
                        c3_val += (self.charges[i][a] *
                                   self.charges[i][b] *
                                   self.charges[i][c] * d[a, b, c])
        return c3_val / 3.0

    def chiral_index(self, geometry: CICYGeometry) -> float:
        """Indice chiral = ½|c₃(V)|. Doit être = 3 pour le Modèle Standard."""
        return abs(self.c3(geometry)) / 2.0

    def __repr__(self):
        charges_str = ", ".join(str(q) for q in self.charges)
        return f"LineBundleSum(rank={self.rank}, charges=[{charges_str}])"


@dataclass
class BundleAnalysis:
    """
    Résultat complet de l'analyse d'un fibré sur un CY.

    Regroupe toutes les informations physiques extraites :
    groupe de jauge, spectre de matière, validité des contraintes.
    """
    geometry: CICYGeometry
    bundle: LineBundleSum

    # Contraintes
    c1_vanishes: bool = False
    anomaly_cancelled: bool = False
    chiral_index_value: float = 0.0
    correct_generations: bool = False

    # Physique
    structure_group_rank: int = 0
    gauge_group: str = ""
    gut_viable: bool = False

    # Score composite
    sm_score: float = 0.0

    def compute(self) -> "BundleAnalysis":
        """Calcule toutes les propriétés physiques du fibré."""
        # 1. Condition SU(n)
        self.c1_vanishes = self.bundle.c1_vanishes

        # 2. Rang du structure group
        self.structure_group_rank = self.bundle.rank

        # 3. Groupe de jauge résiduel
        if self.structure_group_rank in GAUGE_GROUP_TABLE:
            info = GAUGE_GROUP_TABLE[self.structure_group_rank]
            self.gauge_group = info["group"]
            self.gut_viable = info["gut_viable"]
        else:
            self.gauge_group = f"Commutant(SU({self.structure_group_rank}), E₈)"
            self.gut_viable = False

        # 4. Indice chiral (nombre de générations du fibré)
        self.chiral_index_value = self.bundle.chiral_index(self.geometry)
        self.correct_generations = abs(self.chiral_index_value - 3.0) < 0.5

        # 5. Annulation d'anomalie (simplifiée)
        # En toute rigueur : c₂(V) = c₂(TY) - [W] où [W] est la classe des 5-branes
        # Ici on vérifie que c₂(V) ≤ c₂(TY) (condition d'effectivité)
        c2_bundle = self.bundle.c2(self.geometry)
        c2_tangent_scalar = float(np.sum(self.geometry.c2_tangent))
        self.anomaly_cancelled = (c2_bundle <= c2_tangent_scalar + 0.5)

        # 6. Score SM composite
        self._compute_score()

        return self

    def _compute_score(self):
        """
        Score de compatibilité avec le Modèle Standard (0–100).

        Critères pondérés :
          - c₁ = 0                    : 15 pts (condition nécessaire)
          - Groupe GUT viable         : 25 pts
          - 3 générations             : 30 pts (le plus discriminant)
          - Annulation d'anomalie     : 20 pts
          - SU(5) préféré (plus proche du SM) : 10 pts bonus
        """
        score = 0.0
        if self.c1_vanishes:
            score += 15.0
        if self.gut_viable:
            score += 25.0
        if self.correct_generations:
            score += 30.0
        elif self.chiral_index_value > 0:
            # Score partiel si on est proche de 3
            proximity = max(0, 1 - abs(self.chiral_index_value - 3) / 10)
            score += 30.0 * proximity
        if self.anomaly_cancelled:
            score += 20.0
        if self.gauge_group == "SU(5)":
            score += 10.0
        elif self.gauge_group == "SO(10)":
            score += 7.0
        elif self.gauge_group == "E₆":
            score += 5.0

        self.sm_score = min(100.0, score)

    def to_dict(self) -> dict:
        return {
            "h11": int(self.geometry.h11),
            "h21": int(self.geometry.h21),
            "euler": int(self.geometry.euler),
            "bundle_rank": int(self.structure_group_rank),
            "bundle_charges": [[int(x) for x in q] for q in self.bundle.charges],
            "c1_vanishes": bool(self.c1_vanishes),
            "gauge_group": str(self.gauge_group),
            "gut_viable": bool(self.gut_viable),
            "chiral_index": float(round(self.chiral_index_value, 2)),
            "correct_generations": bool(self.correct_generations),
            "anomaly_cancelled": bool(self.anomaly_cancelled),
            "sm_score": float(round(self.sm_score, 1)),
        }

    def summary_line(self) -> str:
        check = lambda b: "✓" if b else "✗"
        return (
            f"CY({self.geometry.h11},{self.geometry.h21}) | "
            f"SU({self.structure_group_rank})→{self.gauge_group:>6} | "
            f"c₁=0:{check(self.c1_vanishes)} | "
            f"N_gen={self.chiral_index_value:>5.1f} {check(self.correct_generations)} | "
            f"anom:{check(self.anomaly_cancelled)} | "
            f"score={self.sm_score:>5.1f}"
        )
