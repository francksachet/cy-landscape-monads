"""
yukawa.py -- Couplages de Yukawa pour les fibres monades sur CICYs

Pour SO(10) avec V de rang 4, c1(V) = 0 :
  16  : H^1(V)       -> fermions chiraux (quarks, leptons)
  16bar: H^1(V*)      -> anti-fermions
  10  : H^1(∧²V)     -> Higgs

  Couplage 16_i · 16_j · 10_α :
    Y_{ijα} = ∫_X ω_i ∪ ω_j ∪ φ_α
  ou ω_i ∈ H^1(V), φ_α ∈ H^1(∧²V)

  Le cup product H^1(V) ⊗ H^1(V) ⊗ H^1(∧²V) -> H^3(∧^4 V) = H^3(O_X) ≅ C
  utilise l'application naturelle V ⊗ V ⊗ ∧²V -> ∧^4 V = O_X.

Pour un monad 0 -> V -> B -> C -> 0 :
  Les elements de H^1(V) vivent dans ker(H^1(B) -> H^1(C)).
  Chaque element a une decomposition en termes des H^1(O(b_i)).
  La TEXTURE du couplage (quels Y_{ijα} sont nuls) est determinee
  par les regles de selection de charges.

Ref: Anderson, Gray, Lukas, Palti "Yukawa Textures From Heterotic Stability Walls" (2011)
     Braun, He, Ovrut, Pantev "The Exact MSSM Spectrum from String Theory" (2006)
"""
import numpy as np
from itertools import product as iprod, combinations
from typing import List, Dict, Optional
from dataclasses import dataclass

from cy_landscape.core.exact_cohomology import koszul_cohomology
from cy_landscape.core.monads import MonadBundle


@dataclass
class YukawaResult:
    """Resultat du calcul des couplages de Yukawa."""
    gauge: str
    n_gen: int           # Nombre de generations
    n_higgs: int         # Nombre de champs Higgs
    texture: np.ndarray  # Matrice de texture Y_{ijα} (n_gen x n_gen x n_higgs)
    rank_per_higgs: List[int]  # Rang de Y_{..α} pour chaque Higgs α
    total_rank: int      # Rang maximal (sur tous les Higgs)
    eigenvalues: Optional[np.ndarray]  # Valeurs propres normalisees
    mass_hierarchy: Optional[float]    # Ratio top/bottom eigenvalue
    physical_summary: str


def _h1_decomposition(ambient, config, b_charges, c_charges):
    """
    Decompose H^1(V) en termes des contributions de chaque O(b_i).

    Pour 0 -> V -> B -> C -> 0 :
      H^0(B) -> H^0(C) -> H^1(V) -> H^1(B) -> H^1(C)

    H^1(V) = coker(H^0(B)->H^0(C)) ⊕ ker(H^1(B)->H^1(C))  (dans le cas generique)

    Chaque element de H^1(V) provenant de H^1(O(b_i)) porte la "charge" b_i.
    """
    m = len(ambient)
    h1_per_b = []  # (index, charges, h1_dim)
    for idx, b in enumerate(b_charges):
        h = koszul_cohomology(ambient, config, b)
        h1 = h.get(1, 0)
        if h1 > 0:
            h1_per_b.append((idx, b, int(h1)))

    # Pour C (rank 1)
    c = c_charges[0]
    h_c = koszul_cohomology(ambient, config, c)
    h1_c = h_c.get(1, 0)

    # H^1(V) = ker(H^1(B) -> H^1(C)) + coker(H^0(B) -> H^0(C))
    # La partie venant de H^1(B) porte les charges des b_i
    # La partie cokernel est plus subtile mais porte aussi des charges

    # Pour la TEXTURE, ce qui compte c'est quelles charges contribuent
    h0_per_b = []
    for idx, b in enumerate(b_charges):
        h = koszul_cohomology(ambient, config, b)
        h0_per_b.append((idx, b, int(h.get(0, 0))))

    h0_c = int(koszul_cohomology(ambient, config, c).get(0, 0))

    return h1_per_b, h0_per_b, h0_c, h1_c


def compute_yukawa_texture(ambient, config, monad, gauge="SO(10)"):
    """
    Calcule la texture des couplages de Yukawa.

    Pour SO(10) : 16_i · 16_j · 10_α
      Regle de selection : charges(ω_i) + charges(ω_j) + charges(φ_α)
      doit etre "compatible" (somme dans le cone effectif de C).

    On construit la matrice de texture T_{ijα} ∈ {0, 1} :
      T_{ijα} = 1 si le couplage est permis par les charges
      T_{ijα} = 0 sinon
    """
    m = monad.m
    b = monad.b_charges
    c = monad.c_charges[0]

    # Decomposer H^1(V)
    h1_per_b, h0_per_b, h0_c, h1_c = _h1_decomposition(
        ambient, config, b, monad.c_charges)

    # Les "charges" des elements de H^1(V)
    # Approximation : chaque element de H^1(O(b_i)) contribue avec charge b_i
    gen_charges = []  # charges des 16 (generatrices)
    for idx, bi, h1_dim in h1_per_b:
        for _ in range(h1_dim):
            gen_charges.append(list(bi))

    # Contributions du cokernel de f: H^0(C) -> H^1(V)
    # Le connecting morphism δ applique a σ ∈ H^0(O(c)) donne
    # un element de H^1(V) qui, dans la suite de Koszul, a une
    # composante de charge effective ~0 (la difference c - c = 0).
    # Plus precisement: l'obstruction vit dans la "direction C" du fibré,
    # ce qui revient a une charge nulle dans le complément V.
    total_h0_B = sum(h for _, _, h in h0_per_b)
    coker_dim = max(0, h0_c - total_h0_B)
    for _ in range(coker_dim):
        gen_charges.append([0] * m)  # Charge effective nulle

    # Charges des 10 (H^1(∧²V))
    # ∧²V vient de 0 -> ∧²V -> ∧²B -> V⊗C -> 0
    # Les elements de H^1(∧²B) portent les charges b_i + b_j
    higgs_charges = []
    for i, j in combinations(range(len(b)), 2):
        bi_bj = [b[i][k] + b[j][k] for k in range(m)]
        h = koszul_cohomology(ambient, config, bi_bj)
        h1 = h.get(1, 0)
        for _ in range(h1):
            higgs_charges.append(bi_bj)

    # Ajouter contributions de V⊗C
    for bi in b:
        bi_c = [bi[k] + c[k] for k in range(m)]
        h = koszul_cohomology(ambient, config, bi_c)
        h1 = h.get(1, 0)
        for _ in range(h1):
            higgs_charges.append(bi_c)

    n_gen = len(gen_charges)
    n_higgs = len(higgs_charges)

    if n_gen == 0 or n_higgs == 0:
        return YukawaResult(
            gauge=gauge, n_gen=0, n_higgs=0,
            texture=np.zeros((0, 0, 0)),
            rank_per_higgs=[], total_rank=0,
            eigenvalues=None, mass_hierarchy=None,
            physical_summary="Pas de couplage (spectre vide)")

    # Matrice de texture : T_{ijα} = 1 si charges compatibles
    # Regle : b_i + b_j + (b_k + b_l) = 2c  (conservation de charge totale)
    # Simplifie : somme des charges des 3 champs = 2c (pour la trace sur ∧^4 V)
    texture = np.zeros((n_gen, n_gen, n_higgs), dtype=float)
    target = [0] * m  # charge totale = 0 (car det V = O_X, trace ∧^4 V -> O)

    for i in range(n_gen):
        for j in range(i, n_gen):  # Symetrique en i,j
            for alpha in range(n_higgs):
                total = [gen_charges[i][k] + gen_charges[j][k] + higgs_charges[alpha][k]
                         for k in range(m)]
                if total == target:
                    texture[i, j, alpha] = 1.0
                    texture[j, i, alpha] = 1.0  # Symetrie

    # Rang par Higgs
    rank_per_higgs = []
    for alpha in range(n_higgs):
        Y_alpha = texture[:, :, alpha]
        rank_per_higgs.append(int(np.linalg.matrix_rank(Y_alpha)))

    total_rank = max(rank_per_higgs) if rank_per_higgs else 0

    # Valeurs propres (de la plus grande matrice Y_α)
    eigenvalues = None
    mass_hierarchy = None
    if total_rank > 0:
        best_alpha = np.argmax(rank_per_higgs)
        Y_best = texture[:, :, best_alpha]
        # Y est symetrique -> valeurs propres reelles
        evals = np.sort(np.abs(np.linalg.eigvalsh(Y_best)))[::-1]
        evals = evals[evals > 1e-10]
        if len(evals) > 0:
            eigenvalues = evals / evals[0]  # Normaliser
            if len(evals) > 1 and evals[-1] > 0:
                mass_hierarchy = float(evals[0] / evals[-1])

    # Resume physique
    n_massive = total_rank
    n_massless = n_gen - n_massive
    summary_parts = [
        f"Spectre {gauge} : {n_gen} generations",
        f"Couplages de Yukawa : matrice {n_gen}x{n_gen} (x{n_higgs} Higgs)",
        f"Rang maximal : {total_rank}/{n_gen}",
        f"  -> {n_massive} generation(s) massive(s)",
        f"  -> {n_massless} generation(s) sans masse au niveau tree",
    ]
    if mass_hierarchy is not None and mass_hierarchy > 1:
        summary_parts.append(f"Hierarchie de masse : {mass_hierarchy:.1f}:1")
        if mass_hierarchy > 100:
            summary_parts.append("  -> Compatible avec t/u ~ 40000:1 du SM")
        elif mass_hierarchy > 10:
            summary_parts.append("  -> Hierarchie moderee")
        else:
            summary_parts.append("  -> Hierarchie faible (quasi-degenere)")
    if eigenvalues is not None and len(eigenvalues) >= 3:
        summary_parts.append(
            f"Masses normalisees : {eigenvalues[0]:.3f}, {eigenvalues[1]:.3f}, {eigenvalues[2]:.3f}")

    summary = "\n".join(summary_parts)

    return YukawaResult(
        gauge=gauge, n_gen=n_gen, n_higgs=n_higgs,
        texture=texture, rank_per_higgs=rank_per_higgs,
        total_rank=total_rank, eigenvalues=eigenvalues,
        mass_hierarchy=mass_hierarchy, physical_summary=summary)
