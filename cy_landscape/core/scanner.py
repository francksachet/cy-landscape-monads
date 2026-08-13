"""
scanner.py — Recherche systématique de fibrés vectoriels viables.

Pour chaque géométrie CY candidate, on explore un espace de fibrés
(sommes de fibrés en droites) en cherchant ceux qui satisfont :
  1. c₁(V) = 0  (structure group SU(n))
  2. Indice chiral = 3  (trois générations)
  3. Annulation d'anomalie  (c₂(V) ≤ c₂(TY))
  4. Groupe GUT viable  (rang 3, 4 ou 5 → E₆, SO(10), SU(5))

La recherche est paramétrique : on balaye les charges entières
dans un intervalle [-charge_max, charge_max] pour chaque fibré en droite.
"""

import itertools
import numpy as np
from typing import List, Tuple, Optional, Dict
from dataclasses import dataclass

from cy_landscape.core.bundles import (
    CICYGeometry, LineBundleSum, BundleAnalysis, GAUGE_GROUP_TABLE,
)


@dataclass
class ScanConfig:
    """Configuration du scan de fibrés."""
    bundle_ranks: List[int] = None       # Rangs à explorer (3, 4, 5)
    charge_max: int = 3                  # Charge maximale par entrée
    max_bundles_per_rank: int = 500      # Limite par rang pour chaque CY
    min_score: float = 50.0              # Score minimum pour retenir un résultat
    require_c1_zero: bool = True         # Exiger c₁ = 0
    target_generations: int = 3          # Nombre de générations cible

    def __post_init__(self):
        if self.bundle_ranks is None:
            self.bundle_ranks = [3, 4, 5]


class BundleScanner:
    """
    Moteur de recherche de fibrés vectoriels viables.

    Pour chaque rang r ∈ {3, 4, 5} et chaque géométrie CY :
      1. Génère des sommes de r fibrés en droites avec c₁ = 0
      2. Calcule l'indice chiral et vérifie les contraintes
      3. Classe les résultats par score SM
    """

    def __init__(self, config: ScanConfig = None):
        self.config = config or ScanConfig()
        self.results: List[BundleAnalysis] = []
        self._scan_stats = {
            "geometries_scanned": 0,
            "bundles_tested": 0,
            "bundles_valid_c1": 0,
            "bundles_3gen": 0,
            "bundles_anomaly_ok": 0,
            "bundles_retained": 0,
        }

    def scan_geometry(self, geometry: CICYGeometry) -> List[BundleAnalysis]:
        """
        Scanne une géométrie CY pour trouver des fibrés viables.

        Returns:
            Liste de BundleAnalysis triés par score décroissant.
        """
        self._scan_stats["geometries_scanned"] += 1
        results = []

        for rank in self.config.bundle_ranks:
            rank_results = self._scan_rank(geometry, rank)
            results.extend(rank_results)

        # Trier par score
        results.sort(key=lambda r: r.sm_score, reverse=True)
        self.results.extend(results)
        return results

    def scan_all(self, geometries: List[CICYGeometry]) -> List[BundleAnalysis]:
        """Scanne toutes les géométries et retourne les résultats globaux."""
        all_results = []
        for geom in geometries:
            results = self.scan_geometry(geom)
            all_results.extend(results)

        all_results.sort(key=lambda r: r.sm_score, reverse=True)
        return all_results

    def _scan_rank(
        self, geometry: CICYGeometry, rank: int,
    ) -> List[BundleAnalysis]:
        """
        Recherche de sommes de `rank` fibrés en droites sur `geometry`.

        Stratégie : on génère les charges des (rank-1) premiers fibrés,
        puis on déduit le dernier par la condition c₁ = 0.
        Cela garantit automatiquement c₁(V) = 0.
        """
        h = geometry.n_kahler
        cmax = self.config.charge_max
        results = []
        count = 0

        # Bornes de recherche
        charge_range = range(-cmax, cmax + 1)

        # Pour réduire l'espace de recherche, on utilise des stratégies
        # adaptées à chaque rang
        if rank <= 3:
            bundles = self._generate_bundles_small_rank(h, rank, cmax)
        else:
            bundles = self._generate_bundles_sampling(h, rank, cmax)

        for charges in bundles:
            if count >= self.config.max_bundles_per_rank:
                break

            bundle = LineBundleSum(charges=charges)
            self._scan_stats["bundles_tested"] += 1

            if not bundle.c1_vanishes:
                continue
            self._scan_stats["bundles_valid_c1"] += 1

            # Analyser le fibré
            analysis = BundleAnalysis(
                geometry=geometry,
                bundle=bundle,
            ).compute()

            if analysis.correct_generations:
                self._scan_stats["bundles_3gen"] += 1
            if analysis.anomaly_cancelled:
                self._scan_stats["bundles_anomaly_ok"] += 1

            if analysis.sm_score >= self.config.min_score:
                results.append(analysis)
                self._scan_stats["bundles_retained"] += 1
                count += 1

        return results

    def _generate_bundles_small_rank(
        self, h: int, rank: int, cmax: int,
    ) -> List[List[List[int]]]:
        """
        Génération exhaustive (modérée) pour rangs 3-5.

        On fixe (rank-1) vecteurs de charges et on déduit le dernier
        par c₁ = 0. Pour réduire l'espace, on limite h ≤ 4 et
        on utilise des charges dans [-cmax, cmax].
        """
        h_eff = min(h, 4)  # Limiter pour la tractabilité
        charge_range = list(range(-cmax, cmax + 1))
        bundles = []

        # Générer les (rank-1) premiers fibrés
        if rank == 3 and h_eff <= 2:
            # Exploration quasi-exhaustive pour le cas le plus petit
            for combo in itertools.product(
                itertools.product(charge_range, repeat=h_eff),
                repeat=rank - 1,
            ):
                # Dernier fibré déduit par c₁ = 0
                partial = list(combo)
                last = [-sum(q[i] for q in partial) for i in range(h_eff)]

                # Vérifier que les charges du dernier sont dans les bornes
                if all(-cmax <= l <= cmax for l in last):
                    charges = [list(q) for q in partial] + [last]
                    # Compléter à la dimension h si nécessaire
                    if h > h_eff:
                        charges = [q + [0] * (h - h_eff) for q in charges]
                    bundles.append(charges)

                if len(bundles) >= 2000:
                    break
        else:
            # Échantillonnage pour les cas plus grands
            bundles = self._generate_bundles_sampling(h, rank, cmax)

        return bundles

    def _generate_bundles_sampling(
        self, h: int, rank: int, cmax: int, n_samples: int = 1000,
    ) -> List[List[List[int]]]:
        """
        Génération par échantillonnage aléatoire pour les grands espaces.

        Stratégie : échantillonner (rank-1) vecteurs, déduire le dernier.
        Inclut aussi des configurations structurées connues pour être
        prometteuses (charges alternées, motifs cyclopermutés, etc.).
        """
        rng = np.random.RandomState(42)
        h_eff = min(h, 6)
        bundles = []

        # 1. Configurations structurées
        structured = self._structured_bundles(h_eff, rank, cmax)
        for s in structured:
            if h > h_eff:
                s = [q + [0] * (h - h_eff) for q in s]
            bundles.append(s)

        # 2. Échantillonnage aléatoire
        for _ in range(n_samples):
            partial = []
            for _ in range(rank - 1):
                q = [int(rng.randint(-cmax, cmax + 1)) for _ in range(h_eff)]
                partial.append(q)
            last = [-sum(q[i] for q in partial) for i in range(h_eff)]
            if all(-cmax - 1 <= l <= cmax + 1 for l in last):
                charges = partial + [last]
                if h > h_eff:
                    charges = [q + [0] * (h - h_eff) for q in charges]
                bundles.append(charges)

        return bundles

    def _structured_bundles(
        self, h: int, rank: int, cmax: int,
    ) -> List[List[List[int]]]:
        """
        Fibrés structurés connus pour être physiquement intéressants.

        Inclut :
        - Charges unitaires alternées (+1, -1, 0, ...)
        - Motifs cycliques
        - Charges "standard embedding" (V = TY)
        """
        bundles = []

        # Pattern 1 : charges unitaires avec alternance de signes
        if rank <= h:
            # Ex pour rank=3, h=6 : [[1,0,0,-1,0,0], [0,1,0,0,-1,0], [0,0,1,0,0,-1]]
            for shift in range(min(h - rank + 1, 3)):
                charges = []
                for r in range(rank - 1):
                    q = [0] * h
                    idx = (r + shift) % h
                    q[idx] = 1
                    charges.append(q)
                last = [-sum(q[i] for q in charges) for i in range(h)]
                charges.append(last)
                bundles.append(charges)

        # Pattern 2 : charges (1, -1) dans des directions différentes
        if h >= 2:
            for i in range(min(h, 3)):
                for j in range(i + 1, min(h, 4)):
                    charges = []
                    for r in range(rank - 1):
                        q = [0] * h
                        if r % 2 == 0:
                            q[i] = 1
                            q[j] = -1
                        else:
                            q[i] = -1
                            q[j] = 1
                        charges.append(q)
                    last = [-sum(q[k] for q in charges) for k in range(h)]
                    charges.append(last)
                    bundles.append(charges)

        # Pattern 3 : charges croissantes
        for scale in [1, 2]:
            charges = []
            for r in range(rank - 1):
                q = [0] * h
                q[r % h] = scale
                q[(r + 1) % h] = -scale
                charges.append(q)
            last = [-sum(q[i] for q in charges) for i in range(h)]
            charges.append(last)
            bundles.append(charges)

        return bundles

    @property
    def stats(self) -> dict:
        return dict(self._scan_stats)

    def top_results(self, n: int = 20) -> List[BundleAnalysis]:
        """Les n meilleurs résultats par score SM."""
        return sorted(self.results, key=lambda r: r.sm_score, reverse=True)[:n]

    def results_by_gauge_group(self) -> Dict[str, List[BundleAnalysis]]:
        """Résultats groupés par groupe de jauge."""
        groups = {}
        for r in self.results:
            g = r.gauge_group
            if g not in groups:
                groups[g] = []
            groups[g].append(r)
        return groups
