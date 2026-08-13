"""
database.py — Classe centrale de la base de données CY.

Charge les données CICY, calcule les invariants topologiques dérivés,
et fournit une interface de requête et de filtrage.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple
import json


@dataclass
class CalabiYauManifold:
    """
    Représentation d'une variété de Calabi-Yau par ses invariants topologiques.

    Attributs principaux :
        h11 : nombre de Hodge h^{1,1} — compte les déformations de Kähler
              (intuition : "taille et forme" des cycles internes)
        h21 : nombre de Hodge h^{2,1} — compte les déformations de structure complexe
              (intuition : "forme fine" de la géométrie)
        count : nombre de matrices de configuration CICY distinctes donnant cette topologie
        source : base de données d'origine ("CICY", "KS", "gCICY")

    Attributs dérivés (calculés automatiquement) :
        euler : caractéristique d'Euler χ = 2(h11 - h21)
        n_gen : nombre de générations de fermions |χ|/2 (compactification hétérotique)
        hodge_total : nombre total de cycles h11 + h21 + 2
                      (le +2 vient de h00 = h33 = 1 pour un CY3)
        is_mirror_self : True si h11 == h21 (variété auto-miroir)
    """
    h11: int
    h21: int
    count: int = 1
    source: str = "CICY"

    # Champs dérivés, calculés dans __post_init__
    euler: int = field(init=False)
    n_gen: int = field(init=False)
    hodge_total: int = field(init=False)
    is_mirror_self: bool = field(init=False)

    def __post_init__(self):
        self.euler = 2 * (self.h11 - self.h21)
        self.n_gen = abs(self.h11 - self.h21)
        self.hodge_total = self.h11 + self.h21 + 2
        self.is_mirror_self = (self.h11 == self.h21)

    def to_dict(self) -> dict:
        return {
            "h11": self.h11,
            "h21": self.h21,
            "euler": self.euler,
            "n_gen": self.n_gen,
            "count": self.count,
            "source": self.source,
            "hodge_total": self.hodge_total,
            "is_mirror_self": self.is_mirror_self,
        }

    def __repr__(self):
        return (f"CY(h11={self.h11}, h21={self.h21}, χ={self.euler}, "
                f"N_gen={self.n_gen}, ×{self.count})")


class CYDatabase:
    """
    Base de données des variétés de Calabi-Yau.

    Charge les données CICY, gère la symétrie miroir,
    et fournit des méthodes de requête et de filtrage.
    """

    def __init__(self):
        self.manifolds: List[CalabiYauManifold] = []
        self._loaded_sources: List[str] = []

    # ─── Chargement ───────────────────────────────────────────────

    def load_cicy(self, include_mirrors: bool = True) -> "CYDatabase":
        """
        Charge la base CICY complète (7890 variétés, ~266 paires de Hodge distinctes).

        Args:
            include_mirrors : si True, ajoute aussi les variétés miroir (h11 ↔ h21)
                             pour compléter le Hodge shield.
        """
        from cy_landscape.data.cicy_hodge_data import CICY_HODGE_PAIRS

        seen_pairs = set()
        loaded = 0

        for h11, h21, count in CICY_HODGE_PAIRS:
            pair = (h11, h21)
            if pair not in seen_pairs:
                self.manifolds.append(
                    CalabiYauManifold(h11=h11, h21=h21, count=count, source="CICY")
                )
                seen_pairs.add(pair)
                loaded += count

            # Ajouter le miroir si demandé
            if include_mirrors and h11 != h21:
                mirror_pair = (h21, h11)
                if mirror_pair not in seen_pairs:
                    self.manifolds.append(
                        CalabiYauManifold(h11=h21, h21=h11, count=count, source="CICY-mirror")
                    )
                    seen_pairs.add(mirror_pair)

        self._loaded_sources.append("CICY")
        return self  # chainage fluide

    # ─── Statistiques ─────────────────────────────────────────────

    @property
    def total_manifolds(self) -> int:
        """Nombre total de variétés (comptées avec multiplicité)."""
        return sum(m.count for m in self.manifolds)

    @property
    def total_topologies(self) -> int:
        """Nombre de paires (h11, h21) distinctes."""
        return len(self.manifolds)

    def summary(self) -> dict:
        """Résumé statistique de la base."""
        if not self.manifolds:
            return {"status": "vide"}

        h11_vals = [m.h11 for m in self.manifolds]
        h21_vals = [m.h21 for m in self.manifolds]
        euler_vals = [m.euler for m in self.manifolds]
        ngen_vals = [m.n_gen for m in self.manifolds]

        return {
            "sources": self._loaded_sources,
            "topologies_distinctes": self.total_topologies,
            "varietes_totales": self.total_manifolds,
            "h11_range": (min(h11_vals), max(h11_vals)),
            "h21_range": (min(h21_vals), max(h21_vals)),
            "euler_range": (min(euler_vals), max(euler_vals)),
            "n_gen_range": (min(ngen_vals), max(ngen_vals)),
            "auto_miroir": sum(1 for m in self.manifolds if m.is_mirror_self),
        }

    # ─── Filtrage topologique ─────────────────────────────────────

    def filter(
        self,
        n_gen: Optional[int] = None,
        h11_range: Optional[Tuple[int, int]] = None,
        h21_range: Optional[Tuple[int, int]] = None,
        euler_range: Optional[Tuple[int, int]] = None,
        source: Optional[str] = None,
        mirror_self_only: bool = False,
    ) -> List[CalabiYauManifold]:
        """
        Filtre les variétés selon des critères topologiques.

        Args:
            n_gen : nombre de générations exact (ex: 3 pour le Modèle Standard)
            h11_range : intervalle [min, max] pour h^{1,1}
            h21_range : intervalle [min, max] pour h^{2,1}
            euler_range : intervalle [min, max] pour χ
            source : filtrer par source ("CICY", "CICY-mirror", etc.)
            mirror_self_only : si True, ne retourner que les variétés auto-miroir

        Returns:
            Liste des CalabiYauManifold satisfaisant tous les critères.
        """
        results = self.manifolds

        if n_gen is not None:
            results = [m for m in results if m.n_gen == n_gen]

        if h11_range is not None:
            lo, hi = h11_range
            results = [m for m in results if lo <= m.h11 <= hi]

        if h21_range is not None:
            lo, hi = h21_range
            results = [m for m in results if lo <= m.h21 <= hi]

        if euler_range is not None:
            lo, hi = euler_range
            results = [m for m in results if lo <= m.euler <= hi]

        if source is not None:
            results = [m for m in results if m.source == source]

        if mirror_self_only:
            results = [m for m in results if m.is_mirror_self]

        return results

    def filter_standard_model_candidates(self) -> List[CalabiYauManifold]:
        """
        Filtre spécialisé : variétés compatibles avec 3 générations de fermions.

        En compactification hétérotique E8×E8, le nombre de générations
        est |χ|/2 = |h11 - h21|.

        Pour reproduire les 3 familles du Modèle Standard : |h11 - h21| = 3.
        """
        return self.filter(n_gen=3)

    # ─── Distribution ─────────────────────────────────────────────

    def generation_distribution(self) -> dict:
        """
        Distribution du nombre de générations dans la base.
        Retourne {n_gen: (nombre_topologies, nombre_varietes)}.
        """
        dist = {}
        for m in self.manifolds:
            if m.n_gen not in dist:
                dist[m.n_gen] = {"topologies": 0, "varietes": 0}
            dist[m.n_gen]["topologies"] += 1
            dist[m.n_gen]["varietes"] += m.count
        return dict(sorted(dist.items()))

    # ─── Export ────────────────────────────────────────────────────

    def to_json(self, filepath: str):
        """Exporte la base en JSON."""
        data = {
            "metadata": self.summary(),
            "manifolds": [m.to_dict() for m in self.manifolds],
        }
        with open(filepath, "w") as f:
            json.dump(data, f, indent=2)

    def to_csv(self, filepath: str):
        """Exporte la base en CSV."""
        import csv
        with open(filepath, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["h11", "h21", "euler", "n_gen", "count", "source",
                             "hodge_total", "is_mirror_self"])
            for m in self.manifolds:
                writer.writerow([m.h11, m.h21, m.euler, m.n_gen, m.count,
                                 m.source, m.hodge_total, m.is_mirror_self])
