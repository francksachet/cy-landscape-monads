"""
cohomology.py -- Calcul de la cohomologie des fibres et extraction du spectre.

En compactification heterotique, le spectre de particules a 4D provient
des groupes de cohomologie du fibre vectoriel V sur le CY Y :

    Cohomologie              ->   Representation du groupe de jauge
    H1(Y, V)                ->   Particules chirales (type 1)
    H1(Y, V*)               ->   Particules anti-chirales (type 1)
    H1(Y, wedge2 V)          ->   Particules chirales (type 2)
    H1(Y, End V)             ->   Moduli du fibre (scalaires neutres)

Le contenu depend du groupe de jauge residuel :

  SU(5) GUT : 10 <- H1(V),  5bar <- H1(wedge2 V),  1 <- H1(End V)
  SO(10)    : 16 <- H1(V),  10 <- H1(wedge2 V),  1 <- H1(End V)
  E6        : 27 <- H1(V),  27bar <- H1(V*),  1 <- H1(End V)
"""

from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
import numpy as np

from cy_landscape.core.bundles import (
    CICYGeometry, LineBundleSum, BundleAnalysis, GAUGE_GROUP_TABLE,
)


# =================================================================
# Contenu en particules du Modele Standard (cible)
# =================================================================

@dataclass
class SMParticleContent:
    """Contenu en particules du Modele Standard, comme reference."""
    quarks_up: int = 3
    quarks_down: int = 3
    leptons_charged: int = 3
    neutrinos: int = 3
    higgs_doublets: int = 1

    @property
    def n_generations(self) -> int:
        return 3

    @property
    def n_10(self) -> int:
        return 3

    @property
    def n_5bar(self) -> int:
        return 3 + 1  # 3 generations + 1 Higgs 5bar

    @property
    def n_5(self) -> int:
        return 1  # 1 Higgs 5


SM_TARGET = SMParticleContent()


# =================================================================
# Cohomologie de fibres en droites
# =================================================================

def line_bundle_cohomology(
    charges: List[int],
    geometry: CICYGeometry,
) -> Dict[int, int]:
    """
    Calcule h^i(Y, L) pour un fibre en droites L.
    Approximation basee sur Hirzebruch-Riemann-Roch + Bott-Borel-Weil.
    """
    h = geometry.n_kahler
    d = geometry.intersection_numbers
    charges = list(charges)

    # Degre
    degree = 0.0
    for i in range(h):
        for j in range(h):
            for k in range(h):
                degree += charges[i] * d[i, j, k]
    degree /= 6.0

    # Index chi(L) par HRR
    chi_L = 0.0
    for i in range(h):
        for j in range(h):
            for k in range(h):
                chi_L += charges[i] * charges[j] * charges[k] * d[i, j, k]
    chi_L /= 6.0

    c2_term = 0.0
    for i in range(h):
        c2_term += geometry.c2_tangent[i] * charges[i]
    chi_L += c2_term / 12.0
    chi_L = round(chi_L)

    n_pos = sum(1 for a in charges if a > 0)
    n_neg = sum(1 for a in charges if a < 0)

    cohom = {0: 0, 1: 0, 2: 0, 3: 0}

    if n_neg == 0 and degree > 0:
        cohom[0] = max(0, int(chi_L))
    elif n_pos == 0 and degree < 0:
        cohom[3] = max(0, int(-chi_L))
    elif degree > 1:
        cohom[0] = max(0, int(chi_L))
        if chi_L < 0:
            cohom[1] = int(-chi_L)
            cohom[0] = 0
    elif degree < -1:
        cohom[3] = max(0, int(-chi_L))
        if chi_L > 0:
            cohom[2] = int(chi_L)
            cohom[3] = 0
    else:
        if chi_L > 0:
            cohom[0] = int(chi_L)
        elif chi_L < 0:
            cohom[1] = int(-chi_L)

    return cohom


# =================================================================
# Cohomologie d'un fibre vectoriel (somme de line bundles)
# =================================================================

def bundle_cohomology(
    bundle: LineBundleSum,
    geometry: CICYGeometry,
) -> Dict[str, Dict[int, int]]:
    """
    Calcule la cohomologie de V, V*, wedge2 V, End V.
    """
    n = bundle.rank
    h = geometry.n_kahler

    # H^i(V) = direct_sum H^i(Lj)
    cohom_V = {0: 0, 1: 0, 2: 0, 3: 0}
    for j in range(n):
        coh_j = line_bundle_cohomology(bundle.charges[j], geometry)
        for i in range(4):
            cohom_V[i] += coh_j[i]

    # H^i(V*) = direct_sum H^i(Lj*)
    cohom_Vdual = {0: 0, 1: 0, 2: 0, 3: 0}
    for j in range(n):
        dual_charges = [-a for a in bundle.charges[j]]
        coh_j = line_bundle_cohomology(dual_charges, geometry)
        for i in range(4):
            cohom_Vdual[i] += coh_j[i]

    # H^i(wedge2 V) = direct_sum_{j<k} H^i(Lj tensor Lk)
    cohom_wedge2 = {0: 0, 1: 0, 2: 0, 3: 0}
    for j in range(n):
        for k in range(j + 1, n):
            prod_charges = [
                bundle.charges[j][a] + bundle.charges[k][a]
                for a in range(len(bundle.charges[j]))
            ]
            coh_jk = line_bundle_cohomology(prod_charges, geometry)
            for i in range(4):
                cohom_wedge2[i] += coh_jk[i]

    # H^i(End V) = direct_sum_{j,k} H^i(Lj tensor Lk*)
    cohom_end = {0: 0, 1: 0, 2: 0, 3: 0}
    for j in range(n):
        for k in range(n):
            end_charges = [
                bundle.charges[j][a] - bundle.charges[k][a]
                for a in range(len(bundle.charges[j]))
            ]
            coh_jk = line_bundle_cohomology(end_charges, geometry)
            for i in range(4):
                cohom_end[i] += coh_jk[i]

    return {
        "V": cohom_V,
        "V_dual": cohom_Vdual,
        "wedge2V": cohom_wedge2,
        "end_V": cohom_end,
    }


# =================================================================
# Spectre de particules
# =================================================================

@dataclass
class ParticleSpectrum:
    """Spectre de particules extrait de la cohomologie."""
    gauge_group: str = ""
    representations: Dict[str, int] = field(default_factory=dict)
    n_generations: int = 0
    n_anti_generations: int = 0
    n_higgs_candidates: int = 0
    n_singlets: int = 0
    n_exotics: int = 0
    generation_match: bool = False
    higgs_present: bool = False
    exotic_free: bool = False
    sm_compatibility: float = 0.0

    def compute_sm_compatibility(self):
        score = 0.0

        # 3 generations (40 pts)
        if self.n_generations == 3:
            score += 40.0
            self.generation_match = True
        elif self.n_generations in [2, 4]:
            score += 15.0
        elif self.n_generations == 1:
            score += 5.0

        # Higgs (20 pts)
        if self.n_higgs_candidates >= 1:
            score += 20.0
            self.higgs_present = True
            if self.n_higgs_candidates == 1:
                score += 5.0

        # Pas d'exotiques (25 pts)
        if self.n_exotics == 0:
            score += 25.0
            self.exotic_free = True
        elif self.n_exotics <= 2:
            score += 10.0

        # Singlets moderes (10 pts)
        if 1 <= self.n_singlets <= 20:
            score += 10.0
        elif self.n_singlets <= 50:
            score += 5.0

        self.sm_compatibility = min(100.0, score)

    def to_dict(self) -> dict:
        return {
            "gauge_group": self.gauge_group,
            "representations": dict(self.representations),
            "n_generations": int(self.n_generations),
            "n_anti_generations": int(self.n_anti_generations),
            "n_higgs_candidates": int(self.n_higgs_candidates),
            "n_singlets": int(self.n_singlets),
            "n_exotics": int(self.n_exotics),
            "generation_match": bool(self.generation_match),
            "higgs_present": bool(self.higgs_present),
            "exotic_free": bool(self.exotic_free),
            "sm_compatibility": float(round(self.sm_compatibility, 1)),
        }


def extract_spectrum_su5(cohom):
    sp = ParticleSpectrum(gauge_group="SU(5)")
    n_10 = cohom["V"].get(1, 0)
    n_10bar = cohom["V_dual"].get(1, 0)
    n_5bar = cohom["wedge2V"].get(1, 0)
    n_5 = cohom["wedge2V"].get(2, 0)
    n_singlets = cohom["end_V"].get(1, 0)

    sp.representations = {
        "10": n_10, "10bar": n_10bar,
        "5bar": n_5bar, "5": n_5, "1": n_singlets,
    }
    sp.n_generations = abs(n_10 - n_10bar)
    sp.n_anti_generations = min(n_10, n_10bar)
    gen_5bar = min(n_5bar, sp.n_generations)
    excess_5bar = max(0, n_5bar - gen_5bar)
    sp.n_higgs_candidates = min(n_5, excess_5bar)
    sp.n_singlets = n_singlets
    sp.n_exotics = max(0, n_10 + n_10bar - sp.n_generations - sp.n_anti_generations * 2)
    sp.compute_sm_compatibility()
    return sp


def extract_spectrum_so10(cohom):
    sp = ParticleSpectrum(gauge_group="SO(10)")
    n_16 = cohom["V"].get(1, 0)
    n_16bar = cohom["V_dual"].get(1, 0)
    n_10 = cohom["wedge2V"].get(1, 0)
    n_singlets = cohom["end_V"].get(1, 0)

    sp.representations = {
        "16": n_16, "16bar": n_16bar, "10": n_10, "1": n_singlets,
    }
    sp.n_generations = abs(n_16 - n_16bar)
    sp.n_anti_generations = min(n_16, n_16bar)
    sp.n_higgs_candidates = n_10
    sp.n_singlets = n_singlets
    sp.n_exotics = 0
    sp.compute_sm_compatibility()
    return sp


def extract_spectrum_e6(cohom):
    sp = ParticleSpectrum(gauge_group="E6")
    n_27 = cohom["V"].get(1, 0)
    n_27bar = cohom["V_dual"].get(1, 0)
    n_singlets = cohom["end_V"].get(1, 0)

    sp.representations = {"27": n_27, "27bar": n_27bar, "1": n_singlets}
    sp.n_generations = abs(n_27 - n_27bar)
    sp.n_anti_generations = min(n_27, n_27bar)
    sp.n_higgs_candidates = max(0, sp.n_generations - 3) + sp.n_anti_generations
    sp.n_singlets = n_singlets + sp.n_generations
    sp.n_exotics = sp.n_anti_generations
    sp.compute_sm_compatibility()
    return sp


def extract_spectrum(bundle_analysis):
    """Point d'entree principal : extrait le spectre selon le groupe de jauge."""
    cohom = bundle_cohomology(bundle_analysis.bundle, bundle_analysis.geometry)
    if bundle_analysis.gauge_group == "SU(5)":
        return extract_spectrum_su5(cohom)
    elif bundle_analysis.gauge_group == "SO(10)":
        return extract_spectrum_so10(cohom)
    elif "E" in bundle_analysis.gauge_group:
        return extract_spectrum_e6(cohom)
    else:
        sp = ParticleSpectrum(gauge_group=bundle_analysis.gauge_group)
        sp.compute_sm_compatibility()
        return sp


# =================================================================
# Analyse complete : fibre + spectre
# =================================================================

@dataclass
class FullAnalysis:
    """Resultat complet : geometrie + fibre + spectre + score final."""
    bundle_analysis: BundleAnalysis
    spectrum: ParticleSpectrum
    final_score: float = 0.0

    def compute_final_score(self) -> float:
        self.final_score = (
            0.4 * self.bundle_analysis.sm_score +
            0.6 * self.spectrum.sm_compatibility
        )
        return self.final_score

    def to_dict(self) -> dict:
        return {
            "bundle": self.bundle_analysis.to_dict(),
            "spectrum": self.spectrum.to_dict(),
            "final_score": float(round(self.final_score, 1)),
        }

    def summary(self) -> str:
        ba = self.bundle_analysis
        sp = self.spectrum
        ck = lambda b: "Y" if b else "N"
        lines = [
            f"=== CY({ba.geometry.h11},{ba.geometry.h21}) + "
            f"SU({ba.structure_group_rank}) -> {ba.gauge_group} ===",
            f"",
            f"  Topologie : h11={ba.geometry.h11}, h21={ba.geometry.h21}, "
            f"chi={ba.geometry.euler}",
            f"  Fibre     : rang {ba.structure_group_rank}, c1=0: {ck(ba.c1_vanishes)}",
            f"  Jauge     : {ba.gauge_group}",
            f"",
            f"  Spectre de matiere :",
        ]
        for rep, count in sp.representations.items():
            if count > 0:
                lines.append(f"    {count:>3} x {rep}")
        lines.extend([
            f"",
            f"  Physique :",
            f"    Generations    : {sp.n_generations} {ck(sp.generation_match)}",
            f"    Anti-gen.      : {sp.n_anti_generations}",
            f"    Higgs candidats: {sp.n_higgs_candidates} {ck(sp.higgs_present)}",
            f"    Singlets       : {sp.n_singlets}",
            f"    Exotiques      : {sp.n_exotics} {ck(sp.exotic_free)}",
            f"",
            f"  Scores :",
            f"    Fibre  (etage 2) : {ba.sm_score:.0f}/100",
            f"    Spectre (etage 3): {sp.sm_compatibility:.0f}/100",
            f"    * Final combine  : {self.final_score:.0f}/100",
        ])
        return "\n".join(lines)
