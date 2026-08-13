"""
extensions.py -- Fibres d'extension sur les CICYs

Construction : 0 -> F1 -> V -> F2 -> 0  (suite exacte courte)
  F1, F2 : line bundles ou sommes de line bundles
  V      : fibre d'extension, rk(V) = rk(F1) + rk(F2)

L'extension existe si Ext^1(F2, F1) = H^1(F1 ⊗ F2*) ≠ 0.

Cohomologie de V : suite exacte longue directe
  0 -> H⁰(F1) -> H⁰(V) -> H⁰(F2) -> H¹(F1) -> H¹(V) -> H¹(F2) -> ...

Avantage par rapport aux monades :
  - Les sous-faisceaux de V sont plus contraints
  - Si mu(F1) = mu(F2) = 0 et l'extension est non-triviale,
    V est souvent semi-stable

Ref: Friedman, Morgan, Witten "Vector Bundles and F Theory" (1997)
     Anderson et al. "Heterotic Compactification" (2011)
"""
import numpy as np
from cy_landscape.core.positive_monads import _mix
from typing import List, Optional
from dataclasses import dataclass

from cy_landscape.core.exact_cohomology import koszul_cohomology


@dataclass
class ExtensionBundle:
    """Un fibre d'extension 0 -> F1 -> V -> F2 -> 0."""
    f1_charges: List[List[int]]  # Charges de F1 (liste de vecteurs)
    f2_charges: List[List[int]]  # Charges de F2

    @property
    def rank_F1(self): return len(self.f1_charges)

    @property
    def rank_F2(self): return len(self.f2_charges)

    @property
    def rank_V(self): return self.rank_F1 + self.rank_F2

    @property
    def m(self): return len(self.f1_charges[0])

    @property
    def c1_F1(self):
        return [sum(f[i] for f in self.f1_charges) for i in range(self.m)]

    @property
    def c1_F2(self):
        return [sum(f[i] for f in self.f2_charges) for i in range(self.m)]

    @property
    def c1_V(self):
        return [self.c1_F1[i] + self.c1_F2[i] for i in range(self.m)]

    @property
    def c1_vanishes(self):
        return all(x == 0 for x in self.c1_V)


def check_extension_exists(ext, ambient, config):
    """
    Verifie que Ext^1(F2, F1) = H^1(F1 ⊗ F2*) ≠ 0.

    Pour F1 = ⊕ O(a_i), F2 = ⊕ O(b_j) :
    F1 ⊗ F2* = ⊕_{i,j} O(a_i - b_j)
    H^1(F1 ⊗ F2*) = ⊕_{i,j} H^1(O(a_i - b_j))
    """
    total_h1 = 0
    for a in ext.f1_charges:
        for b in ext.f2_charges:
            diff = [a[k] - b[k] for k in range(ext.m)]
            h = koszul_cohomology(ambient, config, diff)
            total_h1 += h.get(1, 0)

    return total_h1 > 0, total_h1


def compute_extension_cohomology(ext, ambient, config):
    """
    H^i(V) via la suite exacte longue 0 -> F1 -> V -> F2 -> 0.

    H^i(F1) = ⊕_k H^i(O(a_k))
    H^i(F2) = ⊕_k H^i(O(b_k))

    Suite exacte longue :
    0 -> H⁰(F1) -> H⁰(V) -> H⁰(F2) -> H¹(F1) -> H¹(V) -> H¹(F2)
      -> H²(F1) -> H²(V) -> H²(F2) -> H³(F1) -> H³(V) -> H³(F2) -> 0

    Pour une extension GENERIQUE (non-scindee) :
    Les connecting morphisms sont de rang maximal.
    """
    hF1 = {i: 0 for i in range(4)}
    for a in ext.f1_charges:
        h = koszul_cohomology(ambient, config, a)
        for i in range(4):
            hF1[i] += h.get(i, 0)

    hF2 = {i: 0 for i in range(4)}
    for b in ext.f2_charges:
        h = koszul_cohomology(ambient, config, b)
        for i in range(4):
            hF2[i] += h.get(i, 0)

    # Suite exacte longue (extension generique)
    # 0 -> H⁰(F1) -> H⁰(V) -> H⁰(F2) -delta-> H¹(F1) -> H¹(V) -> H¹(F2) -> ...
    # delta est de rang maximal pour une extension generique

    # H⁰: 0 -> H⁰(F1) -> H⁰(V) -> ker(delta) -> 0
    rank_d0 = min(hF2[0], hF1[1])  # connecting: H⁰(F2) -> H¹(F1)
    h0V = hF1[0] + hF2[0] - rank_d0

    # H³: coker -> H³(V) -> H³(F2) -> 0
    rank_d2 = min(hF2[2], hF1[3])  # connecting: H²(F2) -> H³(F1)
    h3V = hF1[3] - rank_d2 + hF2[3]

    # chi(V) = chi(F1) + chi(F2)
    chiF1 = sum((-1)**i * hF1[i] for i in range(4))
    chiF2 = sum((-1)**i * hF2[i] for i in range(4))
    chiV = chiF1 + chiF2

    # H¹: coker(delta_0) -> H¹(V) -> ker(delta_1)
    coker_d0 = hF1[1] - rank_d0
    rank_d1 = min(hF2[1], hF1[2])
    h1V = coker_d0 + hF2[1] - rank_d1
    h1V = max(0, h1V)

    h2V = h1V - chiV + h0V - h3V
    h2V = max(0, h2V)

    return {0: max(0, int(h0V)), 1: max(0, int(h1V)),
            2: max(0, int(h2V)), 3: max(0, int(h3V)),
            'hF1': dict(hF1), 'hF2': dict(hF2), 'chi_V': int(chiV)}


def generate_extensions(m, rank_V, max_charge=3, n_random=200, rng=None,
                        seed=42):
    """
    Genere des fibres d'extension candidats.

    Types :
    - (1, rk-1) : F1 = O(a), F2 = somme de rk-1 line bundles
    - (2, rk-2) : F1 = 2 line bundles, F2 = rk-2 line bundles
    - (rk-1, 1) : inverse
    """
    # `rng` n'est plus utilise pour les tirages : conserve pour compatibilite
    # d'appel. Voir ci-dessous.
    extensions = []

    for rk1 in range(1, rank_V):
        rk2 = rank_V - rk1
        if rk1 > rk2:
            break  # Symetrie: (rk1, rk2) ~ (rk2, rk1) par dualite

        # RNG PROPRE AU COUPLE (rk1, rk2), derive de (seed, rk1, rk2, m,
        # max_charge). Le generateur puisait dans un RNG partage entre les
        # valeurs de rk1 ET avec l'appelant : allonger la premiere boucle
        # decalait les tirages de la seconde. Mesure de l'audit : passer de
        # n_random = 200 a 800 FAISAIT PERDRE 40 extensions au rang 4 et 31
        # au rang 5. Un generateur doit etre monotone en son propre budget.
        #
        # RESERVE : ceci retablit la monotonie en `n_random`, PAS en
        # `max_charge` -- changer les bornes du tirage change la suite des
        # valeurs tirees, et l'audit mesurait 200 pertes sur 206 en passant
        # de max_charge 2 a 3. Seule une enumeration corrigerait ce
        # second point ; elle reste a faire.
        rng_loc = np.random.RandomState(
            _mix(seed, (m, max_charge, rk1, rk2), rank_V))

        for _ in range(n_random // max(1, rank_V - 1)):
            f1 = []
            for r in range(rk1):
                q = [int(rng_loc.randint(-max_charge, max_charge + 1)) for _ in range(m)]
                f1.append(q)

            # F2 : c1(F2) = -c1(F1)
            c1_f1 = [sum(f[k] for f in f1) for k in range(m)]
            f2 = []
            for r in range(rk2 - 1):
                q = [int(rng_loc.randint(-max_charge, max_charge + 1)) for _ in range(m)]
                f2.append(q)
            # Derniere charge : forcer c1(V) = 0
            last = [-c1_f1[k] - sum(f[k] for f in f2) for k in range(m)]
            if all(abs(l) <= max_charge + 2 for l in last):
                f2.append(last)
                ext = ExtensionBundle(f1, f2)
                if ext.c1_vanishes:
                    extensions.append(ext)

        # Extensions structurees : F1 = O(a), F2 = O(-a) ⊕ O(0)^{rk-2}
        if rk1 == 1:
            for i in range(m):
                for sign in [1, -1]:
                    a = [0] * m; a[i] = sign
                    neg_a = [-x for x in a]
                    f2 = [neg_a] + [[0]*m for _ in range(rk2 - 1)]
                    ext = ExtensionBundle([a], f2)
                    if ext.c1_vanishes:
                        extensions.append(ext)

    return extensions


# ======================================================================
# CHEMIN CORRECT POUR LE FIBRE D'EXTENSION  (repare le defaut 4.7)
# ======================================================================
#
# Le pipeline construisait une PSEUDO-MONADE B = F1 (+) F2, C = F2 pour
# reutiliser le chemin des monades. Le noyau de F1 (+) F2 -> F2 est de rang
# rank(F1) et de caracteristique chi(F1), alors que le fibre d'extension
# 0 -> F1 -> V -> F2 -> 0 est de rang rank(F1) + rank(F2) et de
# caracteristique chi(F1) + chi(F2). Cohomologie, Hoppe et groupe de jauge
# portaient donc sur un AUTRE objet -- 1571 entrees sur 1571 en incoherence
# de rang sur le scan test_v3.
#
# Les deux fonctions ci-dessous n'utilisent plus la pseudo-monade.

def chi_extension(ext, chical):
    """
    chi(V) EXACT pour l'extension. chi est additif sur les suites exactes :
    chi(V) = chi(F1) + chi(F2). Aucune hypothese, aucune degenerescence.

    A comparer a ce que donnait la pseudo-monade : chi(F1 (+) F2) - chi(F2)
    = chi(F1). C'est l'erreur numerique concrete du defaut 4.7.
    """
    return (chical.bundle(ext.f1_charges) + chical.bundle(ext.f2_charges))


def _charges_wedge(charges, a):
    """Charges de wedge^a (+) O(c_i) : sommes sur les a-uplets croissants."""
    from itertools import combinations
    if a == 0:
        return [[0] * len(charges[0])]
    m = len(charges[0])
    return [[sum(charges[i][k] for i in idx) for k in range(m)]
            for idx in combinations(range(len(charges)), a)]


def hoppe_extension(ext, ambient, config):
    """
    Critere de Hoppe pour le fibre d'extension, par BORNE SUPERIEURE.

    ----------------------------------------------------------------------
    Principe
    ----------------------------------------------------------------------
    wedge^p V admet une filtration dont les quotients gradues sont les
    wedge^a F1 (x) wedge^b F2 avec a + b = p. La suite exacte longue donne
    alors

        h^0(wedge^p V)  <=  somme_{a+b=p} h^0(wedge^a F1 (x) wedge^b F2)

    F1 et F2 etant des sommes de fibres en droites, chaque terme est une
    somme de h^0 de fibres en droites, calculables exactement par Koszul.

    Si toutes ces sommes sont NULLES, alors h^0(wedge^p V) = 0 pour tout
    p = 1..rk-1, et le critere de Hoppe est satisfait : V est stable. La
    condition est SUFFISANTE et ne peut pas donner de faux positif.

    Si une somme est non nulle, on ne conclut pas : les sections des
    quotients gradues ne se relevent pas necessairement a V. L'etat est
    alors 'indetermine', jamais 'stable'.

    Chaque h^0 doit etre CERTIFIE par `koszul_cohomology_ex` ; un seul degre
    non certifie rend le verdict indetermine.
    """
    from cy_landscape.core.exact_cohomology import koszul_cohomology_ex
    rk = ext.rank_V
    sommes, non_certifie = {}, []
    for p in range(1, rk):
        total = 0
        for a in range(0, min(p, ext.rank_F1) + 1):
            b = p - a
            if b > ext.rank_F2:
                continue
            for u in _charges_wedge(ext.f1_charges, a):
                for v in _charges_wedge(ext.f2_charges, b):
                    ch = [u[k] + v[k] for k in range(ext.m)]
                    r = koszul_cohomology_ex(ambient, config, ch)
                    if not r['certified_by_degree'][0]:
                        non_certifie.append((p, tuple(ch)))
                        continue
                    total += r[0]
        sommes[p] = total
        if total > 0:
            return {'stable': None, 'etat': f'indetermine : borne h0(w^{p}V) = {total}',
                    'bornes': sommes, 'non_certifie': non_certifie}
    if non_certifie:
        return {'stable': None, 'etat': 'indetermine : degres non certifies',
                'bornes': sommes, 'non_certifie': non_certifie}
    return {'stable': True, 'etat': 'stable (toutes les bornes nulles)',
            'bornes': sommes, 'non_certifie': []}
