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


# ======================================================================
# ENUMERATION EXHAUSTIVE  (remplace l'echantillonnage, defaut 5.11)
# ======================================================================
#
# L'echantillonnage n'etait monotone NI en `n_random` (corrige par un RNG
# derive par couple (rk1, rk2)) NI en `max_charge` : passer de max_charge
# 2 a 3 faisait PERDRE 216 extensions sur 222, soit 97 %. Ce second point
# ne peut pas se corriger par un reglage de RNG -- changer les bornes du
# tirage change la suite tiree.
#
# Un generateur non monotone en son propre parametre de portee interdit
# tout enonce d'ABSENCE : « aucun fibre stable a max_charge 3 » ne dit
# rien tant qu'on ignore si le domaine a max_charge 3 contient celui de
# max_charge 2. L'enumeration rend cette inclusion vraie PAR CONSTRUCTION,
# et non par chance de RNG.
#
# Domaine enumere : toutes les charges de F1 et de F2 dans la boite
# [-max_charge, max_charge]^m, la DERNIERE charge etant determinee par
# c1(V) = 0. C'est un produit de boites, donc croissant en max_charge.
#
# NOTE sur le domaine du tirage. Celui-ci acceptait la derniere charge
# jusqu'a max_charge + 2 alors que toutes les autres etaient tirees dans
# [-max_charge, max_charge]. Cette asymetrie etait un artefact de la
# construction « tirer n-1 vecteurs, deduire le dernier » : elle n'a
# aucune justification geometrique et fait dependre le domaine de l'ORDRE
# des facteurs de F2. Le tirage de secours est aligne sur la boite
# symetrique, faute de quoi l'inclusion « tirage c enumeration » ne
# serait pas testable -- et c'est elle qui autorise a remplacer « aucun
# survivant parmi ce qu'on a tire » par « aucun survivant sur le
# domaine ».


def _splits(rank_V):
    """
    Decoupages (rk1, rk2) avec rk1 + rk2 = rank_V et rk1 <= rk2.

    On s'arrete a rk1 <= rk2 : le dual de 0 -> F1 -> V -> F2 -> 0 est
    0 -> F2* -> V* -> F1* -> 0, et la boite [-q, q]^m est stable par
    negation. Le decoupage (rk2, rk1) est donc exactement l'image du
    decoupage (rk1, rk2) par dualite, et V* est stable si et seulement si
    V l'est.
    """
    return [(rk1, rank_V - rk1) for rk1 in range(1, rank_V)
            if rk1 <= rank_V - rk1]


def _boite(m, q):
    """Tous les vecteurs de charge de [-q, q]^m."""
    from itertools import product
    return [tuple(v) for v in product(range(-q, q + 1), repeat=m)]


def _compte_sommes_bornees(n, q, borne):
    """
    #{(x_1..x_n) dans [-q,q]^n : |x_1 + ... + x_n| <= borne}, par
    convolution. Ne construit AUCUN tuple : c'est une reference
    independante du chemin d'enumeration, exactement comme `compte_B`
    l'est de `enumerer_positive_B`.
    """
    dist = {0: 1}
    for _ in range(n):
        suiv = {}
        for s, k in dist.items():
            for x in range(-q, q + 1):
                suiv[s + x] = suiv.get(s + x, 0) + k
        dist = suiv
    return sum(k for s, k in dist.items() if abs(s) <= borne)


def compte_extensions(m, rank_V, max_charge):
    """
    Nombre de tuples ORDONNES du domaine, sans rien enumerer.

    C'est un MAJORANT du nombre d'extensions distinctes : deux tuples qui
    ne different que par l'ordre des facteurs de F1 (ou de F2) donnent le
    MEME fibre. Sert de garde-fou avant d'enumerer, comme `compte_B`.

    Les m composantes sont independantes -- la contrainte c1(V) = 0 se lit
    composante par composante -- d'ou la puissance m-ieme.
    """
    n_libres = rank_V - 1          # le dernier vecteur est determine
    par_split = _compte_sommes_bornees(n_libres, max_charge, max_charge) ** m
    return len(_splits(rank_V)) * par_split


def _multisets_somme(boite, taille, cible, q, i_min=0):
    """
    GENERATEUR des multiensembles de `taille` vecteurs de `boite` dont la
    somme vaut `cible`. Chaque multiensemble est produit UNE SEULE FOIS
    (indices non decroissants), donc aucun ensemble de deduplication n'est
    conserve en memoire -- c'est exactement le piege qui avait produit un
    MemoryError apres 153 min de scan sur `enumerer_positive_B` (§5.11).

    Elagage EXACT : les `taille` vecteurs restants contribuent au plus
    taille*q par composante ; si |cible_k| depasse ce montant, aucune
    completion n'existe et la branche est coupee sans perte.
    """
    m = len(cible)
    if taille == 0:
        if all(x == 0 for x in cible):
            yield []
        return
    if any(abs(x) > taille * q for x in cible):
        return
    for i in range(i_min, len(boite)):
        v = boite[i]
        reste = tuple(cible[k] - v[k] for k in range(m))
        for suite in _multisets_somme(boite, taille - 1, reste, q, i):
            yield [v] + suite


def _enumerer_extensions(m, rank_V, q):
    """
    GENERATEUR de toutes les extensions du domaine, sans doublon.

    F1 est libre (multiensemble de rk1 vecteurs de la boite) ; F2 est
    contraint par c1(F2) = -c1(F1). Deux decoupages (rk1, rk2) distincts
    donnent des fibres de rangs (rk1, rk2) distincts, donc aucun doublon
    inter-decoupage n'est possible : rien a dedupliquer globalement.
    """
    from itertools import combinations_with_replacement
    boite = _boite(m, q)
    for rk1, rk2 in _splits(rank_V):
        for f1 in combinations_with_replacement(boite, rk1):
            cible = tuple(-sum(v[k] for v in f1) for k in range(m))
            for f2 in _multisets_somme(boite, rk2, cible, q):
                yield ExtensionBundle([list(v) for v in f1],
                                      [list(v) for v in f2])


def enumerer_extensions(m, rank_V, max_charge=3, plafond=200_000):
    """
    TOUTES les extensions du domaine [-max_charge, max_charge]^m, ou None
    si le nombre de tuples ordonnes depasse `plafond` -- a charge de
    l'appelant de retomber sur l'echantillonnage.

    Renvoie un GENERATEUR : le plafond borne le NOMBRE d'extensions, pas
    la memoire, et ce sont deux choses differentes des lors qu'on conserve
    chaque objet.
    """
    n_ordonnes = compte_extensions(m, rank_V, max_charge)
    if n_ordonnes == 0 or n_ordonnes > plafond:
        return None
    return _enumerer_extensions(m, rank_V, max_charge)


def generate_extensions(m, rank_V, max_charge=3, n_random=200, rng=None,
                        seed=42, exhaustif_max=200_000, stats=None):
    """
    Genere des fibres d'extension candidats.

    Types :
    - (1, rk-1) : F1 = O(a), F2 = somme de rk-1 line bundles
    - (2, rk-2) : F1 = 2 line bundles, F2 = rk-2 line bundles

    Par DEFAUT le domaine est ENUMERE (exhaustif_max > 0) ; l'echantillonnage
    n'est plus qu'un recours quand le domaine depasse le plafond. Ce recours
    est non monotone en `max_charge` (§5.11) : aucun enonce d'absence ne peut
    s'appuyer dessus, d'ou le compteur ci-dessous.

    `stats`, si un dict est fourni, recoit le mode reellement employe --
    indispensable pour savoir sur quel domaine un resultat est DEMONTRE et
    sur quel domaine il n'est qu'un sondage.

    `rng` n'est plus utilise (RNG derive de `seed`) ; conserve pour
    compatibilite d'appel.

    Renvoie un ITERABLE : generateur en mode exhaustif, liste en mode tirage.
    """
    if stats is not None:
        stats.setdefault('ext_exhaustifs', 0)
        stats.setdefault('ext_echantillonnes', 0)

    if exhaustif_max > 0:
        enum = enumerer_extensions(m, rank_V, max_charge,
                                   plafond=exhaustif_max)
        if enum is not None:
            if stats is not None:
                stats['ext_exhaustifs'] += 1
            return enum

    if stats is not None:
        stats['ext_echantillonnes'] += 1
    return _echantillonner_extensions(m, rank_V, max_charge, n_random, seed)


def _echantillonner_extensions(m, rank_V, max_charge=3, n_random=200,
                               seed=42):
    """
    Ancien chemin, conserve comme SECOURS pour les domaines hors plafond.
    NON MONOTONE en `max_charge` : ne jamais en tirer un enonce d'absence.
    """
    extensions = []

    for rk1, rk2 in _splits(rank_V):
        # RNG PROPRE AU COUPLE (rk1, rk2), derive de (seed, rk1, rk2, m,
        # max_charge). Le generateur puisait dans un RNG partage entre les
        # valeurs de rk1 ET avec l'appelant : allonger la premiere boucle
        # decalait les tirages de la seconde. Mesure de l'audit : passer de
        # n_random = 200 a 800 FAISAIT PERDRE 40 extensions au rang 4 et 31
        # au rang 5. Un generateur doit etre monotone en son propre budget.
        #
        # RESERVE : ceci retablit la monotonie en `n_random`, PAS en
        # `max_charge` -- changer les bornes du tirage change la suite des
        # valeurs tirees, et l'audit mesurait 216 pertes sur 222 en passant
        # de max_charge 2 a 3. C'est `enumerer_extensions` qui regle ce
        # second point ; ce chemin-ci ne sert plus que hors plafond.
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
            # Derniere charge : forcer c1(V) = 0. Borne alignee sur la
            # boite symetrique (etait max_charge + 2, cf. la note ci-dessus)
            # pour que tirage et enumeration portent sur le meme domaine.
            last = [-c1_f1[k] - sum(f[k] for f in f2) for k in range(m)]
            if all(abs(l) <= max_charge for l in last):
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


# ======================================================================
# PENTE : les sous-faisceaux DESTABILISANTS, lus sur la filtration
# ======================================================================
#
# Pourquoi ce test existe, et pourquoi il vient AVANT Hoppe.
#
# Le critere de Hoppe s'enonce « c1(V) = 0 => V stable <=> h0(w^p V) = 0
# pour p = 1..rk-1 ». Cette EQUIVALENCE suppose Pic(X) de rang 1. Sur une
# CICY a m > 1 facteurs, la stabilite de pente depend de la classe de
# Kahler J choisie, et la condition de Hoppe reste NECESSAIRE sans etre
# suffisante : elle est aveugle a la polarisation.
#
# Pour une monade, aucun sous-faisceau ne se presente commodement et l'on
# ne peut pas faire mieux a bon compte. Pour une extension, si :
# 0 -> F1 -> V -> F2 -> 0 exhibe F1 comme sous-faisceau de V, et avec lui
# toutes ses sous-sommes, ainsi que les preimages des sous-sommes de F2.
# La stabilite exige deg_J(W) < 0 pour chacun -- de l'arithmetique pure
# sur les d_ijk, au meme cout que le prefiltre chi.
#
# ----------------------------------------------------------------------
# UN PIEGE, RENCONTRE EN ECRIVANT CE MODULE  -- a ne pas refaire
# ----------------------------------------------------------------------
# Premiere version : « cherchons un J sur la grille [1,4]^m ; si aucun ne
# rend tous les degres negatifs, le fibre est instable ». Elle annoncait
# 635 extensions destabilisees sur 2 647, soit 24 %.
#
# CE CHIFFRE N'EXISTE PAS. Il mesurait la grille, pas les fibres :
#
#   J_max  |  3     6     12    24
#   sans temoin | 1748  1299  1042   925
#
# Aucune saturation. Le « sans temoin » suit le budget de recherche, comme
# le faux lieu de base du §5.4 suivait la dimension de la source. Un echec
# de recherche sur une grille finie ne demontre rien, et ne doit JAMAIS
# etre inscrit comme une elimination.
#
# D'ou la separation stricte ci-dessous entre ce qui est DEMONTRE et ce
# qui n'est que cherche. Trois issues :
#   False -> instable, DEMONTRE par un certificat exact ;
#   True  -> un temoin J rend tous ces degres negatifs -- condition
#            NECESSAIRE satisfaite pour ce J, ce qui ne prouve pas la
#            stabilite ;
#   None  -> ni certificat, ni temoin. Ne dit rien, et n'elimine pas.
#
# Tally sur le lot de controle (2 647 extensions Hoppe-stables) :
# 105 instables DEMONTRES (4,0 %), 1 722 avec temoin a J_max = 24, le
# reste indetermine.


def _sous_faisceaux(ext):
    """
    c1 de tous les sous-faisceaux propres de V lisibles sur la filtration.

    - toute sous-somme non vide de F1 est un sous-faisceau de F1 c V ;
    - pour toute sous-somme PROPRE S de F2, la preimage de S par V ->> F2
      est un sous-faisceau de V, de c1 = c1(F1) + c1(S).

    Le RANG n'intervient pas : mu(V) = 0, donc mu(W) < 0 <=> deg(W) < 0.
    """
    from itertools import combinations
    m = ext.m
    out = []
    for r in range(1, ext.rank_F1 + 1):
        for idx in combinations(range(ext.rank_F1), r):
            out.append(tuple(sum(ext.f1_charges[i][k] for i in idx)
                             for k in range(m)))
    c1f1 = ext.c1_F1
    for r in range(1, ext.rank_F2):        # S propre : |S| < rk(F2)
        for idx in combinations(range(ext.rank_F2), r):
            out.append(tuple(c1f1[k] + sum(ext.f2_charges[i][k] for i in idx)
                             for k in range(m)))
    return out


def degre(d_ijk, v, J):
    """deg_J(O(v)) = sum_ijk d_ijk v_i J_j J_k. Entier, exact."""
    v = np.asarray(v, dtype=np.int64)
    J = np.asarray(J, dtype=np.int64)
    return int(np.einsum('ijk,i,j,k->', np.asarray(d_ijk, dtype=np.int64),
                         v, J, J))


def _candidats_J(m, J_max=3, budget=200, seed=0):
    """
    Classes de Kahler testees, dans l'orthant positif.

    Renvoie (liste, exhaustive). `exhaustive` dit si la liste EPUISE
    [1, J_max]^m : un « aucun temoin » y est un enonce sur toute la boite,
    alors que sur un echantillon ce n'est qu'un sondage. Meme distinction
    que `ext_mode` pour le domaine des extensions (§5.11) -- ne pas la
    tracer reviendrait a confondre « pas de solution » et « pas cherche ».
    """
    from itertools import product
    if J_max ** m <= budget:
        return [tuple(j) for j in product(range(1, J_max + 1), repeat=m)], True
    base = [tuple([1] * m)]
    for i in range(m):
        for val in range(2, J_max + 1):
            J = [1] * m
            J[i] = val
            base.append(tuple(J))
    rng = np.random.RandomState(seed + m)
    while len(base) < budget:
        base.append(tuple(int(x) for x in rng.randint(1, J_max + 1, size=m)))
    return base[:budget], False


def certificat_instabilite(vs, coeff_max=3, taille_max=3):
    """
    PREUVE d'instabilite, ou None. Ne cherche pas : demontre.

    deg_J(v) = sum_i v_i * D_i(J)  avec  D_i(J) = sum_jk d_ijk J_j J_k.
    Les d_ijk d'une CICY dans un produit d'espaces projectifs sont positifs
    ou nuls -- verifie sur les 7 890 entrees de `cicylist.txt`, D_i(J) >= 0
    et D_i(1,..,1) > 0 avec. Donc D(J) est un point de l'orthant positif,
    quelle que soit la classe de Kahler.

    Si une combinaison a coefficients y_k >= 0 des c1 des sous-faisceaux
    verifie  sum_k y_k v_k >= 0  composante par composante, alors

        sum_k y_k * deg_J(v_k) = (sum_k y_k v_k) . D(J) >= 0

    donc au moins un deg_J(v_k) est >= 0 : le sous-faisceau correspondant
    ne peut pas etre de pente strictement negative. Et ceci pour TOUT J de
    l'orthant. C'est le sens facile du theoreme de transposition de
    Motzkin ; il suffit d'exhiber le y.

    RESERVE : la reciproque -- l'absence d'un tel y prouverait l'existence
    d'un p >= 0 avec v.p < 0 partout -- demanderait un solveur de
    programmation lineaire, et le depot ne depend que de numpy. On
    n'explore donc que les y entiers a petits coefficients : trouver un
    certificat DEMONTRE l'instabilite, ne pas en trouver ne demontre rien.

    Mesure de saturation sur le lot de controle : 105 certificats a
    taille <= 2, et pas un de plus en montant a taille <= 3 et
    coefficients <= 3. La granularite n'est donc pas le facteur limitant
    sur ce lot.
    """
    from itertools import combinations, product
    n = len(vs)
    m = len(vs[0])
    for r in range(1, min(taille_max, n) + 1):
        for idx in combinations(range(n), r):
            for co in product(range(1, coeff_max + 1), repeat=r):
                s = [sum(c * vs[i][k] for c, i in zip(co, idx))
                     for k in range(m)]
                if all(x >= 0 for x in s):
                    return {'y': {int(i): int(c) for i, c in zip(idx, co)},
                            'somme': s}
    return None


class ContextePente:
    """
    Classes de Kahler candidates et leurs D_i(J), pour UNE CICY.

    A construire une fois par CICY et a reutiliser : sans cela le test
    recalculerait D(J) = sum_jk d_ijk J_j J_k pour chaque extension, alors
    que D ne depend que de la geometrie. Avec la mise en cache, le cout
    par extension tombe a un produit matriciel (n_sous_faisceaux x m) par
    (m x n_J), soit quelques centaines de microsecondes.
    """

    __slots__ = ('D', 'cands', 'exhaustif', 'm')

    def __init__(self, d_ijk, m, J_max=24, budget=20000, seed=0):
        self.m = m
        self.cands, self.exhaustif = _candidats_J(m, J_max, budget, seed)
        d = np.asarray(d_ijk, dtype=np.int64)
        Js = np.asarray(self.cands, dtype=np.int64)
        # D[t, i] = sum_jk d_ijk J^t_j J^t_k
        self.D = np.einsum('ijk,tj,tk->ti', d, Js, Js)


def pente_extension(ext, d_ijk=None, J_max=24, budget=20000, seed=0,
                    coeff_max=3, taille_max=3, ctx=None):
    """
    Le fibre d'extension peut-il etre de pente stable pour une classe de
    Kahler de l'orthant positif ?

    Deux etapes, dans cet ordre, et surtout pas confondues :

    (a) CERTIFICAT. `certificat_instabilite` demontre, quand il aboutit,
        qu'aucune polarisation ne convient. Verdict `False`.

    (b) TEMOIN. Sinon on cherche un J rendant tous les degres strictement
        negatifs. Trouve : `True`, condition NECESSAIRE satisfaite pour ce
        J -- cela ne prouve pas la stabilite, les sous-faisceaux non
        lisibles sur la filtration n'etant pas testes. Non trouve :
        `None`, JAMAIS `False` : le tableau de saturation en tete de
        section montre que ce chiffre-la suit le budget de recherche.

    `J_exhaustif` dit si la grille des J epuise [1, J_max]^m. Un `None`
    n'a pas la meme portee selon ce drapeau -- et n'en a de toute facon
    aucune comme elimination.

    `ctx` : un `ContextePente` deja construit pour cette CICY. Fortement
    recommande en boucle de scan ; sinon il en est construit un, ce qui
    recalcule tous les D(J).
    """
    vs = _sous_faisceaux(ext)
    cert = certificat_instabilite(vs, coeff_max, taille_max)
    if cert is not None:
        return {'stable_possible': False,
                'etat': 'instable : certificat de Motzkin '
                        f'y = {cert["y"]}, somme = {cert["somme"]} >= 0',
                'temoin': None, 'certificat': cert, 'J_exhaustif': True}
    if ctx is None:
        ctx = ContextePente(d_ijk, ext.m, J_max, budget, seed)
    cands, exhaustif = ctx.cands, ctx.exhaustif
    # degres[k, t] = deg_{J_t}(v_k) = v_k . D[t]
    degres = np.asarray(vs, dtype=np.int64) @ ctx.D.T
    bons = np.flatnonzero((degres < 0).all(axis=0))
    if bons.size:
        J = cands[int(bons[0])]
        return {'stable_possible': True,
                'etat': f'temoin de polarisation J = {list(J)}',
                'temoin': list(J), 'certificat': None,
                'J_exhaustif': exhaustif}
    return {'stable_possible': None,
            'etat': ('indetermine : ni certificat, ni temoin sur toute la '
                     f'grille [1,{J_max}]^m' if exhaustif else
                     'indetermine : ni certificat, ni temoin sur '
                     f'{len(cands)} classes de Kahler tirees'),
            'temoin': None, 'certificat': None, 'J_exhaustif': exhaustif}


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
    p = 1..rk-1, et le critere de Hoppe est satisfait. La BORNE est
    suffisante et ne peut pas donner de faux positif.

    ATTENTION a ce que « satisfait » veut dire. Le critere de Hoppe est une
    equivalence avec la stabilite SEULEMENT si Pic(X) est de rang 1 ; sur
    une CICY a m > 1 il reste necessaire sans etre suffisant, parce qu'il
    ne voit pas la classe de Kahler. `stable: True` se lit donc « non
    elimine par Hoppe ». C'est `pente_extension` qui teste ce que Hoppe ne
    voit pas -- et il elimine 105 des 2 647 extensions declarees
    Hoppe-stables sur le scan de controle.

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


def cohomology_extension_ex(ext, ambient, config, chical=None):
    """
    H^i(V) par BORNES RIGOUREUSES, sans hypothese de rang maximal.

    ----------------------------------------------------------------------
    Pourquoi ne pas utiliser `compute_extension_cohomology`
    ----------------------------------------------------------------------
    Celle-ci pose les morphismes de liaison « de rang maximal pour une
    extension generique ». C'est exactement l'hypothese qui avait fausse
    `monads.compute_monad_cohomology` (§4.4) et `monad_wedge` v1 (§4.3) :
    n_gen restait juste -- il vaut |chi(V)| -- mais la REPARTITION entre h1
    et h2, donc n_anti, donc le classement, etait fausse (accord 0/11 avec
    la version rigoureuse). On ne refait pas l'erreur ici.

    ----------------------------------------------------------------------
    Ce qui est reellement connu
    ----------------------------------------------------------------------
    A n'appeler qu'apres que `hoppe_extension` a PROUVE la stabilite : on
    s'appuie sur h0(V) = h3(V) = 0, qui en decoule (une section globale
    donnerait un sous-faisceau O de meme pente ; h3(V) = h0(w^{rk-1}V) par
    det V = O et Serre, et hoppe_extension annule ce terme aussi).

    La suite exacte longue de 0 -> F1 -> V -> F2 -> 0 s'ecrit alors, avec
    r1 = rg( H1(F2) -> H2(F1) ) et r2 = rg( H2(F2) -> H3(F1) ) :

        h1(V) = h1(F1) + h1(F2) - r1
        h2(V) = h2(F1) + h2(F2) - r1 - r2
        h3(V) = h3(F1) + h3(F2) - r2      ( = 0 par stabilite)

    h0(F1) = h0(F2) = 0 est garanti par la borne p = 1 de hoppe_extension,
    donc le premier morphisme de liaison est nul et n'introduit aucune
    inconnue.

    r1 et r2 sont INCONNUS ; on ne renvoie donc que leurs consequences
    certaines. On verifie au passage que chi ne depend ni de r1 ni de r2 --
    il retombe sur chi(F1) + chi(F2), calcule par un tout autre chemin
    (Riemann-Roch, `chi_extension`).

    Renvoie None si un seul h^i de F1 ou F2 n'est pas certifie : sans cela
    les bornes porteraient sur des nombres inventes.
    """
    from cy_landscape.core.exact_cohomology import koszul_cohomology_ex

    # On n'exige la certification QUE des degres reellement utilises : 1 et
    # 2, dont dependent h1(V) et r1. Meme politique que le chemin des
    # monades, ou exiger les quatre degres rejetterait trois fois plus de
    # cas sans rien apporter au critere de selection. Les degres 0 et 3
    # servent, quand ils sont certifies, de CONTROLES de coherence.
    hF, cert = [], []
    for charges in (ext.f1_charges, ext.f2_charges):
        h = {i: 0 for i in range(4)}
        cd = {i: True for i in range(4)}
        for ch in charges:
            r = koszul_cohomology_ex(ambient, config, ch)
            if not (r['certified_by_degree'][1] and r['certified_by_degree'][2]):
                return None
            for i in range(4):
                h[i] += r[i]
                cd[i] = cd[i] and bool(r['certified_by_degree'][i])
        hF.append(h)
        cert.append(cd)
    h1, h2 = hF

    # chi vient de Riemann-Roch (arithmetique pure, toujours exact), et NON
    # de la somme alternee des h^i, qui ne le serait que si les quatre
    # degres etaient certifies. Quand ils le sont, les deux chemins sont
    # confrontes : c'est le seul controle croise disponible ici.
    if chical is None:
        from cy_landscape.core.chi_exact import make_calculator
        chical = make_calculator(ambient, config)
    chi = chi_extension(ext, chical)
    tous_certifies = all(cert[j][i] for j in (0, 1) for i in range(4))
    if tous_certifies:
        chi_alt = sum((-1) ** i * (h1[i] + h2[i]) for i in range(4))
        if chi_alt != chi:
            return None

    # h0(V) = h3(V) = 0 sont IMPOSES par la stabilite (deja prouvee).
    # Quand h0 est certifie, il doit valoir 0 : sinon H0(F1) s'injecterait
    # dans H0(V) = 0. Une violation signale un appel hors contrat.
    if cert[0][0] and h1[0]:
        return None
    if cert[1][0] and h2[0]:
        return None

    r1_max = min(h2[1], h1[2])          # rg( H1(F2) -> H2(F1) )
    if cert[0][3] and cert[1][2] and cert[1][3]:
        # h3(V) = h3(F1) + h3(F2) - r2 doit pouvoir s'annuler.
        if h1[3] + h2[3] - min(h2[2], h1[3]) > 0:
            return None

    lo1 = max(0, h1[1] + h2[1] - r1_max, -chi)
    hi1 = h1[1] + h2[1]
    if lo1 > hi1:
        return None
    bounds = {0: (0, 0), 1: (lo1, hi1), 2: (lo1 + chi, hi1 + chi), 3: (0, 0)}
    out = {'chi': int(chi), 'bounds': bounds,
           'determine': {i: bounds[i][0] == bounds[i][1] for i in range(4)},
           'chi_recoupe': bool(tous_certifies),
           'hF1': dict(h1), 'hF2': dict(h2)}
    for i in range(4):
        out[i] = int(bounds[i][0]) if bounds[i][0] == bounds[i][1] else None
    return out
