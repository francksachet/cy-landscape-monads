"""
positive_monads.py -- Monades positives sur les CICYs

Une monade 0 -> V -> B -> C -> 0 est POSITIVE si pour tout (j,i),
les composantes de (c_j - b_i) sont >= 0 et au moins une > 0.

Proprietes :
  - La map f: B -> C est automatiquement surjective (sections abondantes)
  - V est un fibre (pas un faisceau)
  - Theoreme (Flenner, 1984) : sur P^n, les monades positives
    sont stables pour n >= 2. Sur les CICYs, la positivite
    est une forte indication de stabilite.

Ref: Okonek, Schneider, Spindler "Vector Bundles on Complex Projective Spaces"
     Anderson et al. "Monad Bundles in Heterotic String Compactifications" (2008)

--------------------------------------------------------------------------
VERSION 2 du generateur (GENERATOR_VERSION = 2)
--------------------------------------------------------------------------
Trois changements par rapport a la v1 :

1. ELAGAGE EXACT (aucune perte mathematique)
   Pour r_C = 1, `_find_positive_B` impose b_i[k] <= c[k] - 1 pour toute
   composante k avec c[k] > 0, et sum_i b_i[k] = c[k]. Ce systeme est
   infaisable des que c[k] = 1 : tous les b_i[k] valent alors 0 et leur
   somme ne peut pas valoir 1. Tout vecteur c ayant au moins une
   composante egale a 1 est donc rejete d'emblee, sans aucun tirage.
   Mesure sur m = 6, max_charge = 4 : 15 624 vecteurs c -> 4 095 retenus
   (facteur 3,8) pour un espace de recherche STRICTEMENT identique.

2. ECHANTILLONNAGE HYBRIDE (opt-in, desactive par defaut)
   Au-dela d'un seuil sur la charge totale sum(c), on ne parcourt plus
   tous les vecteurs c : on en tire un sous-ensemble aleatoire uniforme
   plafonne a `sampling_budget` par valeur de total (tirage sans remise
   dans TOUT l'espace restant, pas une troncature des premiers). En
   dessous du seuil, l'enumeration reste exhaustive. Desactive par
   defaut (sampling_threshold=None).

3. ALEA DETERMINISTE PAR VECTEUR c
   En v1, `_find_positive_B` puisait dans un RNG partage : l'ordre
   d'enumeration influencait les tirages. En v2 chaque vecteur c a son
   propre RNG derive de (seed, c, r_B). Le resultat pour un c donne ne
   depend plus de ce qui a ete enumere avant : l'elagage et
   l'echantillonnage sont sans effet de bord, et le scan est
   reproductible independamment de l'ordre et du nombre de workers.

Tracabilite : GENERATOR_VERSION est enregistre dans chaque ligne de
results.jsonl (champ "gen"), afin de distinguer a posteriori les CICYs
scannees avec la v1 de celles scannees avec la v2.
"""
import numpy as np
from itertools import combinations, product as iprod
from typing import List, Optional

from itertools import product as _iproduct

from cy_landscape.core.monads import MonadBundle

GENERATOR_VERSION = 2

_UINT32 = 2 ** 31 - 1
# Au-dela de cette taille, on n'enumere plus l'espace des c en memoire.
_MAX_ENUMERABLE = 2_000_000


def _mix(seed, c_tuple, r_B):
    """
    Graine deterministe derivee de (seed, c, r_B).

    N'utilise PAS hash() : le hachage des tuples Python est randomise par
    processus (PYTHONHASHSEED), ce qui rendrait les tirages differents
    d'un worker a l'autre et non reproductibles entre deux lancements.
    """
    h = (seed * 1_000_003 + r_B * 7919 + 2_166_136_261) & 0xFFFFFFFF
    for v in c_tuple:
        h = ((h ^ (int(v) + 0x9E3779B9)) * 16_777_619) & 0xFFFFFFFF
    return h % _UINT32


def is_positive_monad(monad):
    """
    Verifie si la monade est positive :
    pour tout j, i : c_j - b_i >= 0 composante par composante,
    avec au moins une composante > 0.
    """
    for cj in monad.c_charges:
        for bi in monad.b_charges:
            diff = [cj[k] - bi[k] for k in range(monad.m)]
            if any(d < 0 for d in diff):
                return False
            if all(d == 0 for d in diff):
                return False  # Doit etre strictement positif
    return True


def _c_is_feasible(c_ch, r_B):
    """
    Condition NECESSAIRE ET SUFFISANTE, composante par composante, pour
    qu'il existe r_B vecteurs b_i >= 0 avec b_i[k] <= c[k]-1 (si c[k]>0)
    et sum_i b_i[k] = c[k] -- c'est exactement la contrainte imposee par
    `_find_positive_B`.

    Pour une composante k :
      - c[k] = 0 : tous les b_i[k] = 0, faisable.
      - c[k] > 0 : il faut c[k] <= r_B * (c[k] - 1), soit
                   c[k] * (r_B - 1) >= r_B. Avec r_B >= 2 cela revient a
                   c[k] >= 2.
    Aucune approximation : un c rejete ici ne pouvait produire aucun B.
    """
    if all(x == 0 for x in c_ch):
        return False
    for x in c_ch:
        if x == 0:
            continue
        if x * (r_B - 1) < r_B:
            return False
    return True


def _c_candidates(m, max_charge, r_B, sampling_threshold=None,
                  sampling_budget=400, seed=0):
    """
    Liste des vecteurs c candidats pour r_C = 1, par charge totale croissante.

    - Tout c infaisable (cf. `_c_is_feasible`) est ecarte : elagage exact.
    - `sampling_threshold` None -> enumeration exhaustive.
    - Sinon, pour les totaux strictement superieurs au seuil, tirage
      uniforme sans remise d'au plus `sampling_budget` vecteurs.
    """
    space = (max_charge + 1) ** m
    rng = np.random.RandomState(_mix(seed, (m, max_charge), r_B))

    by_total = {}
    if space <= _MAX_ENUMERABLE:
        for c in iprod(range(max_charge + 1), repeat=m):
            if not _c_is_feasible(c, r_B):
                continue
            by_total.setdefault(sum(c), []).append(c)
    else:
        # Espace trop grand pour etre materialise : tirage direct.
        target = max(sampling_budget, 1) * (max_charge * m + 1)
        seen = set()
        for _ in range(target * 20):
            if len(seen) >= target:
                break
            c = tuple(int(v) for v in rng.randint(0, max_charge + 1, size=m))
            if c in seen or not _c_is_feasible(c, r_B):
                continue
            seen.add(c)
            by_total.setdefault(sum(c), []).append(c)

    out = []
    for total in sorted(by_total):
        bucket = by_total[total]
        if (sampling_threshold is None or total <= sampling_threshold
                or len(bucket) <= sampling_budget):
            out.extend(bucket)
        else:
            idx = rng.choice(len(bucket), size=sampling_budget, replace=False)
            out.extend(bucket[i] for i in sorted(idx))
    return out


def generate_positive_monads(m, rank_V, max_charge=4, n_systematic=500, rng=None,
                             sampling_threshold=None, sampling_budget=400,
                             seed=42, exhaustif_max=0, stats=None):
    """
    Genere des monades positives c1(V) = 0.

    Strategie :
    1) Fixer C (charges positives)
    2) Choisir B avec b_i <= c_j pour tout j, composante par composante
    3) Verifier c1(V) = c1(B) - c1(C) = 0

    Parametres ajoutes en v2 :
      sampling_threshold : None -> enumeration exhaustive des c (defaut).
                           Entier t -> exhaustif pour sum(c) <= t, puis
                           echantillonnage aleatoire plafonne au-dela.
      sampling_budget    : nombre max de vecteurs c tires par valeur de
                           charge totale, au-dela du seuil.
      seed               : graine des RNG derives (deterministe par c).

    Parametre ajoute en v3 :
      exhaustif_max : 0 (defaut) -> comportement inchange, 50 tirages
                      aleatoires de B par vecteur c. Entier P > 0 -> pour
                      tout c dont le nombre de B admissibles vaut au plus P,
                      on les ENUMERE TOUS ; au-dela seulement, on echantillonne.
                      `stats`, si un dict est fourni, recoit le decompte des
                      c traites de chaque facon -- indispensable pour savoir
                      sur quel domaine un resultat d'absence est demontre et
                      sur quel domaine il n'est qu'un sondage.
    """
    if stats is not None:
        stats.setdefault('c_exhaustifs', 0)
        stats.setdefault('c_echantillonnes', 0)
        stats.setdefault('B_exhaustifs', 0)
    if rng is None:
        rng = np.random.RandomState(42)

    monads = []

    for r_C in [1, 2]:
        r_B = rank_V + r_C

        if r_C == 1:
            # C = O(c) avec c > 0
            # B = O(b_1), ..., O(b_{rB}) avec b_i < c (composante par composante)
            # et sum(b_i) = c (pour c1 = 0)
            candidates = _c_candidates(
                m, max_charge, r_B,
                sampling_threshold=sampling_threshold,
                sampling_budget=sampling_budget,
                seed=seed)

            for c_tuple in candidates:
                c_ch = list(c_tuple)
                c_rng = np.random.RandomState(_mix(seed, c_tuple, r_B))

                b_list = None
                if exhaustif_max > 0:
                    b_list = enumerer_positive_B(c_ch, r_B,
                                                 plafond=exhaustif_max)
                if b_list is None:
                    b_list = _find_positive_B(c_ch, r_B, c_rng,
                                              max_attempts=50)
                    if stats is not None:
                        stats['c_echantillonnes'] += 1
                elif stats is not None:
                    stats['c_exhaustifs'] += 1
                for b_charges in b_list:
                    if stats is not None:
                        stats['B_exhaustifs'] += 1
                    monad = MonadBundle(b_charges, [c_ch])
                    if monad.c1_vanishes and is_positive_monad(monad):
                        monads.append(monad)

        elif r_C == 2:
            # C = O(c1) + O(c2), B = r_B line bundles
            #
            # RNG PROPRE A LA BRANCHE. Elle puisait dans le `rng` partage,
            # donc ses tirages dependaient de tout ce qui avait ete consomme
            # avant -- c'est-a-dire de la branche r_C = 1. La v2 avait corrige
            # cet effet de bord pour r_C = 1 (un RNG derive par vecteur c) mais
            # avait laisse celle-ci en l'etat. Consequence mesuree : activer
            # `exhaustif_max`, qui ne touche pourtant QUE la branche r_C = 1,
            # faisait disparaitre 20 candidats de rank_C = 2 sur 42 -- le
            # generateur n'etait pas monotone en son propre parametre.
            rng2 = np.random.RandomState(_mix(seed, (m, max_charge, 2), r_B))
            for _ in range(n_systematic):
                c1 = [int(rng2.randint(0, max_charge + 1)) for _ in range(m)]
                c2 = [int(rng2.randint(0, max_charge + 1)) for _ in range(m)]
                if all(x == 0 for x in c1) or all(x == 0 for x in c2):
                    continue

                c_charges = [c1, c2]
                c1_total = [c1[k] + c2[k] for k in range(m)]

                # B : rB vecteurs avec b_i <= min(c1, c2) et sum = c1_total
                c_min = [min(c1[k], c2[k]) for k in range(m)]
                b_charges = _find_positive_B_multi(c_min, c1_total, r_B, rng2)
                if b_charges is None:
                    continue

                monad = MonadBundle(b_charges, c_charges)
                if monad.c1_vanishes and is_positive_monad(monad):
                    monads.append(monad)

    return monads


def _partitions_bounded(total, m, max_val):
    """
    Genere les partitions de 'total' en m parts avec 0 <= part <= max_val.

    Conservee pour compatibilite et tests : le generateur principal passe
    desormais par `_c_candidates`, qui enumere directement le produit
    cartesien avec elagage de faisabilite.
    """
    if m == 1:
        if 0 <= total <= max_val:
            yield (total,)
        return
    for v in range(min(total, max_val) + 1):
        for rest in _partitions_bounded(total - v, m - 1, max_val):
            yield (v,) + rest


def _compositions_bornees(total, parts, borne):
    """
    Listes de `parts` entiers de [0, borne] de somme `total`. Generateur.
    """
    if parts == 1:
        if 0 <= total <= borne:
            yield (total,)
        return
    lo = max(0, total - borne * (parts - 1))
    hi = min(borne, total)
    for x in range(lo, hi + 1):
        for reste in _compositions_bornees(total - x, parts - 1, borne):
            yield (x,) + reste


def _compte_compositions(total, parts, borne):
    """Nombre de telles listes, par programmation dynamique -- sans enumerer."""
    if borne < 0:
        return 0
    dp = [1] + [0] * total
    for _ in range(parts):
        ndp = [0] * (total + 1)
        cum = 0
        for s in range(total + 1):
            cum += dp[s]
            if s - borne - 1 >= 0:
                cum -= dp[s - borne - 1]
            ndp[s] = cum
        dp = ndp
    return dp[total]


def compte_B(c_ch, r_B):
    """
    Nombre de tuples ORDONNES (b_1, ..., b_rB) admissibles, SANS les enumerer.

    Contraintes, identiques a celles de `_find_positive_B` :
      - b_i[k] = 0 si c[k] = 0 ;
      - 0 <= b_i[k] <= c[k] - 1 sinon ;
      - somme_i b_i[k] = c[k].

    Le produit sur les composantes suffit : les contraintes ne couplent pas
    les composantes entre elles. Sert de garde AVANT enumeration, et de
    reference independante APRES (le nombre enumere doit coincider).
    """
    # c = 0 : le seul B possible est nul, et c_j - b_i vaut alors 0 partout,
    # ce que `is_positive_monad` rejette. Zero, donc, et non un. Sans cette
    # ligne l'enumeration produisait une monade degeneree que
    # `_c_is_feasible` ecarte pourtant depuis toujours -- les deux chemins
    # doivent s'accorder.
    if all(x == 0 for x in c_ch):
        return 0
    total = 1
    for ck in c_ch:
        if ck == 0:
            continue
        n = _compte_compositions(ck, r_B, ck - 1)
        if n == 0:
            return 0
        total *= n
    return total


def enumerer_positive_B(c_ch, r_B, plafond=200000):
    """
    TOUS les B admissibles pour un c donne, a permutation pres des b_i.

    Renvoie None si le nombre de tuples ordonnes depasse `plafond` -- charge
    a l'appelant de retomber sur l'echantillonnage.

    ----------------------------------------------------------------------
    Pourquoi cette fonction existe
    ----------------------------------------------------------------------
    `_find_positive_B` tire 50 B au hasard parmi N possibles. Mesure sur les
    vecteurs c reellement rencontres : N va de 2.10^4 a 3.10^5, et 50 graines
    successives ne decouvrent que ~120 B distincts, soit 0,04 % -- sans le
    moindre signe de saturation, la croissance restant lineaire. Un resultat
    d'ABSENCE obtenu ainsi ne porte donc sur rien.

    Enumerer exhaustivement change la nature de l'enonce : on passe de
    « aucun survivant parmi ce qu'on a tire » a « aucun survivant, sur tout
    le domaine ou N <= plafond ». C'est le seul moyen d'obtenir un resultat
    negatif defendable.

    La deduplication par multiensemble est indispensable : deux B qui ne
    different que par l'ordre des b_i donnent le MEME fibre, et les compter
    deux fois gonflerait artificiellement la couverture d'un facteur r_B!
    """
    n_ordonnes = compte_B(c_ch, r_B)
    if n_ordonnes == 0 or n_ordonnes > plafond:
        return None

    m = len(c_ch)
    par_composante = []
    for k in range(m):
        if c_ch[k] == 0:
            par_composante.append([(0,) * r_B])
        else:
            par_composante.append(
                list(_compositions_bornees(c_ch[k], r_B, c_ch[k] - 1)))

    return _enumerer(par_composante, m, r_B)


def _enumerer(par_composante, m, r_B):
    """
    GENERATEUR. Ne materialise que l'ensemble de deduplication.

    La premiere version renvoyait une liste : elle construisait, en plus des
    cles de deduplication, une liste complete de listes de listes. Avec un
    plafond a 500 000 et sept workers en parallele, cela a produit un
    MemoryError apres 153 minutes de scan (scan_su5). Le plafond bornait le
    NOMBRE de B, pas la memoire -- deux choses differentes des lors qu'on
    conserve chaque B.
    """
    vus = set()
    for combi in _iproduct(*par_composante):
        cle = tuple(sorted(tuple(combi[k][i] for k in range(m))
                           for i in range(r_B)))
        if cle in vus:
            continue
        vus.add(cle)
        yield [list(b) for b in cle]


def _find_positive_B(c_ch, r_B, rng, max_attempts=50):
    """
    Trouve des ensembles de r_B vecteurs b_i tels que :
    - 0 <= b_i[k] < c_ch[k] pour tout k (ou b_i[k] = 0 si c_ch[k] = 0)
    - sum(b_i) = c_ch
    """
    m = len(c_ch)
    results = []

    for _ in range(max_attempts):
        b_charges = []
        remaining = list(c_ch)

        for r in range(r_B - 1):
            b = [0] * m
            for k in range(m):
                if remaining[k] > 0:
                    # max pour ce b[k] : min(remaining[k], c_ch[k] - 1)
                    upper = min(remaining[k], c_ch[k] - 1)
                    if upper > 0:
                        b[k] = int(rng.randint(0, upper + 1))
                    remaining[k] -= b[k]
            b_charges.append(b)

        # Dernier bundle = remaining
        last = remaining
        # Verifier positivite : last < c
        if all(last[k] < c_ch[k] or c_ch[k] == 0 for k in range(m)):
            if any(last[k] > 0 for k in range(m)) or all(c_ch[k] == 0 for k in range(m)):
                b_charges.append(last)
                results.append(b_charges)
                if len(results) >= 3:
                    break

    return results


def _find_positive_B_multi(c_min, c1_total, r_B, rng, max_attempts=30):
    """Trouve B pour r_C = 2 : b_i <= c_min et sum = c1_total."""
    m = len(c_min)
    for _ in range(max_attempts):
        b_charges = []
        remaining = list(c1_total)
        ok = True
        for r in range(r_B - 1):
            b = [0] * m
            for k in range(m):
                upper = min(remaining[k], c_min[k])
                if upper > 0:
                    b[k] = int(rng.randint(0, upper + 1))
                elif upper == 0:
                    b[k] = 0
                else:
                    ok = False; break
                remaining[k] -= b[k]
            if not ok: break
            b_charges.append(b)
        if not ok: continue
        last = remaining
        if all(last[k] <= c_min[k] for k in range(m)):
            b_charges.append(last)
            return b_charges
    return None
