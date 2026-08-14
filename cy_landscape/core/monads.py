"""
monads.py -- Fibrés monades sur les CICYs

Construction : 0 -> V -> B -> C -> 0
  B = ⊕ O(b_i)  somme de line bundles (rang r_B)
  C = ⊕ O(c_j)  somme de line bundles (rang r_C)
  V = ker(f)     fibre monade (rang r_V = r_B - r_C)

La map f: B -> C est une matrice de sections f_{ji} ∈ H⁰(O(c_j - b_i)).
Pour que V soit un fibre (pas un faisceau), f doit etre surjective.

Cohomologie de V : suite exacte longue
  0 -> H⁰(V) -> H⁰(B) -> H⁰(C) -> H¹(V) -> H¹(B) -> H¹(C) 
    -> H²(V) -> H²(B) -> H²(C) -> H³(V) -> H³(B) -> H³(C) -> 0

Ref: Anderson, Gray, Lukas, Palti — "Two Hundred Heterotic Standard Models on Smooth CY Threefolds"
"""
import numpy as np
from itertools import product
from dataclasses import dataclass, field
from typing import List, Optional

from cy_landscape.core.exact_cohomology import bundle_cohomology_exact


_UINT32 = 2 ** 31 - 1


def _mix(seed, c_tuple, r_B):
    """
    Graine deterministe derivee de (seed, c, r_B).

    N'utilise PAS hash() : le hachage des tuples Python est randomise par
    processus (PYTHONHASHSEED), ce qui rendrait les tirages differents
    d'un worker a l'autre et non reproductibles entre deux lancements.

    (Definie ici plutot que dans positive_monads.py, qui importe ce
    module : `generate_monads` en a besoin, et l'inverse serait circulaire.
    `positive_monads` la reexporte, la signature est inchangee.)
    """
    h = (seed * 1_000_003 + r_B * 7919 + 2_166_136_261) & 0xFFFFFFFF
    for v in c_tuple:
        h = ((h ^ (int(v) + 0x9E3779B9)) * 16_777_619) & 0xFFFFFFFF
    return h % _UINT32


@dataclass
class MonadBundle:
    """Un fibre monade defini par ses charges B et C."""
    b_charges: List[List[int]]   # r_B vecteurs de charges, chacun de dimension m
    c_charges: List[List[int]]   # r_C vecteurs de charges, chacun de dimension m

    @property
    def rank_B(self): return len(self.b_charges)

    @property
    def rank_C(self): return len(self.c_charges)

    @property
    def rank_V(self): return self.rank_B - self.rank_C

    @property
    def m(self): return len(self.b_charges[0])

    @property
    def c1_B(self):
        return [sum(b[i] for b in self.b_charges) for i in range(self.m)]

    @property
    def c1_C(self):
        return [sum(c[i] for c in self.c_charges) for i in range(self.m)]

    @property
    def c1_V(self):
        """c1(V) = c1(B) - c1(C). Doit etre 0 pour SU(n)."""
        return [self.c1_B[i] - self.c1_C[i] for i in range(self.m)]

    @property
    def c1_vanishes(self):
        return all(x == 0 for x in self.c1_V)

    def map_degrees(self):
        """
        Retourne la matrice des degres c_j - b_i (r_C x r_B).
        f_{ji} ∈ H⁰(O(c_j - b_i)) doit etre non-trivial.
        """
        degs = []
        for cj in self.c_charges:
            row = []
            for bi in self.b_charges:
                row.append(tuple(cj[k] - bi[k] for k in range(self.m)))
            degs.append(row)
        return degs


def check_map_exists(monad, ambient_dims, config_matrix):
    """
    Verifie que la map f: B -> C peut exister.

    Pour chaque entree (j,i) de la matrice, on a besoin de
    H⁰(O_X(c_j - b_i)) > 0 pour au moins assez d'entrees.

    Condition necessaire (faible) : pour chaque j, il existe au moins
    un i tel que H⁰(O(c_j - b_i)) > 0.

    Condition suffisante pour la surjectivite generique (forte) :
    pour chaque j, il faut que la somme des H⁰ sur les i soit >= 1.
    En pratique, c_j - b_i >= 0 composante par composante suffit
    pour avoir des sections.
    """
    degs = monad.map_degrees()
    h0_matrix = []

    for j, row in enumerate(degs):
        h0_row = []
        for i, deg in enumerate(row):
            # Condition simple : toutes les composantes >= 0
            # => H⁰(O(deg)) > 0 sur l'espace ambiant
            all_nonneg = all(d >= 0 for d in deg)
            if all_nonneg and any(d > 0 for d in deg):
                h0_row.append(True)
            elif all(d == 0 for d in deg):
                h0_row.append(True)  # O(0) a une section constante
            else:
                # Calcul exact via cohomologie
                try:
                    cohom_full = bundle_cohomology_exact(
                        ambient_dims, config_matrix, [list(deg)])
                    cohom = cohom_full['V'] if isinstance(cohom_full, dict) and 'V' in cohom_full else cohom_full
                    h0_row.append(cohom.get(0, 0) > 0)
                except Exception:
                    h0_row.append(False)
        h0_matrix.append(h0_row)

    # Condition : chaque ligne (chaque c_j) doit avoir au moins un H⁰ > 0
    for j in range(monad.rank_C):
        if not any(h0_matrix[j]):
            return False, h0_matrix

    return True, h0_matrix


def compute_monad_cohomology(monad, ambient_dims, config_matrix):
    """
    Calcule H^i(V) via la suite exacte longue de 0 -> V -> B -> C -> 0.

    H^i(B) = ⊕ H^i(O(b_k))
    H^i(C) = ⊕ H^i(O(c_k))

    Suite exacte longue :
    0 -> H⁰(V) -> H⁰(B) -f-> H⁰(C) -> H¹(V) -> H¹(B) -> H¹(C)
      -> H²(V) -> H²(B) -> H²(C) -> H³(V) -> H³(B) -> H³(C) -> 0

    En pratique, pour un monad generique, f est surjective sur H⁰,
    donc H⁰(V) = ker(f|_{H⁰}) et la suite se simplifie.

    Approximation standard (monad generique) :
      h⁰(V) = max(0, h⁰(B) - h⁰(C))   [si f est injective sur H⁰]
      h¹(V) = h¹(B) - h¹(C) + max(0, h⁰(C) - h⁰(B))
      h²(V) = h²(B) - h²(C) + max(0, h³(C) - h³(B))  [par Serre]
      h³(V) = max(0, h³(B) - h³(C))   [si f* est injective]

    Mais on utilise la formule exacte via chi(V) = chi(B) - chi(C)
    et les bornes de la suite exacte longue.
    """
    # Calculer H^i pour chaque line bundle dans B et C
    hB = {0: 0, 1: 0, 2: 0, 3: 0}
    for b in monad.b_charges:
        try:
            h_full = bundle_cohomology_exact(ambient_dims, config_matrix, [b])
            h = h_full['V'] if isinstance(h_full, dict) and 'V' in h_full else h_full
            for i in range(4):
                hB[i] += h.get(i, 0)
        except Exception:
            return None

    hC = {0: 0, 1: 0, 2: 0, 3: 0}
    for c in monad.c_charges:
        try:
            h_full = bundle_cohomology_exact(ambient_dims, config_matrix, [c])
            h = h_full['V'] if isinstance(h_full, dict) and 'V' in h_full else h_full
            for i in range(4):
                hC[i] += h.get(i, 0)
        except Exception:
            return None

    # Suite exacte longue (approximation monad generique)
    # La map f: H⁰(B) -> H⁰(C) est generiquement de rang maximal
    # La map f*: H³(B) -> H³(C) aussi (par Serre duality)

    # Rang de f sur H⁰ : min(h⁰(B), h⁰(C))
    rank_f0 = min(hB[0], hC[0])
    h0V = hB[0] - rank_f0   # ker(f) sur H⁰
    coker_f0 = hC[0] - rank_f0   # coker(f) sur H⁰

    # Rang de f sur H³ (par dualite de Serre)
    rank_f3 = min(hB[3], hC[3])
    h3V = hB[3] - rank_f3
    coker_f3 = hC[3] - rank_f3  # unused but for completeness

    # chi(V) = chi(B) - chi(C) donne une contrainte
    chiB = sum((-1)**i * hB[i] for i in range(4))
    chiC = sum((-1)**i * hC[i] for i in range(4))
    chiV = chiB - chiC

    # h¹(V) - h²(V) = chiV - h⁰(V) + h³(V)
    # De la suite exacte :
    # h¹(V) >= max(0, h¹(B) - h¹(C) + coker_f0)
    # On utilise l'approximation generique :
    h1V = hB[1] - hC[1] + coker_f0
    if h1V < 0:
        h1V = 0  # borne inferieure

    # h²(V) par la contrainte chi
    h2V = h1V - chiV + h0V - h3V
    if h2V < 0:
        # Ajuster h1V
        h2V = 0
        h1V = chiV - h0V + h3V

    return {0: max(0, h0V), 1: max(0, h1V), 2: max(0, h2V), 3: max(0, h3V),
            'hB': dict(hB), 'hC': dict(hC), 'chi_V': chiV, 'exact': False}


def check_monad_stability(monad, ambient_dims, config_matrix, intersection_numbers):
    """
    Verifie la semi-stabilite du fibre monade.

    Pour un monad 0 -> V -> B -> C -> 0 :
    Critere de Hoppe : V est stable si H⁰(Λ^p V(-1)) = 0
    pour tout 1 <= p <= rk(V)-1 et tout diviseur ample.

    En pratique, on verifie une condition necessaire :
    tout sous-faisceau L ↪ V de rang 1 doit avoir mu(L) <= mu(V) = 0.

    Pour un monad, L ↪ V ↪ B, donc L ↪ B = ⊕ O(b_i).
    Cela signifie L = O(l) avec l <= b_i (composante par composante)
    pour au moins un i. Puis on verifie que la composition
    L -> B -> C est nulle, ie L ∈ ker(f).

    Condition suffisante simplifiee (Anderson et al.) :
    V est stable si pour tout l tel que :
      1) l <= b_i pour un i
      2) c_j - l n'est PAS effectif pour tout j
    on a mu(O(l)) <= 0.

    Et si aucun tel l avec mu(O(l)) > 0 n'existe, V est stable.
    """
    d = intersection_numbers
    m = monad.m
    rk = monad.rank_V

    if rk <= 0:
        return {"semi_stable": False, "reason": "rank <= 0"}

    # mu(V) = 0 car c1(V) = 0

    # Generer les sous-line-bundles potentiels
    # L = O(l) avec l <= b_i pour au moins un i
    # On verifie mu(L) = deg(L) / 1

    # Calculer le degre d'un line bundle O(l) via intersection
    # deg(L) = sum_{ijk} l_i * d_{ijk} * J_j * J_k (pour une polarisation J)
    # On utilise J = (1,...,1) comme polarisation test

    def degree(l):
        """Degre de O(l) pour la polarisation J=(1,...,1)."""
        deg = 0.0
        shape = d.shape
        for idx in np.ndindex(*shape):
            val = d[idx]
            if val != 0:
                deg += val * l[idx[0]]
        return deg

    # Bornes inferieures depuis les b_i
    b_min = [min(b[k] for b in monad.b_charges) for k in range(m)]
    b_max = [max(b[k] for b in monad.b_charges) for k in range(m)]

    destabilizing = []

    # Scanner les l possibles dans [b_min, b_max]
    ranges = [range(b_min[k], b_max[k] + 1) for k in range(m)]
    for l in product(*ranges):
        l = list(l)

        # Condition 1 : l <= b_i pour au moins un i
        sub_of_B = any(all(l[k] <= b[k] for k in range(m))
                       for b in monad.b_charges)
        if not sub_of_B:
            continue

        # Condition 2 : la composition L -> B -> C doit pouvoir etre nulle
        # ie pour tout j, c_j - l ne doit PAS donner une section non-nulle
        # qui forcerait la composition a etre non-nulle.
        # Condition suffisante pour L ⊂ ker(f) : 
        # c_j - l a toutes ses composantes < 0 pour tout j
        # (pas de section O(c_j - l) => la map L -> C_j est forcement 0)
        in_kernel = True
        for cj in monad.c_charges:
            diff = [cj[k] - l[k] for k in range(m)]
            # Si toutes >= 0, il y a une section, la composition pourrait etre non-nulle
            if all(d >= 0 for d in diff):
                in_kernel = False
                break

        if not in_kernel:
            continue

        # Ce l est un sous-line-bundle potentiel de V
        deg_l = degree(l)
        if deg_l > 0:
            destabilizing.append({'l': l, 'degree': deg_l})

    if destabilizing:
        return {
            "semi_stable": False,
            "reason": f"{len(destabilizing)} sous-faisceaux destabilisants",
            "examples": destabilizing[:3]
        }

    return {"semi_stable": True, "reason": "Aucun sous-faisceau destabilisant trouve"}


# ======================================================================
# La famille « sommes de vecteurs unite » -- enumeree, plus tiree
# ======================================================================
#
# CE QUE CE BLOC CORRIGE
# ----------------------
# L'ancien `generate_monads` produisait ses monades « anti-symetriques »
# par TIRAGE : `for _ in range(min(m*3, 10))`, soit DIX tirages, pris sur
# le RNG PARTAGE avec le generateur positif. Deux consequences.
#
#   1. Couverture. La famille visee -- des b_i de la forme e_a, ou
#      e_a +/- e_b -- compte 101 multiensembles valides rien que pour
#      m = 5, r_B = 5 dans sa strate pure, et 1 911 en autorisant un
#      vecteur perturbe. Dix tirages en voyaient au mieux dix.
#
#   2. Reproductibilite. Le RNG etant partage, le nombre de tirages
#      consommes en amont par `generate_positive_monads` decidait de la
#      valeur des dix tirages. Toute modification du generateur positif
#      -- y compris une correction sans rapport, comme celle du §5.11 --
#      redistribuait la loterie.
#
# C'est ce mecanisme qui a fait DISPARAITRE les deux candidats phares du
# projet (#6890 et #6947) entre deux scans, sans qu'aucun filtre ne les
# elimine : ils n'ont simplement plus ete engendres. Or tous deux sont
# des sommes de vecteurs unite PURS :
#
#      #6890 : B = O(e_1) + O(e_2)^3 + O(e_3)     -> multiensemble {1,2,2,2,3}
#      #6947 : B = O(e_0)^3 + O(e_1) + O(e_3)     -> multiensemble {0,0,0,1,3}
#      #6715 : B = O(e_3)^3 + O(e_0) + O(e_0+e_2) -> UN vecteur perturbe
#
# Les enumerer les rend DEMONTRES et non plus chanceux.
#
# CE QUI RESTE TIRE, ET POURQUOI
# ------------------------------
# La famille complete (tout b_i de la forme e_a +/- e_b) n'est pas
# enumerable : pour m = 5, r_B = 6, elle depasse deja 2 000 000 de
# multiensembles a |c1(B)|_inf <= 3. On enumere donc par STRATES, selon
# le nombre de vecteurs perturbes :
#
#      k = 0 : tous les b_i sont des vecteurs unite   (toujours enumeree)
#      k = 1 : au plus un b_i perturbe                (si sous le plafond)
#      k >= 2 : jamais enumere -- echantillonne, et DECLARE comme tel.
#
# Conformement a la regle des filtres (§8), `stats` recoit pour chaque
# strate le nombre enumere, le nombre total, et le mode. Un resultat
# d'absence sur cette branche n'est interpretable qu'avec ces trois
# nombres sous les yeux.

def _vecteurs_unite(m):
    """Les m vecteurs e_a."""
    return [tuple(1 if k == a else 0 for k in range(m)) for a in range(m)]


def _vecteurs_perturbes(m):
    """Les e_a + eps*e_b, a != b, eps = +/-1. Deduplique et trie : l'ordre
    ne doit dependre ni de PYTHONHASHSEED ni de l'ordre d'insertion."""
    out = set()
    for a in range(m):
        for b in range(m):
            if b == a:
                continue
            for eps in (-1, 1):
                q = [0] * m
                q[a] = 1
                q[b] += eps
                out.add(tuple(q))
    return sorted(out)


def _multisets_unite(m, taille, borne_mult):
    """
    Enumere les multiensembles de `taille` vecteurs unite dont chaque
    multiplicite est <= `borne_mult`.

    Rend directement les multiplicites (un vecteur de m entiers >= 0 de
    somme `taille`), ce qui evite de construire puis rejeter : la borne
    est appliquee pendant la descente, pas apres.
    """
    mult = [0] * m
    reste = taille

    def descente(i, reste):
        if i == m - 1:
            if reste <= borne_mult:
                mult[i] = reste
                yield tuple(mult)
                mult[i] = 0
            return
        haut = min(borne_mult, reste)
        for v in range(haut + 1):
            mult[i] = v
            yield from descente(i + 1, reste - v)
        mult[i] = 0

    yield from descente(0, reste)


def _compte_multisets_unite(m, taille, borne_mult):
    """Cardinal de `_multisets_unite`, par recurrence -- sans rien construire.

    Sert a decider AVANT enumeration si la strate tient sous le plafond.
    """
    dp = [1] + [0] * taille
    for _ in range(m):
        nouveau = [0] * (taille + 1)
        for s in range(taille + 1):
            if not dp[s]:
                continue
            for v in range(min(borne_mult, taille - s) + 1):
                nouveau[s + v] += dp[s]
        dp = nouveau
    return dp[taille]


def _b_depuis_mult(mult, m):
    """Multiplicites -> liste de vecteurs de charges."""
    b = []
    for a in range(m):
        for _ in range(mult[a]):
            q = [0] * m
            q[a] = 1
            b.append(q)
    return b


def familles_unite(m, r_B, max_charge, plafond=200_000, plafond_perturbe=20_000,
                   rng=None, n_echantillon=2000, stats=None, cle=None):
    """
    La famille « sommes de vecteurs unite, au plus un vecteur perturbe »,
    filtree par |c1(B)|_inf <= max_charge.

    Deux plafonds, parce que les deux strates n'ont ni la meme taille ni
    la meme importance :

    - `plafond` (strate pure, k = 0). Genereux : 200 000 suffit a rendre
      cette strate EXHAUSTIVE pour tout m <= 12 et tout r_B <= 7, donc
      pour la totalite des 194 CICYs a quotient libre (m <= 10). C'est la
      famille structurellement naturelle -- B = somme directe de O(e_a) --
      et celle d'ou sortent #6890 et #6947.
    - `plafond_perturbe` (strate k = 1). Plus serre : cette strate croit
      comme |Q| = 2m(m-1) fois la precedente et depasse le million des
      m = 8. 20 000 la garde exhaustive jusqu'a m = 6 (158 des 194 CICYs)
      et la rabat sur un echantillonnage DECLARE au-dela.

    Retourne la liste des b_charges. Renseigne `stats[cle]` avec, par
    strate : mode ('exhaustif' / 'echantillonne'), nombre produit, et
    nombre total quand il est connu.
    """
    infos = {}
    sortie = []

    # ---- strate k = 0 : que des vecteurs unite ------------------------
    # |c1(B)|_inf <= max_charge equivaut ici a « chaque multiplicite
    # <= max_charge », car c1(B)_k EST la multiplicite de e_k. Le filtre
    # est donc exact et applique dans la descente : rien n'est produit
    # pour etre ensuite rejete.
    n0 = _compte_multisets_unite(m, r_B, max_charge)
    if n0 <= plafond:
        for mult in _multisets_unite(m, r_B, max_charge):
            sortie.append(_b_depuis_mult(mult, m))
        infos['k0'] = {'mode': 'exhaustif', 'produit': n0, 'total': n0}
    else:
        # Jamais atteint pour les CICYs de la liste d'Oxford (m <= 19,
        # r_B <= 7 donne au plus quelques dizaines de milliers), mais on
        # ne laisse pas le cas silencieux.
        vus = set()
        r = rng if rng is not None else np.random.RandomState(0)
        for _ in range(n_echantillon):
            mult = [0] * m
            for _ in range(r_B):
                mult[int(r.randint(0, m))] += 1
            if max(mult) <= max_charge:
                vus.add(tuple(mult))
        for mult in sorted(vus):
            sortie.append(_b_depuis_mult(mult, m))
        infos['k0'] = {'mode': 'echantillonne', 'produit': len(vus),
                       'total': n0}

    # ---- strate k = 1 : un seul vecteur perturbe ----------------------
    Q = _vecteurs_perturbes(m)
    if m > 1 and r_B >= 1:
        # Un vecteur perturbe deplace chaque composante de +1 et de eps ;
        # la partie pure doit donc rester dans [-(max_charge+1), max_charge+1].
        n_pur = _compte_multisets_unite(m, r_B - 1, max_charge + 1)
        travail = n_pur * len(Q)
        if travail <= plafond_perturbe:
            for mult in _multisets_unite(m, r_B - 1, max_charge + 1):
                base = _b_depuis_mult(mult, m)
                for q in Q:
                    c1 = [mult[k] + q[k] for k in range(m)]
                    if all(abs(x) <= max_charge for x in c1):
                        sortie.append(base + [list(q)])
            infos['k1'] = {'mode': 'exhaustif', 'produit': travail,
                           'total': travail}
        else:
            vus = set()
            r = rng if rng is not None else np.random.RandomState(0)
            for _ in range(n_echantillon):
                mult = [0] * m
                for _ in range(r_B - 1):
                    mult[int(r.randint(0, m))] += 1
                q = Q[int(r.randint(0, len(Q)))]
                c1 = [mult[k] + q[k] for k in range(m)]
                if all(abs(x) <= max_charge for x in c1):
                    vus.add((tuple(mult), q))
            for mult, q in sorted(vus):
                sortie.append(_b_depuis_mult(list(mult), m) + [list(q)])
            infos['k1'] = {'mode': 'echantillonne', 'produit': len(vus),
                           'total': travail}

    # ---- strates k >= 2 : jamais enumerees ----------------------------
    infos['k2+'] = {'mode': 'non_couvert', 'produit': 0, 'total': None}

    if stats is not None:
        stats.setdefault('familles_unite', {})[cle or (m, r_B)] = infos
    return sortie


def _c_depuis_c1B(c1B, r_C, m):
    """Repartit c1(B) sur r_C fibres en droites, de sorte que c1(V) = 0."""
    if r_C == 1:
        return [list(c1B)]
    c1 = [c1B[k] // 2 for k in range(m)]
    c2 = [c1B[k] - c1[k] for k in range(m)]
    return [c1, c2]


def generate_monads(m, rank_V, max_charge=3, n_random=100, rng=None,
                    seed=42, plafond_unite=200_000, plafond_perturbe=20_000,
                    stats=None):
    """
    Genere des fibres monades candidats.

    Strategie :
    - r_C = 1 ou 2 (monades simples)
    - r_B = r_V + r_C
    - c1(V) = c1(B) - c1(C) = 0
    - Charges dans [-max_charge, max_charge]

    Deux proprietes acquises ici et absentes de la version precedente :

    - la famille des sommes de vecteurs unite est ENUMEREE (voir le bloc
      ci-dessus) au lieu d'etre tiree dix fois ;
    - le tirage residuel utilise un RNG DERIVE de (seed, m, max_charge,
      rank_V, r_C) et non plus le RNG partage avec le generateur positif.
      Une modification du generateur positif ne peut donc plus deplacer
      ce que celui-ci produit. C'est la lecon du §5.11, qui n'avait ete
      appliquee qu'a `generate_positive_monads`.

    `rng` est conserve dans la signature pour compatibilite avec les
    appelants existants, mais n'est plus utilise.
    """
    monads = []

    for r_C in [1, 2]:
        r_B = rank_V + r_C
        # RNG PROPRE a (seed, m, max_charge, rank_V, r_C).
        rng_local = np.random.RandomState(
            _mix(seed, (m, max_charge, rank_V, r_C), r_B))

        # Monades structurees
        # Type 1 : B = O(e_1) ⊕ ... ⊕ O(e_{rB}), C = O(sum_B)
        # Sous-cas des multiensembles de vecteurs unite ci-dessous ; garde
        # pour ne pas changer l'ordre de sortie des premieres monades.
        if r_C == 1:
            for shift in range(min(m, 3)):
                b_charges = []
                for r in range(r_B):
                    q = [0] * m
                    q[(r + shift) % m] = 1
                    b_charges.append(q)
                c1B = [sum(b[k] for b in b_charges) for k in range(m)]
                if all(abs(x) <= max_charge for x in c1B):
                    monads.append(MonadBundle(b_charges, [c1B]))

        # Type 2 : sommes de vecteurs unite -- ENUMEREES
        for b_charges in familles_unite(
                m, r_B, max_charge, plafond=plafond_unite,
                plafond_perturbe=plafond_perturbe, rng=rng_local,
                stats=stats, cle=(m, rank_V, r_C)):
            c1B = [sum(b[k] for b in b_charges) for k in range(m)]
            if all(abs(x) <= max_charge for x in c1B):
                monads.append(MonadBundle(
                    b_charges, _c_depuis_c1B(c1B, r_C, m)))

        # Monades aleatoires (charges quelconques dans [-max_charge, max_charge])
        for _ in range(n_random):
            b_charges = []
            for r in range(r_B):
                q = [int(rng_local.randint(-max_charge, max_charge + 1))
                     for _ in range(m)]
                b_charges.append(q)
            c1B = [sum(b[k] for b in b_charges) for k in range(m)]
            if any(abs(x) > max_charge * r_B for x in c1B):
                continue
            monads.append(MonadBundle(
                b_charges, _c_depuis_c1B(c1B, r_C, m)))

    return monads


# ======================================================================
# Version rigoureuse de la cohomologie d'une monade
# ======================================================================

def compute_monad_cohomology_ex(monad, ambient_dims, config_matrix,
                                stable_hypothesis=True):
    """
    H^i(V) pour 0 -> V -> B -> C -> 0, SANS hypothese de rang maximal.

    ----------------------------------------------------------------------
    Pourquoi cette fonction
    ----------------------------------------------------------------------
    `compute_monad_cohomology` ci-dessus deduit h^i(V) en supposant chaque
    application induite de rang maximal (rank_f0 = min(h0(B), h0(C)), etc.)
    puis rattrape le resultat avec chi. C'est exactement l'hypothese qui a
    ete retiree de `exact_cohomology` et de `monad_wedge`, et qui y produisait
    60 a 79 % de valeurs fausses. Elle etait restee ici.

    Consequence concrete : n_gen = |h1 - h2| est sur, puisqu'il vaut |chi(V)|
    des lors que h0 = h3 = 0. Mais la REPARTITION entre h1 et h2 -- donc
    n_anti = min(h1, h2), donc le tri des candidats « propres » -- venait de
    l'hypothese de rang maximal.

    ----------------------------------------------------------------------
    Ce que fait cette version
    ----------------------------------------------------------------------
    1. chi(V) = chi(B) - chi(C), EXACT (champ 'chi' de koszul_cohomology_ex).
    2. Bornes rigoureuses sur chaque h^i via la suite exacte longue :
           h^i(V) = coker(f_(i-1)) + ker(f_i)
           ker(f_i)   dans [max(0, h^i(B) - h^i(C)), h^i(B)]
           coker(f_j) dans [max(0, h^j(C) - h^j(B)), h^j(C)]
    3. Si `stable_hypothesis`, on impose h0(V) = h3(V) = 0 -- vrai pour tout
       fibre stable de pente nulle, et c'est le seul cas qui nous interesse.
       f_0 est alors injective (rang h0(B)) et coker(f_0) = h0(C) - h0(B)
       devient DETERMINE ; idem en degre 3. La contrainte chi relie ensuite
       h1 et h2, ce qui suffit souvent a fixer les deux.
    4. Certification : un degre n'est retenu que si son intervalle est un
       point ET que tous les fibres en droites intervenant sont certifies
       aux degres concernes.

    Retour : {0..3: h^i, 'chi': int, 'certified_by_degree': {i: bool},
              'bounds': {i: (lo, hi)}}
    Un degre non certifie porte la borne inferieure comme valeur ; se fier
    a `certified_by_degree` avant d'utiliser un h^i.
    """
    from cy_landscape.core.exact_cohomology import koszul_cohomology_ex

    def _somme(charges):
        h = {i: 0 for i in range(4)}
        cert = {i: True for i in range(4)}
        chi = 0
        for ch in charges:
            r = koszul_cohomology_ex(ambient_dims, config_matrix, list(ch))
            cd = r.get('certified_by_degree') or {}
            chi += int(r.get('chi', 0))
            for i in range(4):
                h[i] += int(r.get(i, 0))
                if not cd.get(i, False):
                    cert[i] = False
        return h, cert, chi

    try:
        hB, certB, chiB = _somme(monad.b_charges)
        hC, certC, chiC = _somme(monad.c_charges)
    except Exception:
        return None

    chiV = chiB - chiC

    lo = {}
    hi = {}
    for i in range(4):
        ker_lo, ker_hi = max(0, hB[i] - hC[i]), hB[i]
        if i == 0:
            cok_lo = cok_hi = 0
        else:
            cok_lo, cok_hi = max(0, hC[i - 1] - hB[i - 1]), hC[i - 1]
        lo[i], hi[i] = ker_lo + cok_lo, ker_hi + cok_hi

    if stable_hypothesis:
        # h0(V) = h3(V) = 0 : f_0 injective, f_3 injective.
        # Realisable seulement si h0(B) <= h0(C) et h3(B) <= h3(C).
        if hB[0] > hC[0] or hB[3] > hC[3]:
            return None                      # incompatible avec la stabilite
        lo[0] = hi[0] = 0
        lo[3] = hi[3] = 0
        # coker(f_0) = h0(C) - h0(B), determine
        cok0 = hC[0] - hB[0]
        lo[1] = cok0 + max(0, hB[1] - hC[1])
        hi[1] = cok0 + hB[1]

    # Contrainte chi : -h1 + h2 = chi(V) - h0 + h3 = chi(V) sous l'hypothese.
    # h2 = h1 + (h0 - h3 - chi) ; on croise les intervalles.
    d = (lo[0] - lo[3] - chiV) if stable_hypothesis else None
    if d is not None:
        n_lo = max(lo[1], lo[2] - d)
        n_hi = min(hi[1], hi[2] - d)
        if n_lo <= n_hi:
            lo[1], hi[1] = n_lo, n_hi
            lo[2], hi[2] = n_lo + d, n_hi + d

    cert = {}
    for i in range(4):
        ok = (lo[i] == hi[i])
        if ok and not (stable_hypothesis and i in (0, 3)):
            for src, idx in ((certB, i), (certC, i), (certB, i - 1), (certC, i - 1)):
                if idx >= 0 and not src.get(idx, False):
                    ok = False
                    break
        cert[i] = ok

    return {0: lo[0], 1: lo[1], 2: lo[2], 3: lo[3],
            'chi': chiV,
            'certified_by_degree': cert,
            'bounds': {i: (lo[i], hi[i]) for i in range(4)}}


def check_monad_nondegenerate(monad):
    """
    Rejette les monades structurellement degenerees, AVANT tout calcul.

    ----------------------------------------------------------------------
    Pourquoi
    ----------------------------------------------------------------------
    `check_map_exists` ne verifie qu'une chose : pour chaque c_j, il existe
    au moins un b_i avec H^0(O(c_j - b_i)) > 0. C'est trop faible. Trois
    situations passent ce filtre alors que V n'est pas un fibre stable :

    1. UN b_i EGAL A UN c_j
       f_i : O(b_i) -> O(c_j) est alors un isomorphisme, la monade se scinde
       et V devient une somme directe de fibres en droites. Une somme directe
       n'est jamais stable. Exemple reel rencontre sur la CICY 7669 :
       C = [-2, 0, -3] et b_1 = [-2, 0, -3], d'ou V = O(b_2) + O(b_3) + O(b_4).

    2. TROP PEU D'ENTREES NON NULLES DANS f
       f est une matrice rank_C x rank_B dont l'entree (j,i) vit dans
       H^0(O(c_j - b_i)). Pour que f soit surjective il faut au minimum
       rank_C entrees non nulles reparties sur des colonnes distinctes --
       sinon le rang de f est structurellement inferieur a rank_C. Sur le
       meme exemple, une seule entree sur quatre etait non nulle.

    3. UNE COLONNE ENTIEREMENT NULLE
       Si aucune entree de la colonne i n'est non nulle, O(b_i) est un
       facteur direct du noyau : V se scinde a nouveau.

    Ces trois tests ne coutent que des comparaisons de charges. Le critere
    de Hoppe ne les rattrape pas : pour des charges a signes melangees, tous
    les h^0(wedge^p V) peuvent s'annuler et le critere conclut a tort a la
    stabilite -- il est suffisant relativement a une classe de Kahler, et une
    somme directe est deja destabilisee par l'un de ses facteurs.

    Renvoie (ok, motif).
    """
    b = [list(x) for x in monad.b_charges]
    c = [list(x) for x in monad.c_charges]
    m = len(b[0])

    # 1. b_i == c_j
    for cj in c:
        for bi in b:
            if all(cj[k] == bi[k] for k in range(m)):
                return False, "b_i = c_j (monade scindee)"

    # matrice des entrees potentiellement non nulles : c_j - b_i >= 0
    nz = [[all(cj[k] - bi[k] >= 0 for k in range(m)) for bi in b] for cj in c]

    # 3. colonne nulle
    for i in range(len(b)):
        if not any(nz[j][i] for j in range(len(c))):
            return False, f"colonne {i} de f nulle (facteur direct)"

    # 2. rang structurel insuffisant : couplage maximal < rank_C
    if _couplage_max(nz, len(c), len(b)) < len(c):
        return False, "rang structurel de f < rank_C"

    return True, None


def _couplage_max(nz, nl, nc):
    """Couplage biparti maximal (Kuhn) sur le motif des entrees non nulles."""
    apparie = [-1] * nc

    def tente(j, vus):
        for i in range(nc):
            if nz[j][i] and not vus[i]:
                vus[i] = True
                if apparie[i] == -1 or tente(apparie[i], vus):
                    apparie[i] = j
                    return True
        return False

    total = 0
    for j in range(nl):
        if tente(j, [False] * nc):
            total += 1
    return total
