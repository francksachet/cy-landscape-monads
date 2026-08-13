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


def generate_monads(m, rank_V, max_charge=3, n_random=100, rng=None):
    """
    Genere des fibres monades candidats.

    Strategie :
    - r_C = 1 ou 2 (monades simples)
    - r_B = r_V + r_C
    - c1(V) = c1(B) - c1(C) = 0
    - Charges dans [-max_charge, max_charge]
    """
    if rng is None:
        rng = np.random.RandomState(42)

    monads = []

    for r_C in [1, 2]:
        r_B = rank_V + r_C

        # Monades structurees
        # Type 1 : B = O(e_1) ⊕ ... ⊕ O(e_{rB}), C = O(sum_B)
        if r_C == 1:
            # Quelques configurations simples
            for shift in range(min(m, 3)):
                b_charges = []
                for r in range(r_B):
                    q = [0] * m
                    q[(r + shift) % m] = 1
                    b_charges.append(q)
                c1B = [sum(b[k] for b in b_charges) for k in range(m)]
                c_charges = [c1B]  # c1(C) = c1(B) => c1(V) = 0
                if all(abs(x) <= max_charge for x in c1B):
                    monads.append(MonadBundle(b_charges, c_charges))

        # Type 2 : anti-symetriques
        for _ in range(min(m * 3, 10)):
            b_charges = []
            for r in range(r_B):
                q = [0] * m
                i1 = rng.randint(0, m)
                q[i1] = 1
                if m > 1:
                    i2 = (i1 + 1 + rng.randint(0, m - 1)) % m
                    q[i2] = rng.choice([-1, 0, 1])
                b_charges.append(q)
            c1B = [sum(b[k] for b in b_charges) for k in range(m)]
            if r_C == 1:
                c_charges = [c1B]
            else:
                # Repartir c1B sur 2 line bundles
                c1 = [c1B[k] // 2 for k in range(m)]
                c2 = [c1B[k] - c1[k] for k in range(m)]
                c_charges = [c1, c2]
            if all(abs(x) <= max_charge for x in c1B):
                monads.append(MonadBundle(b_charges, c_charges))

        # Monades aleatoires
        for _ in range(n_random):
            b_charges = []
            for r in range(r_B):
                q = [int(rng.randint(-max_charge, max_charge + 1)) for _ in range(m)]
                b_charges.append(q)
            c1B = [sum(b[k] for b in b_charges) for k in range(m)]
            if any(abs(x) > max_charge * r_B for x in c1B):
                continue
            if r_C == 1:
                c_charges = [c1B]
            else:
                c1 = [c1B[k] // 2 for k in range(m)]
                c2 = [c1B[k] - c1[k] for k in range(m)]
                c_charges = [c1, c2]
            monads.append(MonadBundle(b_charges, c_charges))

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
