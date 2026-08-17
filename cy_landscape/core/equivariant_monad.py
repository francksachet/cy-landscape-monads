"""
equivariant_monad.py -- Equivariance de f : B -> C. Le dernier verrou.

--------------------------------------------------------------------------
Ce qui est teste
--------------------------------------------------------------------------
Pour V = ker(f : B -> C), B = (+)_i O(b_i), C = (+)_j O(c_j), l'action de
Gamma sur X se releve a V si et seulement si f peut etre choisie
equivariante. La substitution S_g : x -> g.x envoie H^0(O(a)) sur
H^0(O(a o sigma^-1)), donc

    S_g( f_{j,i} )  =  lambda_g * f_{rho(j), pi(i)}                    (*)

ou pi et rho sont les permutations induites sur les facteurs de B et de C
(`permutation_charges`), et lambda_g un scalaire.

Les DEGRES se recollent : deg f_{j,i} = c_j - b_i, et
(c_j - b_i) o sigma^-1 = c_{rho(j)} - b_{pi(i)} = deg f_{rho(j), pi(i)}.
L'ancienne version de `gamma_action.espace_f_equivariant` ecrivait la
relation dans l'autre sens -- elle appliquait S_g a f_{rho(j),pi(i)} pour la
comparer a f_{j,i} -- ce qui donne un degre en sigma^-2 et non en sigma^-1.
Son propre garde-fou `if deg_img != degres[j][i]: return None` s'en
apercevait et faisait abandonner le calcul des que sigma n'etait pas une
involution. C'est pour cela qu'elle ne concluait jamais.

--------------------------------------------------------------------------
Le relevement projectif, traite et non contourne
--------------------------------------------------------------------------
lambda_g n'est PAS libre, et ce n'est pas non plus forcement une racine de
l'unite d'ordre |g| : le relevement de Gamma aux fibres en droites n'est
defini qu'a une phase pres, donc T_g^n peut valoir c.Id avec c != 1. On ne
suppose donc rien : on CALCULE le plus petit n tel que T_g^n soit scalaire,
on lit la constante c, et lambda_g parcourt les n racines n-iemes de c dans
GF(p). C'est exactement la bonne enumeration, et elle englobe le cas
Heisenberg de #7669 au lieu de le supposer absent.

Si aucune puissance de T_g jusqu'a `n_max` n'est scalaire, le module ne
conclut pas : il renvoie l'etat 'non scalaire', qui signale que l'espace
n'est pas isotypique et qu'une vraie decomposition spectrale serait
necessaire. Un cas non conclu n'est jamais compte comme un succes.

--------------------------------------------------------------------------
Prerequis : l'anneau doit etre covariant
--------------------------------------------------------------------------
S_g ne descend au quotient R_a = S_a / I_a que si S_g(I_a) = I_{a o sigma^-1}.
C'est faux pour les polynomes aleatoires de `sections.Ring`, vrai pour ceux
de `covariant_ring`. `verifier_descente` le controle explicitement plutot
que de le supposer, et le test refuse de tourner si le controle echoue.
"""
import numpy as np

from cy_landscape.core.sections import basis_multi, dim_multi
from cy_landscape.core.covariant_ring import (matrice_substitution, noyau_mod,
                                              rref_complet,
                                              permutation_facteurs_numerique)
from cy_landscape.core.gamma_action import permutation_charges


# ======================================================================
# Descente de S_g au quotient
# ======================================================================

def matrice_quotient(anneau, M, ambient, sigma, degre, p):
    """
    Matrice de S_g : R_degre -> R_{degre o sigma^-1}.

    Colonnes indexees par les monomes libres de S_degre, lignes par ceux de
    S_{deg_img} : on substitue puis on reduit modulo l'ideal d'arrivee.
    """
    A, deg_img = matrice_substitution(M, ambient, sigma, list(degre), p)
    S, idx, free, piv, Mred = anneau.quotient(list(degre))
    cols = []
    for k in free:
        v = np.zeros(A.shape[1], dtype=np.int64)
        v[k] = 1
        cols.append(anneau.reduce_vec(list(deg_img), (A @ v) % p))
    if not cols:
        return np.zeros((anneau.dimY(list(deg_img)), 0), dtype=np.int64), deg_img
    return (np.array(cols, dtype=np.int64).T % p), deg_img


def verifier_descente(anneau, M, ambient, sigma, degre, p):
    """
    S_g(I_degre) inclus dans I_{degre o sigma^-1} ?

    Controle direct : on prend une base de I_degre, on substitue, on reduit
    dans l'anneau d'arrivee, et le reste doit etre nul. C'est la condition
    qui rend `matrice_quotient` legitime ; sans polynomes covariants elle
    est fausse.
    """
    A, deg_img = matrice_substitution(M, ambient, sigma, list(degre), p)
    S, idx, free, piv, Mred = anneau.quotient(list(degre))
    if Mred.shape[0] == 0:
        return True, 0
    ecart = 0
    for r in range(Mred.shape[0]):
        img = (A @ (Mred[r] % p)) % p
        red = anneau.reduce_vec(list(deg_img), img)
        if red.size:
            ecart = max(ecart, int(red.max(initial=0)))
    return ecart == 0, ecart


# ======================================================================
# Ordre projectif d'un operateur
# ======================================================================

def _est_scalaire(T, p):
    """(True, c) si T = c.Id, sinon (False, None)."""
    n = T.shape[0]
    if T.shape[0] != T.shape[1] or n == 0:
        return False, None
    c = int(T[0, 0]) % p
    D = (T - c * np.eye(n, dtype=np.int64)) % p
    return (not D.any()), (c if not D.any() else None)


def ordre_projectif(T, p, n_max=64):
    """
    Plus petit n >= 1 tel que T^n = c.Id, et la constante c. (None, None) si
    aucun n <= n_max ne convient.
    """
    ok, c = _est_scalaire(T, p)
    if ok:
        return 1, c
    A = T % p
    for n in range(2, n_max + 1):
        A = (A @ T) % p
        ok, c = _est_scalaire(A, p)
        if ok:
            return n, c
    return None, None


def racines_niemes(c, n, p):
    """
    Toutes les x de GF(p) telles que x^n = c.

    On passe par le logarithme discret en base d'un generateur : avec
    c = g^t, l'equation devient n*k = t (mod p-1), soluble si et seulement si
    pgcd(n, p-1) divise t, et donnant alors exactement pgcd(n, p-1) racines.

    Une version anterieure fabriquait une « solution particuliere » par
    x0 = c^(inverse de n/d), ce qui n'est licite que si n est inversible
    modulo p-1 : elle rendait la liste vide sur des cas ou la racine existe
    (7^3 = 343 avec p = 30013). Les p en jeu valent quelques dizaines de
    milliers, donc la table des logarithmes est bon marche et exacte.
    """
    c = int(c) % p
    if c == 0:
        return [0]
    from math import gcd
    log = _table_log(p)
    t = log.get(c)
    if t is None:
        return []
    d = gcd(n, p - 1)
    if t % d:
        return []                       # pas de racine
    m = (p - 1) // d
    k0 = (t // d) * pow((n // d) % m, -1, m) % m
    g = _generateur(p)
    return sorted({pow(g, (k0 + j * m) % (p - 1), p) for j in range(d)})


_LOG = {}


def _table_log(p):
    """{g^k mod p : k}. Memoisee par p."""
    if p not in _LOG:
        g = _generateur(p)
        tab, x = {}, 1
        for k in range(p - 1):
            tab[x] = k
            x = x * g % p
        _LOG[p] = tab
    return _LOG[p]


_GEN = {}


def _generateur(p):
    if p in _GEN:
        return _GEN[p]
    from cy_landscape.core.gamma_action import _facteurs
    fac = _facteurs(p - 1)
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in fac):
            _GEN[p] = g
            return g
    raise RuntimeError("generateur introuvable")


# ======================================================================
# Espace des f equivariantes
# ======================================================================

def espace_f_equivariant(anneau, ambient, b_charges, c_charges,
                         coord_mats, p, maxdim=4000, n_max=64):
    """
    Espace des f = (f_{j,i}) verifiant (*) pour tous les generateurs.

    `anneau`  : CovariantRing sur le meme p (obligatoire -- voir en-tete).
    `coord_mats` : matrices numeriques des generateurs sur les coordonnees.

    Renvoie un dict :
      {'etat': 'ok' | motif d'echec,
       'dim_totale': N,
       'lambdas': [[valeurs possibles par generateur]],
       'solutions': [ {'lambda': (...), 'dim': d} ...],
       'dim_max': plus grande dimension trouvee,
       'existe': bool}
    """
    nB, nC = len(b_charges), len(c_charges)
    m = len(b_charges[0])

    degres = [[[c_charges[j][k] - b_charges[i][k] for k in range(m)]
               for i in range(nB)] for j in range(nC)]
    actif = [[all(x >= 0 for x in degres[j][i]) for i in range(nB)]
             for j in range(nC)]
    cases = [(j, i) for j in range(nC) for i in range(nB) if actif[j][i]]
    if not cases:
        return {'etat': 'aucune entree active', 'existe': False}

    dims, offs, acc = {}, {}, 0
    for (j, i) in cases:
        d = anneau.dimY(degres[j][i])
        dims[(j, i)] = d
        offs[(j, i)] = acc
        acc += d
    N = acc
    if N == 0:
        return {'etat': 'espace nul', 'existe': False}
    if N > maxdim:
        return {'etat': f'espace trop grand ({N})', 'existe': False}

    # sigma, pi, rho par generateur
    donnees = []
    for M in coord_mats:
        sigma = permutation_facteurs_numerique(M, ambient, p)
        if sigma is None:
            return {'etat': 'sigma non extractible', 'existe': False}
        pi = permutation_charges(b_charges, sigma)
        rho = permutation_charges(c_charges, sigma)
        if pi is None or rho is None:
            return {'etat': 'charges non permutees', 'existe': False}
        donnees.append((M, sigma, pi, rho))

    # controle : S_g descend bien au quotient
    for (M, sigma, pi, rho) in donnees:
        for (j, i) in cases:
            ok, ec = verifier_descente(anneau, M, ambient, sigma,
                                       degres[j][i], p)
            if not ok:
                return {'etat': f'S_g ne descend pas au quotient (ecart {ec})',
                        'existe': False}

    # operateur global T_g sur (+) R_{c_j - b_i}
    T = []
    for (M, sigma, pi, rho) in donnees:
        Tg = np.zeros((N, N), dtype=np.int64)
        possible = True
        for (j, i) in cases:
            jj, ii = rho[j], pi[i]
            A, deg_img = matrice_quotient(anneau, M, ambient, sigma,
                                          degres[j][i], p)
            if not actif[jj][ii]:
                if A.any():
                    possible = False       # image non nulle dans une case morte
                    break
                continue
            if list(deg_img) != list(degres[jj][ii]):
                return {'etat': 'degres incoherents', 'existe': False}
            Tg[offs[(jj, ii)]:offs[(jj, ii)] + dims[(jj, ii)],
               offs[(j, i)]:offs[(j, i)] + dims[(j, i)]] = A
        if not possible:
            return {'etat': 'case active envoyee sur une case nulle',
                    'existe': False}
        T.append(Tg % p)

    # valeurs propres admissibles, par ordre projectif
    lambdas = []
    for Tg in T:
        n, c = ordre_projectif(Tg, p, n_max=n_max)
        if n is None:
            return {'etat': 'operateur sans puissance scalaire', 'existe': False,
                    'dim_totale': N}
        rac = racines_niemes(c, n, p)
        if not rac:
            return {'etat': 'aucune racine dans GF(p)', 'existe': False,
                    'dim_totale': N}
        lambdas.append(rac)

    from itertools import product
    solutions = []
    for lam in product(*lambdas):
        blocs = [(T[g] - int(lam[g]) * np.eye(N, dtype=np.int64)) % p
                 for g in range(len(T))]
        base = noyau_mod(np.vstack(blocs), p)
        if base.shape[0]:
            solutions.append({'lambda': tuple(int(x) for x in lam),
                              'dim': int(base.shape[0]), 'base': base})
    return {
        'etat': 'ok',
        'dim_totale': N,
        'lambdas': lambdas,
        'solutions': solutions,
        'dim_max': max([s['dim'] for s in solutions], default=0),
        'existe': bool(solutions),
        'cases': cases, 'offsets': offs, 'dims': dims, 'degres': degres,
    }


# ======================================================================
# Stabilite RESTREINTE au sous-espace equivariant
# ======================================================================

def h0_V_sur_espace(anneau, ambient, b_charges, c_charges, base, cases,
                    offsets, dims, degres, p, rng, n_essais=5):
    """
    h^0(V) pour un f tire dans un sous-espace donne de l'espace des f.

    C'est le test qui MORD. L'existence d'un f equivariant est presque
    automatique des que Gamma est un Z2 agissant par phases : le sous-espace
    propre represente alors la moitie de l'espace, il n'est jamais vide, et
    en conclure « le fibre descend » serait une erreur. La question utile est
    de savoir si un f equivariant donne encore un fibre STABLE, c'est-a-dire
    si h^0(V) = 0 reste vrai quand f est contraint.

    On prend le MINIMUM sur plusieurs tirages : h^0 est semi-continu
    superieurement, donc le minimum observe est la meilleure approximation
    par le haut de la valeur generique du sous-espace.

    Renvoie (h0_min, dim_source).
    """
    from cy_landscape.core.sections import _mult_matrix, rref_mod
    dsrc = sum(anneau.dimY(list(b)) for b in b_charges)
    meilleur = None
    for _ in range(n_essais):
        v = (rng.randint(0, p, size=base.shape[0]) @ base) % p
        lignes = []
        for j, cj in enumerate(c_charges):
            ddst = anneau.dimY(list(cj))
            ligne = []
            for i, b in enumerate(b_charges):
                if (j, i) not in offsets:
                    ligne.append(np.zeros((ddst, anneau.dimY(list(b))),
                                          dtype=np.int64))
                    continue
                deg = degres[j][i]
                S, idx, free, piv, Mred = anneau.quotient(deg)
                coeffs = np.zeros(len(S), dtype=np.int64)
                for t, k in enumerate(free):
                    coeffs[k] = v[offsets[(j, i)] + t]
                ligne.append(_mult_matrix(anneau, list(b), deg, (S, coeffs),
                                          list(cj)))
            lignes.append(np.hstack(ligne))
        M = np.vstack(lignes)
        rang, _ = rref_mod(M.T.copy(), p)
        h0 = dsrc - rang
        meilleur = h0 if meilleur is None else min(meilleur, h0)
    return meilleur, dsrc


def _f_depuis_vecteur(anneau, v, offsets, dims, degres, j, i):
    """
    Coefficients de f_{j,i} sur la base MONOMIALE complete de S_{c_j - b_i}.

    Le vecteur `v` porte les coordonnees sur les monomes LIBRES du quotient ;
    les monomes pivots sont mis a zero, ce qui choisit un representant du
    meme element de R_a. Le resultat ne depend du representant que par un
    element de l'ideal, donc pas du tout apres multiplication et reduction.
    """
    deg = degres[j][i]
    S, idx, free, piv, Mred = anneau.quotient(deg)
    coeffs = np.zeros(len(S), dtype=np.int64)
    for t, k in enumerate(free):
        coeffs[k] = v[offsets[(j, i)] + t]
    return deg, S, coeffs


def h0_wedge2_V_sur_espace(anneau, b_charges, c_charges, base, offsets, dims,
                           degres, p, rng, n_essais=5, maxdim=6000):
    """
    h^0(wedge^2 V) pour un f tire dans un sous-espace donne de l'espace des f.

    Meme construction que `sections.h0_wedge2_V_explicit` -- wedge^2 V est le
    noyau de

        wedge^2 B --> B (x) C,   b1 ^ b2 |-> b1 (x) f(b2) - b2 (x) f(b1)

    -- mais f n'est plus tire au hasard dans tout H^0(O(c - b_i)) : il est
    tire dans le sous-espace equivariant. C'est la seconde moitie du critere
    de Hoppe au rang 4 (h^0(V) = h^0(w2V) = h^3(V) = 0), la premiere etant
    `h0_V_sur_espace`. Sans elle, un candidat de rang 4 ou 5 declare
    « survit » ne l'est que sur un tiers du critere.

    Restreint a rank_C = 1, comme tout le chemin wedge^2 du pipeline
    (`hoppe_fast` renvoie explicitement « non testable » au-dela).

    Renvoie (h0_min, dim_source) ou (None, dim_source) si la taille depasse
    `maxdim` -- auquel cas l'appelant ne doit PAS conclure.
    """
    from cy_landscape.core.sections import _mult_matrix, rref_mod
    if len(c_charges) != 1:
        return None, 0
    m = len(b_charges[0])
    c = list(c_charges[0])
    n = len(b_charges)
    paires = [(i, j) for i in range(n) for j in range(i + 1, n)]

    src = {(i, j): [b_charges[i][k] + b_charges[j][k] for k in range(m)]
           for (i, j) in paires}
    dims_src = {ij: anneau.dimY(src[ij]) for ij in paires}
    dsrc = sum(dims_src.values())
    dst = [[b_charges[k][t] + c[t] for t in range(m)] for k in range(n)]
    ddst = sum(anneau.dimY(d) for d in dst)
    if dsrc == 0:
        return 0, 0
    if max(dsrc, ddst) > maxdim:
        return None, dsrc

    offs_d, acc = {}, 0
    for k in range(n):
        offs_d[k] = acc
        acc += anneau.dimY(dst[k])

    meilleur = None
    for _ in range(n_essais):
        v = (rng.randint(0, p, size=base.shape[0]) @ base) % p
        M = np.zeros((ddst, dsrc), dtype=np.int64)
        oc = 0
        for (i, j) in paires:
            w = dims_src[(i, j)]
            for (k, l, sgn) in ((i, j, 1), (j, i, -1)):
                if (0, l) not in offsets:
                    continue                      # f_l nul (degre negatif)
                deg, S, coeffs = _f_depuis_vecteur(anneau, v, offsets, dims,
                                                   degres, 0, l)
                blk = _mult_matrix(anneau, src[(i, j)], deg, (S, coeffs),
                                   dst[k])
                if blk.size:
                    a, b = offs_d[k], offs_d[k] + blk.shape[0]
                    M[a:b, oc:oc + w] = (M[a:b, oc:oc + w] + sgn * blk) % p
            oc += w
        rang, _ = rref_mod(M.T.copy(), p)
        h0 = dsrc - rang
        meilleur = h0 if meilleur is None else min(meilleur, h0)
    return meilleur, dsrc


def h0_wedgep_V_sur_espace(anneau, b_charges, c_charges, p_ext, base, offsets,
                           dims, degres, p, rng, n_essais=5, maxdim=6000,
                           twist=None):
    """
    h^0(wedge^p V) pour un f tire dans un sous-espace donne. rank_C = 1.

    ----------------------------------------------------------------------
    Pourquoi une version generale en p
    ----------------------------------------------------------------------
    Le critere de Hoppe pour c1(V) = 0 s'enonce exactement

        V stable  <=>  h^0(wedge^p V) = 0  pour p = 1 .. rk-1

    et TOUTES ces quantites se calculent par le meme noyau. Pour C de rang 1,
    la resolution de wedge^p V donne

        0 -> wedge^p V -> wedge^p B -> wedge^{p-1} B (x) C -> ...

    d'ou, en sections globales,

        h^0(wedge^p V) = dim ker( H^0(wedge^p B) -> H^0(wedge^{p-1}B (x) C) )

    L'application est la contraction par f :

        e_{i_1} ^ ... ^ e_{i_p}  |->  sum_k (-1)^{k-1} f_{i_k} . e_{I \\ i_k}

    ----------------------------------------------------------------------
    Ce que cela debloque : h^3
    ----------------------------------------------------------------------
    det V = O donne wedge^{rk-1} V = V*, et la dualite de Serre sur un CY3
    donne h^0(V*) = h^3(V). Donc

        h^3(V)        = h^0(wedge^{rk-1} V)
        h^3(wedge^2 V) = h^0(wedge^{rk-2} V)      (rang 5)

    Le dernier tiers du critere devient donc calculable par ce chemin, y
    compris RESTREINT au sous-espace equivariant -- ce que `hoppe_fast` ne
    pouvait pas faire, traitant h^3 par les bornes de
    `compute_monad_cohomology_ex`, sans f explicite.

    ATTENTION a l'identite utilisee : wedge^{rk-1}V = V* vaut a TOUT rang,
    mais wedge^3 V = V* seulement au rang 4. C'est la confusion de
    `_wedge3_h0_twisted` signalee au §4.5. Ici on indexe par p_ext et on ne
    suppose rien : l'appelant choisit p_ext = rk-1 pour h^3(V).

    ----------------------------------------------------------------------
    Controles disponibles
    ----------------------------------------------------------------------
      p_ext = 1  -> doit redonner h0_V_sur_espace
      p_ext = 2  -> doit redonner h0_wedge2_V_sur_espace
      p_ext = rk -> wedge^{rk} V = det V = O, donc h^0 = h^0(O_Y) = 1
    Le dernier teste la construction entiere -- sources, cibles et SIGNES --
    contre une valeur connue d'avance.

    ----------------------------------------------------------------------
    Le parametre `twist`
    ----------------------------------------------------------------------
    `twist = H` calcule h^0(wedge^p V (-H)) au lieu de h^0(wedge^p V). La
    resolution est la meme tordue par O(-H) :

        0 -> wedge^p V(-H) -> wedge^p B(-H) -> wedge^{p-1}B (x) C(-H) -> ...

    et la contraction par f est INCHANGEE -- f_{i} est de degre c - b_i, que
    le decalage de la source et de la cible par le meme -H ne modifie pas.
    Il suffit donc de retrancher H aux deux degres.

    A quoi cela sert : le critere de Hoppe sous la forme h^0(wedge^p V) = 0
    (twist nul) est une EQUIVALENCE avec la stabilite seulement si Pic(X) est
    de rang 1. En rang superieur la forme suffisante demande
    h^0(wedge^p V(-H)) = 0 pour tout H de degre positif ou nul (§5.14).

    Renvoie (h0_min, dim_source), ou (None, ...) si hors taille : l'appelant
    ne doit alors PAS conclure.
    """
    from itertools import combinations
    from cy_landscape.core.sections import _mult_matrix, rref_mod
    if len(c_charges) != 1:
        return None, 0
    n = len(b_charges)
    m = len(b_charges[0])
    c = list(c_charges[0])
    if not (1 <= p_ext <= n):
        return None, 0

    src_idx = list(combinations(range(n), p_ext))
    dst_idx = list(combinations(range(n), p_ext - 1))
    tw = [0] * m if twist is None else [int(x) for x in twist]
    src_deg = {I: [sum(b_charges[i][k] for i in I) - tw[k] for k in range(m)]
               for I in src_idx}
    dst_deg = {J: [sum(b_charges[j][k] for j in J) + c[k] - tw[k]
                   for k in range(m)]
               for J in dst_idx}

    dsrc_par = {I: anneau.dimY(src_deg[I]) for I in src_idx}
    dsrc = sum(dsrc_par.values())
    ddst_par = {J: anneau.dimY(dst_deg[J]) for J in dst_idx}
    ddst = sum(ddst_par.values())
    if dsrc == 0:
        return 0, 0
    if max(dsrc, ddst) > maxdim:
        return None, dsrc

    off_s, a = {}, 0
    for I in src_idx:
        off_s[I] = a
        a += dsrc_par[I]
    off_d, a = {}, 0
    for J in dst_idx:
        off_d[J] = a
        a += ddst_par[J]

    meilleur = None
    for _ in range(n_essais):
        v = (rng.randint(0, p, size=base.shape[0]) @ base) % p
        fpol = {}
        for i in range(n):
            if (0, i) not in offsets:
                continue
            deg = degres[0][i]
            S, idx, free, piv, Mred = anneau.quotient(deg)
            coeffs = np.zeros(len(S), dtype=np.int64)
            for t, k in enumerate(free):
                coeffs[k] = v[offsets[(0, i)] + t]
            fpol[i] = (deg, S, coeffs)

        M = np.zeros((ddst, dsrc), dtype=np.int64)
        for I in src_idx:
            for k, ik in enumerate(I):
                if ik not in fpol:
                    continue                      # f_{ik} nul (degre negatif)
                J = tuple(x for x in I if x != ik)
                deg, S, coeffs = fpol[ik]
                blk = _mult_matrix(anneau, src_deg[I], deg, (S, coeffs),
                                   dst_deg[J])
                if not blk.size:
                    continue
                signe = 1 if k % 2 == 0 else -1
                r0, c0 = off_d[J], off_s[I]
                M[r0:r0 + blk.shape[0], c0:c0 + blk.shape[1]] = (
                    M[r0:r0 + blk.shape[0], c0:c0 + blk.shape[1]]
                    + signe * blk) % p
        rang, _ = rref_mod(M.T.copy(), p)
        h0 = dsrc - rang
        meilleur = h0 if meilleur is None else min(meilleur, h0)
    return meilleur, dsrc


def hoppe_suffisant_generique(ambient, config, b_charges, c_charges, D,
                              seed=2, rng_seed=5, maxdim=6000,
                              plafond_twists=200000):
    """
    Critere de Hoppe suffisant pour un f GENERIQUE de l'anneau, sans
    contrainte d'equivariance. Monte l'anneau et la base pleine, puis
    delegue a `hoppe_suffisant_sur_espace`.

    Portee de l'enonce : `stable = True` signifie que le membre GENERIQUE
    de la famille de monades est un fibre mu_J-stable. Il ne dit rien d'un
    f particulier -- notamment pas d'un f equivariant, h^0 etant
    semi-continu superieurement : la valeur ne peut que MONTER en un point
    special. Pour un candidat porteur d'une structure equivariante, c'est
    `hoppe_suffisant_sur_espace` avec la base contrainte qui fait foi
    (§5.14).

    Renvoie None si la monade est hors du domaine du modele S/I, ou si
    rank_C != 1 -- dans les deux cas on ne conclut pas.
    """
    from cy_landscape.core.sections import Ring, P, domaine_valide
    if len(c_charges) != 1:
        return None
    if not domaine_valide(ambient, config, b_charges, c_charges):
        return None
    m = len(ambient)
    degres = [[[c_charges[0][k] - b_charges[i][k] for k in range(m)]
               for i in range(len(b_charges))]]
    R = Ring(ambient, config, seed=seed)
    cases = [(0, i) for i in range(len(b_charges))
             if all(x >= 0 for x in degres[0][i])]
    dims, offsets, a = {}, {}, 0
    for k in cases:
        dims[k] = R.dimY(degres[0][k[1]])
        offsets[k] = a
        a += dims[k]
    if a == 0:
        return None
    base = np.eye(a, dtype=np.int64)
    return hoppe_suffisant_sur_espace(R, b_charges, c_charges, base, offsets,
                                      dims, degres, P,
                                      np.random.RandomState(rng_seed), D,
                                      maxdim=maxdim,
                                      plafond_twists=plafond_twists)


def hoppe_suffisant_sur_espace(anneau, b_charges, c_charges, base, offsets,
                               dims, degres, p, rng, D, maxdim=6000,
                               plafond_twists=200000):
    """
    Critere de Hoppe sous sa forme SUFFISANTE, restreint a un sous-espace.

    ----------------------------------------------------------------------
    Ce que `hoppe_sur_espace` ne dit pas
    ----------------------------------------------------------------------
    « c1(V) = 0 => V stable <=> h0(wedge^p V) = 0 » est une EQUIVALENCE
    seulement si Pic(X) est de rang 1. Sur une CICY a m > 1 -- et #6890
    comme #6947 ont h11 = 5 -- elle reste NECESSAIRE sans etre suffisante :
    elle ne voit pas la classe de Kahler. Un verdict `stable: True` de
    `hoppe_sur_espace` se lit donc « non elimine ».

    La forme suffisante ajoute les torsions : V est mu_J-stable des que

        h0(wedge^p V(-H)) = 0   pour p = 1..rk-1 et tout H avec
                                deg_J(H) >= 0

    car un sous-faisceau destabilisant de wedge^p V fournirait une section
    non nulle de wedge^p V(-H) pour un tel H. L'ensemble des H est fini et
    petit (`hoppe_fast.polytope_twists`).

    Le twist H = 0 est inclus dans le polytope : ce critere CONTIENT donc
    `hoppe_sur_espace`, il ne s'y substitue pas.

    ----------------------------------------------------------------------
    Trois issues, jamais confondues
    ----------------------------------------------------------------------
      True  -> stable pour CETTE classe de Kahler J. DEMONTRE.
      False -> un H de degre >= 0 donne une section : V n'est pas stable
               pour J, et le H en question est le temoin.
      None  -> au moins un h0 hors de portee, ou polytope non borne /
               au-dela du plafond. On ne conclut pas -- et surtout on ne
               compte pas cela comme un succes : un critere suffisant n'a
               de valeur que s'il est verifie EN ENTIER.
    """
    from cy_landscape.core.hoppe_fast import polytope_twists
    rk = len(b_charges) - len(c_charges)
    detail, incomplet = {}, []
    for q in range(1, rk):
        Hs = polytope_twists(b_charges, q, D, plafond=plafond_twists)
        if Hs is None:
            return {'stable': None, 'motif': f'polytope des twists non borne '
                                             f'ou hors plafond a p = {q}',
                    'detail': detail, 'n_twists': None}
        # `sources_non_vides` mesure si le verdict a du CONTENU. Un twist
        # dont la source est vide donne h0 = 0 sans qu'aucun rang ne soit
        # calcule : un critere dont tous les twists seraient vides serait
        # vrai sans rien prouver. L'appelant doit pouvoir le voir.
        # `source_H0` sert de temoin que le twist AGIT : si le parametre
        # etait ignore, toutes les sources vaudraient celle de H = 0 et
        # `source_max` retomberait dessus. Le test de regression l'exige
        # strictement superieur -- sans quoi un twist neutralise passerait
        # inapercu, tous les h0 valant 0 de toute facon sur les candidats.
        detail[q] = {'n_twists': len(Hs), 'non_nuls': [], 'indetermines': 0,
                     'sources_non_vides': 0, 'source_max': 0,
                     'source_H0': None}
        for H in Hs:
            h, dsrc = h0_wedgep_V_sur_espace(anneau, b_charges, c_charges, q,
                                             base, offsets, dims, degres, p,
                                             rng, maxdim=maxdim, twist=H)
            if dsrc > 0:
                detail[q]['sources_non_vides'] += 1
                detail[q]['source_max'] = max(detail[q]['source_max'], dsrc)
            if not any(H):
                detail[q]['source_H0'] = dsrc
            if h is None:
                detail[q]['indetermines'] += 1
                incomplet.append((q, tuple(H)))
            elif h > 0:
                detail[q]['non_nuls'].append((list(H), int(h)))
                return {'stable': False,
                        'motif': f'h0(wedge^{q} V(-H)) = {h} avec H = {H}, '
                                 f'deg_J(H) >= 0',
                        'detail': detail, 'temoin': list(H)}
    if incomplet:
        return {'stable': None,
                'motif': f'{len(incomplet)} twists hors de portee '
                         f'(maxdim) -- critere incomplet',
                'detail': detail, 'incomplet': incomplet[:5]}
    n = sum(v['n_twists'] for v in detail.values())
    nv = sum(v['sources_non_vides'] for v in detail.values())
    return {'stable': True,
            'motif': f'stable pour cette classe de Kahler : {n} twists de '
                     f'degre >= 0 testes, tous a h0 nul, dont {nv} a source '
                     f'non vide',
            'detail': detail, 'n_twists': n, 'sources_non_vides': nv}


def hoppe_sur_espace(anneau, b_charges, c_charges, base, offsets, dims,
                     degres, p, rng, maxdim=6000):
    """
    Critere de Hoppe COMPLET, restreint a un sous-espace de f. rank_C = 1.

    Teste h^0(wedge^p V) = 0 pour p = 1 .. rk-1, ce qui est l'enonce meme du
    critere pour c1(V) = 0 -- et non un sous-ensemble comme les phases 1 et 2
    de `hoppe_fast` (§4.5). Le tableau des isomorphismes du §4.5 n'est plus
    necessaire : au lieu de traduire h^3(V) en h^0(wedge^{rk-1}V) puis de
    chercher un chemin de calcul, on calcule directement tous les
    h^0(wedge^p V).

    Renvoie {'stable': bool | None, 'valeurs': {p: h0}, 'motif': str}.
    `stable = None` signifie qu'au moins un p n'a pas pu etre calcule : on ne
    conclut pas, et ce n'est jamais compte comme un succes.
    """
    rk = len(b_charges) - len(c_charges)
    valeurs, incomplet = {}, False
    for q in range(1, rk):
        h, _ = h0_wedgep_V_sur_espace(anneau, b_charges, c_charges, q, base,
                                      offsets, dims, degres, p, rng,
                                      maxdim=maxdim)
        valeurs[q] = h
        if h is None:
            incomplet = True
        elif h > 0:
            return {'stable': False, 'valeurs': valeurs,
                    'motif': f'h0(wedge^{q} V) = {h}'}
    if incomplet:
        return {'stable': None, 'valeurs': valeurs,
                'motif': 'au moins un degre non calculable'}
    return {'stable': True, 'valeurs': valeurs, 'motif': ''}


def multiplicites_propres(T, p, valeurs, base=None):
    """
    Multiplicites des valeurs propres de T, eventuellement restreint au
    sous-espace engendre par les lignes de `base`.

    p etant impair et l'ordre du groupe fini, la representation est
    semi-simple : les multiplicites sont donc ADDITIVES le long de
    0 -> W -> R -> R/W -> 0, ce qui evite d'avoir a choisir un supplementaire
    pour decrire le quotient. Le controle que la somme des multiplicites vaut
    la dimension totale valide cette semi-simplicite au lieu de la supposer.
    """
    if base is None:
        A = np.asarray(T, dtype=np.int64) % p
        n = A.shape[0]
    else:
        R, piv = rref_complet(np.asarray(base, dtype=np.int64), p)
        n = R.shape[0]
        if n == 0:
            return {v: 0 for v in valeurs}
        img = (R @ np.asarray(T, dtype=np.int64).T) % p
        # R etant en forme echelonnee REDUITE, la coordonnee d'un vecteur sur
        # la ligne r est simplement sa valeur a la colonne pivot piv[r].
        A = np.array([[int(img[r, cpiv]) % p for cpiv in piv] for r in range(n)],
                     dtype=np.int64)
    out = {}
    for lam in valeurs:
        M = (A - (int(lam) % p) * np.eye(n, dtype=np.int64)) % p
        R2, _ = rref_complet(M, p)
        out[lam] = n - R2.shape[0]
    return out


def decomposition_h1_V_abelien(anneau, ambient, b_charges, c_charges,
                               coord_mats, base, offsets, dims, degres, p, rng):
    """
    Decomposition de H^1(V) sous un Gamma ABELIEN a plusieurs generateurs.

    ----------------------------------------------------------------------
    Pourquoi cette fonction existe
    ----------------------------------------------------------------------
    `decomposition_h1_V` commence par `if len(coord_mats) != 1: return None`
    -- elle ne traite que Gamma cyclique. Or les seuls candidats du scan
    qui echappent a l'argument de rang du §5.8 sont ceux a Gamma = Z2 x Z2 :
    |Gamma| = 4 donne DEUX lignes de Wilson, la ou une seule plafonne a
    Pati-Salam. Sans cette extension, ce sont precisement les candidats
    interessants dont le spectre reste inconnu.

    ----------------------------------------------------------------------
    Ce qui est calcule
    ----------------------------------------------------------------------
    Les generateurs agissent sur R_c par des operateurs T_g qui COMMUTENT
    (a un scalaire pres) et sont d'ordre projectif fini. On decompose donc
    R_c en sous-espaces propres SIMULTANES, un par caractere de Gamma, au
    moyen des projecteurs

        P_chi = prod_g ( sum_{k} chi(g)^{-k} T_g^k / n_g )

    puis on soustrait les multiplicites de l'image de f. La difference est
    la decomposition de H^1(V) = coker(H^0(B) -> H^0(C)).

    ----------------------------------------------------------------------
    Ce qui est VERIFIE plutot que suppose
    ----------------------------------------------------------------------
    - que les T_g commutent EXACTEMENT (et non a un scalaire pres) : sinon
      la decomposition simultanee n'existe pas, et on rend `None` plutot
      qu'un tableau de nombres sans signification ;
    - que chaque T_g est d'ordre projectif fini, avec sa constante ;
    - que la somme des multiplicites vaut la dimension -- sur R_c ET sur
      l'image. C'est le controle de semi-simplicite, et c'est lui qui a
      detecte le cocycle du relevement projectif sur (1,1,1).

    Reference independante : Gamma agissant librement,
    n_gen(X/Gamma) = |chi(V)| / |Gamma|. Avec h^1(V) = 12 et |Gamma| = 4,
    la partie INVARIANTE doit valoir exactement 3. Une decomposition en
    4+4+2+2 signalerait une erreur, pas une decouverte.
    """
    from cy_landscape.core.sections import _mult_matrix
    if len(c_charges) != 1 or not coord_mats:
        return None
    c = list(c_charges[0])

    # --- image de f dans R_c ------------------------------------------
    v = (rng.randint(0, p, size=base.shape[0]) @ base) % p
    blocs = []
    for i, bb in enumerate(b_charges):
        if (0, i) not in offsets:
            blocs.append(np.zeros((anneau.dimY(c), anneau.dimY(list(bb))),
                                  dtype=np.int64))
            continue
        deg = degres[0][i]
        S, idx, free, piv, Mred = anneau.quotient(deg)
        coeffs = np.zeros(len(S), dtype=np.int64)
        for t, k in enumerate(free):
            coeffs[k] = v[offsets[(0, i)] + t]
        blocs.append(_mult_matrix(anneau, list(bb), deg, (S, coeffs), c))
    Mf = np.hstack(blocs) % p
    W, pivW = rref_complet(Mf.T.copy(), p)
    dimC = anneau.dimY(c)
    h1 = dimC - W.shape[0]

    # --- operateurs des generateurs sur R_c ---------------------------
    Ts, ordres, consts = [], [], []
    for M in coord_mats:
        sigma = permutation_facteurs_numerique(M, ambient, p)
        T, deg_img = matrice_quotient(anneau, M, ambient, sigma, c, p)
        if list(deg_img) != c:
            return None
        n, cst = ordre_projectif(T, p)
        if n is None:
            return None
        Ts.append(np.asarray(T, dtype=np.int64) % p)
        ordres.append(n)
        consts.append(cst)

    # Commutation EXACTE. A un scalaire pres ne suffit pas : un cocycle non
    # trivial interdit la diagonalisation simultanee, et rendre malgre tout
    # des multiplicites donnerait des nombres qui ne decrivent rien.
    for a in range(len(Ts)):
        for b_ in range(a + 1, len(Ts)):
            if not np.array_equal((Ts[a] @ Ts[b_]) % p, (Ts[b_] @ Ts[a]) % p):
                return {'h1': h1, 'coherent': False,
                        'motif': 'les generateurs ne commutent pas sur R_c '
                                 '(cocycle du relevement projectif)'}

    # --- projecteurs simultanes ---------------------------------------
    valeurs = [racines_niemes(cst, n, p) for cst, n in zip(consts, ordres)]
    if any(not v_ for v_ in valeurs):
        return None
    dim = dimC
    I = np.eye(dim, dtype=np.int64)

    def projecteur(T, n, lam):
        """(1/n) sum_k lam^-k T^k -- projecteur sur le sous-espace propre."""
        inv_lam = pow(int(lam), -1, p)
        P = np.zeros((dim, dim), dtype=np.int64)
        A = I.copy()
        for k in range(n):
            P = (P + pow(inv_lam, k, p) * A) % p
            A = (A @ T) % p
        return (P * pow(n, -1, p)) % p

    import itertools as _it
    sur_Rc, sur_W = {}, {}
    for combo in _it.product(*valeurs):
        P = I.copy()
        for T, n, lam in zip(Ts, ordres, combo):
            P = (P @ projecteur(T, n, lam)) % p
        R, _ = rref_complet(P.T.copy(), p)
        sur_Rc[combo] = R.shape[0]
        if W.shape[0]:
            RW, _ = rref_complet((W @ P.T) % p, p)
            sur_W[combo] = RW.shape[0]
        else:
            sur_W[combo] = 0

    coherent = (sum(sur_Rc.values()) == dimC
                and sum(sur_W.values()) == W.shape[0])
    return {
        'h1': h1,
        'ordres': ordres, 'constantes': consts,
        'sur_Rc': sur_Rc, 'sur_image': sur_W,
        'sur_H1': {k: sur_Rc[k] - sur_W[k] for k in sur_Rc},
        'coherent': coherent,
    }


def decomposition_h1_V(anneau, ambient, b_charges, c_charges, coord_mats,
                       base, offsets, dims, degres, p, rng, valeurs=None):
    """
    Decomposition de H^1(V) en representations de Gamma.

    ----------------------------------------------------------------------
    Pourquoi H^1(V) est ici un conoyau explicite
    ----------------------------------------------------------------------
    La suite longue de 0 -> V -> B -> C -> 0 donne

        H^0(B) -> H^0(C) -> H^1(V) -> H^1(B)

    Si H^1(B) = 0 -- condition a VERIFIER par `koszul_cohomology_ex`, pas a
    supposer parce que les charges sont positives -- alors

        H^1(V) = coker( H^0(B) -> H^0(C) )

    c'est-a-dire R_c modulo l'image de f. Gamma agit sur R_c par S_g, et
    l'image est stable puisque f est equivariante : l'action descend donc au
    conoyau.

    ----------------------------------------------------------------------
    Reference independante
    ----------------------------------------------------------------------
    Pour Gamma agissant librement, n_gen(X/Gamma) = n_gen(X)/|Gamma|. Avec
    h^1(V) = 6 et |Gamma| = 2, la partie invariante doit valoir exactement 3.
    C'est une valeur connue D'AVANCE : une decomposition en 4+2 ou 5+1
    signalerait une erreur dans l'action, dans le conoyau, ou dans le fibre.

    Deux controles internes accompagnent le calcul : la somme des
    multiplicites doit valoir dim R_c d'une part, dim de l'image d'autre part
    (semi-simplicite effective).

    Renvoie {'h1': , 'sur_Rc': , 'sur_image': , 'sur_H1': , 'coherent': }.
    """
    from cy_landscape.core.sections import _mult_matrix
    if len(coord_mats) != 1:
        return None                    # ce chemin ne traite que Gamma cyclique
    if valeurs is None:
        valeurs = [1, p - 1]           # Z2

    # RANK_C QUELCONQUE. La cible n'est plus R_c mais la SOMME DIRECTE
    # (+)_j R_{c_j}, et Gamma peut PERMUTER les facteurs : l'operateur n'est
    # alors pas diagonal par blocs. Cette generalisation n'est pas cosmetique
    # -- les quatre candidats Z4 du §5.28, seule route vers le rang 4 exempte
    # du cocycle du §5.27, ont tous rank_C = 2, et c'est le `!= 1` d'origine
    # qui les laissait sans nombre de generations.
    cs = [list(x) for x in c_charges]
    nC = len(cs)
    dimc = [anneau.dimY(cj) for cj in cs]
    depart = [sum(dimc[:j]) for j in range(nC)]
    dimC = sum(dimc)

    v = (rng.randint(0, p, size=base.shape[0]) @ base) % p
    rangees = []
    for j, cj in enumerate(cs):
        blocs = []
        for i, bb in enumerate(b_charges):
            if (j, i) not in offsets:
                # case de degre negatif : f_{j,i} = 0 (§5.28)
                blocs.append(np.zeros((dimc[j], anneau.dimY(list(bb))),
                                      dtype=np.int64))
                continue
            deg = degres[j][i]
            S, idx, free, piv, Mred = anneau.quotient(deg)
            coeffs = np.zeros(len(S), dtype=np.int64)
            for t, k in enumerate(free):
                coeffs[k] = v[offsets[(j, i)] + t]
            blocs.append(_mult_matrix(anneau, list(bb), deg, (S, coeffs), cj))
        rangees.append(np.hstack(blocs) % p)
    Mf = np.vstack(rangees) % p

    W, pivW = rref_complet(Mf.T.copy(), p)          # lignes = image dans (+)R_c
    h1 = dimC - W.shape[0]

    sigma = permutation_facteurs_numerique(coord_mats[0], ambient, p)
    rho = permutation_charges(cs, sigma)
    if rho is None:
        return None                    # Gamma ne permute pas les c_j
    T = np.zeros((dimC, dimC), dtype=np.int64)
    for j, cj in enumerate(cs):
        A, deg_img = matrice_quotient(anneau, coord_mats[0], ambient, sigma,
                                      cj, p)
        jj = rho[j]
        if list(deg_img) != cs[jj]:
            return None
        T[depart[jj]:depart[jj] + dimc[jj],
          depart[j]:depart[j] + dimc[j]] = A
    T %= p

    tot = multiplicites_propres(T, p, valeurs)
    sur_W = multiplicites_propres(T, p, valeurs, base=W)
    coherent = (sum(tot.values()) == dimC
                and sum(sur_W.values()) == W.shape[0])
    return {
        'h1': h1,
        'sur_Rc': tot,
        'sur_image': sur_W,
        'sur_H1': {k: tot[k] - sur_W[k] for k in valeurs},
        'coherent': coherent,
    }


def decomposition_h1_wedge2_V(anneau, ambient, b_charges, c_charges,
                              coord_mats, base, offsets, dims, degres, p, rng,
                              valeurs=None):
    """
    H^1(wedge^2 V) et sa decomposition sous Gamma. rank_C = 1, Gamma cyclique.

    ----------------------------------------------------------------------
    Deux suites courtes, pas une resolution longue
    ----------------------------------------------------------------------
    La filtration de wedge^2 B et la suite de la monade tensorisee par C
    donnent

        0 -> wedge^2 V -> wedge^2 B -> V (x) C -> 0
        0 -> V (x) C   -> B (x) C   -> C^2     -> 0

    -- la seconde identifiant V (x) C au noyau de B (x) C -> C^2. D'ou

        H^0(V (x) C)   = ker( H^0(B(x)C) -> H^0(C^2) )
        H^1(wedge^2 V) = coker( H^0(wedge^2 B) -> H^0(V (x) C) )

    la seconde egalite exigeant H^1(wedge^2 B) = 0, A CERTIFIER par
    `koszul_cohomology_ex` sur chaque b_i + b_j -- pas a supposer.

    Physiquement, pour V de rang 4 (groupe de structure SU(4), commutant
    SO(10)), h^1(wedge^2 V) compte les **10** de SO(10), d'ou viennent les
    doublets de Higgs apres brisure.

    ----------------------------------------------------------------------
    Controles inclus
    ----------------------------------------------------------------------
      - beta o alpha = 0 : c'est un complexe, sinon la construction est fausse ;
      - h^0(wedge^2 V) = dim ker(alpha) doit retomber sur la valeur donnee par
        `h0_wedge2_V_sur_espace`, calculee par un tout autre chemin ;
      - additivite des multiplicites (semi-simplicite).

    Renvoie {'h1', 'h0_controle', 'complexe_ok', 'sur_H1', 'coherent', ...}.
    """
    from itertools import combinations
    from cy_landscape.core.sections import _mult_matrix
    if len(c_charges) != 1 or len(coord_mats) != 1:
        return None
    c = list(c_charges[0])
    m = len(b_charges[0])
    n = len(b_charges)
    if valeurs is None:
        valeurs = [1, p - 1]

    v = (rng.randint(0, p, size=base.shape[0]) @ base) % p
    fpol = {}
    for i in range(n):
        if (0, i) not in offsets:
            continue
        deg = degres[0][i]
        S, idx, free, piv, Mred = anneau.quotient(deg)
        cf = np.zeros(len(S), dtype=np.int64)
        for t, k in enumerate(free):
            cf[k] = v[offsets[(0, i)] + t]
        fpol[i] = (deg, S, cf)

    paires = list(combinations(range(n), 2))
    dA = {I: [b_charges[I[0]][k] + b_charges[I[1]][k] for k in range(m)]
          for I in paires}
    dBC = [[b_charges[i][k] + c[k] for k in range(m)] for i in range(n)]
    dC2 = [2 * c[k] for k in range(m)]

    oA, a = {}, 0
    for I in paires:
        oA[I] = a
        a += anneau.dimY(dA[I])
    nA = a
    oB, a = {}, 0
    for i in range(n):
        oB[i] = a
        a += anneau.dimY(dBC[i])
    nBC = a
    nC2 = anneau.dimY(dC2)

    alpha = np.zeros((nBC, nA), dtype=np.int64)
    for I in paires:
        for k, ik in enumerate(I):
            if ik not in fpol:
                continue
            J = [x for x in I if x != ik][0]
            deg, S, cf = fpol[ik]
            blk = _mult_matrix(anneau, dA[I], deg, (S, cf), dBC[J])
            if blk.size:
                sg = 1 if k % 2 == 0 else -1
                r0, c0 = oB[J], oA[I]
                alpha[r0:r0 + blk.shape[0], c0:c0 + blk.shape[1]] = (
                    alpha[r0:r0 + blk.shape[0], c0:c0 + blk.shape[1]]
                    + sg * blk) % p

    beta = np.zeros((nC2, nBC), dtype=np.int64)
    for i in range(n):
        if i not in fpol:
            continue
        deg, S, cf = fpol[i]
        blk = _mult_matrix(anneau, dBC[i], deg, (S, cf), dC2)
        if blk.size:
            beta[:, oB[i]:oB[i] + blk.shape[1]] = blk % p

    complexe_ok = not ((beta @ alpha) % p).any()
    K = noyau_mod(beta, p)
    imA, _ = rref_complet(alpha.T.copy(), p)
    rg = imA.shape[0]
    h1 = K.shape[0] - rg
    h0_ctrl = nA - rg

    sigma = permutation_facteurs_numerique(coord_mats[0], ambient, p)
    T = np.zeros((nBC, nBC), dtype=np.int64)
    for i in range(n):
        Ti, dimg = matrice_quotient(anneau, coord_mats[0], ambient, sigma,
                                    dBC[i], p)
        if list(dimg) != list(dBC[i]):
            return None
        T[oB[i]:oB[i] + Ti.shape[0], oB[i]:oB[i] + Ti.shape[1]] = Ti

    tK = multiplicites_propres(T, p, valeurs, base=K)
    tI = multiplicites_propres(T, p, valeurs, base=imA)
    coherent = (sum(tK.values()) == K.shape[0] and sum(tI.values()) == rg)
    return {
        'h1': h1, 'h0_controle': h0_ctrl, 'complexe_ok': complexe_ok,
        'dim_VC': K.shape[0], 'rang_alpha': rg,
        'sur_VC': tK, 'sur_image': tI,
        'sur_H1': {k: tK[k] - tI[k] for k in valeurs},
        'coherent': coherent,
    }


def espace_total(anneau, b_charges, c_charges, cases, dims):
    """
    Base de l'espace ENTIER des f, au format attendu par les fonctions
    `..._sur_espace`. Sert de reference : le meme code doit redonner la
    valeur generique de `sections.h0_*_explicit` quand on ne contraint rien.
    """
    N = sum(dims[k] for k in cases)
    return np.eye(N, dtype=np.int64)


def _marge_predite(ambient, b_charges, c, d):
    """
    dim source - dim cible, estimee par `dim_multi` (formule fermee).

    Le critere J_d = R_d ne peut aboutir que si la source couvre la cible :
    le rang est majore par la source. On le PREDIT ici a bon compte, dans
    l'ambiant, pour ne soumettre que des degres viables -- `dimY`
    construirait le quotient, c'est-a-dire une rref, soit l'operation qu'on
    cherche justement a eviter.

    Renvoie None si la cible est vide.
    """
    m = len(c)
    cible = dim_multi(ambient, d)
    if cible == 0:
        return None
    src = 0
    for b in b_charges:
        # deg f_i = c - b_i, donc la source du terme i est R_{d - c + b_i}
        src += dim_multi(ambient, [d[k] - c[k] + b[k] for k in range(m)])
    return src - cible


def _degres_a_essayer(ambient, c, maxdim, n_degres, t_max=8, b_charges=None):
    """
    Multidegres d >= c a soumettre au critere J_d = R_d.

    ----------------------------------------------------------------------
    Ce que la version precedente ne balayait pas -- et le chiffre que cela
    fabriquait
    ----------------------------------------------------------------------
    Elle engendrait deux familles : les marches LONGUES sur un seul axe
    (c + t.e_k) et les pas COURTS mixtes (c + v, |v|_1 <= 2), plus la
    croissance uniforme c + t.(1,...,1). Elle ne combinait jamais longueur
    et mixite. D'ou le constat du §5.4 : « au rang 5 la source ne rattrape
    jamais la cible », et 449 couples classes `indetermine : surjectivite
    de f non certifiee`, dont 420 de rang 5.

    Ce constat mesurait la liste des directions essayees, pas la geometrie.
    Balayage en directions mixtes longues sur #21 (rang 5, m = 5) :

        d = [4, 4, 2, 5, 4]   source 2670   cible 1278   marge +1392
        d = [2, 4, 4, 2, 7]   source 2446   cible 1074   marge +1372
        d = [2, 2, 8, 5, 4]   source 2642   cible 1398   marge +1244

    Trois cibles bien sous le plafond de 6000, et des marges de l'ordre de
    +1300. Le critere etait atteignable ; c'est l'ensemble balaye qui ne
    l'atteignait pas.

    ----------------------------------------------------------------------
    Deux corrections
    ----------------------------------------------------------------------
    (a) MONTEE ANISOTROPE. Aux familles precedentes s'ajoutent
        c + t.(1,...,1) + s.e_k et une recherche locale gloutonne sur la
        marge predite : depuis les meilleurs germes, on essaie d +/- e_k et
        on garde ce qui ameliore, sous le plafond. Cela trouve des optima a
        support plein sans enumerer {0..q}^m, inabordable des que m grandit.

    (b) SELECTION PAR VIABILITE, si `b_charges` est fourni. L'ancienne
        version triait par cout croissant et gardait les `n_degres` PREMIERS
        -- exactement le mauvais sens, la marge s'ameliorant avec la taille :
        elle depensait son budget sur les degres les moins susceptibles
        d'aboutir. On ne garde desormais que les degres a marge predite
        positive, puis les moins chers PARMI EUX.

    Sans `b_charges`, le comportement historique est conserve.
    """
    from itertools import combinations_with_replacement
    m = len(ambient)
    vecteurs = [(0,) * m]
    for k in range(m):
        for t in range(1, t_max + 1):
            v = [0] * m
            v[k] = t
            vecteurs.append(tuple(v))
    for k, l in combinations_with_replacement(range(m), 2):
        v = [0] * m
        v[k] += 1
        v[l] += 1
        vecteurs.append(tuple(v))
    # Croissance UNIFORME. Le long d'un seul axe, la pente de dim R_{a+t.e_k}
    # depend du degre de base a : la cible (base c) croit plus vite que chaque
    # source (base b_i), et la source ne rattrape jamais -- mesure sur #7300,
    # ecart constant de 32 jusqu'a t = 11. En revanche pour d = c + t.(1,...,1)
    # le terme dominant (1/6) sum d_ijk t^3 ne depend PAS du degre de base :
    # source et cible ont le meme terme principal, et la source l'emporte d'un
    # facteur egal au nombre de f_i.
    for t in range(1, t_max + 1):
        vecteurs.append(tuple([t] * m))
    # Uniforme PLUS une anisotropie : c'est la famille qui manquait.
    for t in range(1, t_max + 1):
        for k in range(m):
            for sgn in (1, 2, 3):
                v = [t] * m
                v[k] += sgn
                vecteurs.append(tuple(v))

    def cout(d):
        return dim_multi(ambient, d)

    cands = []
    for v in vecteurs:
        d = [c[k] + v[k] for k in range(m)]
        taille = cout(d)
        if 0 < taille <= maxdim:
            cands.append(tuple(d))
    cands = sorted(set(cands), key=cout)

    if b_charges is not None:
        # (a) recherche locale gloutonne sur la marge predite, depuis les
        #     germes les plus prometteurs.
        def _mg(d):
            v = _marge_predite(ambient, b_charges, c, list(d))
            return -10 ** 9 if v is None else v
        germes = sorted(cands, key=lambda d: -_mg(d))[:6]
        vus = set(cands)
        for g in germes:
            d = list(g)
            mg = _marge_predite(ambient, b_charges, c, d)
            for _ in range(40):
                meilleur = None
                for k in range(m):
                    for pas in (1, -1):
                        e = list(d)
                        e[k] += pas
                        if e[k] < c[k]:
                            continue
                        t = cout(e)
                        if not (0 < t <= maxdim):
                            continue
                        mm = _marge_predite(ambient, b_charges, c, e)
                        if mm is None:
                            continue
                        if mg is None or mm > mg:
                            if meilleur is None or mm > meilleur[0]:
                                meilleur = (mm, e)
                if meilleur is None:
                    break
                mg, d = meilleur
                if tuple(d) not in vus:
                    vus.add(tuple(d))
                    cands.append(tuple(d))
        # (b) selection par VIABILITE : marge predite >= 0 d'abord, puis le
        #     reste. On ne SUPPRIME pas les non viables, on les repousse --
        #     la marge est estimee dans l'AMBIANT alors que le test exact
        #     porte sur R = S/I, et les deux peuvent differer.
        #
        #     Le test est `mg >= 0`, ecrit explicitement : une marge
        #     EXACTEMENT NULLE est le cas des certificats du §5.4 (source =
        #     cible = 24 sur #6890). Un `mg or -1` la rendrait falsy et
        #     l'ecarterait -- soit precisement les degres qui certifient.
        cands = sorted(cands, key=lambda d: (0 if _mg(d) >= 0 else 1, cout(d)))

    res = []
    for d in cands:
        res.append(list(d))
        if len(res) >= n_degres:
            break
    return res


def f_sans_point_base(anneau, b_charges, c_charges, base, offsets, dims,
                      degres, p, rng, n_essais=3, n_degres=4, maxdim=6000):
    """
    CERTIFICAT de surjectivite de f : B -> C, pour rank_C = 1.

    ----------------------------------------------------------------------
    Pourquoi ce test manquait, et pourquoi il porte precisement sur les
    survivants
    ----------------------------------------------------------------------
    Tout le raisonnement suppose que V = ker(f) est un FIBRE, ce qui exige f
    surjective en tout point de Y. Pour un f generique de l'espace entier,
    c'est acquis. Rien ne le garantit sur le SOUS-ESPACE equivariant, qui
    pourrait etre contenu tout entier dans le lieu ou f chute de rang --
    exactement le piege deja rencontre sur h^0(V), ou l'existence d'un f
    equivariant ne disait rien de la stabilite.

    ----------------------------------------------------------------------
    Le critere
    ----------------------------------------------------------------------
    Pour rank_C = 1, f = (f_1, ..., f_n) est surjective si et seulement si
    les f_i n'ont aucun zero commun sur Y. Soit J l'ideal qu'ils engendrent
    dans R = S/I_Y. S'il existe UN multidegre d >= 0 tel que

        J_d = R_d

    alors il n'y a pas de zero commun. Demonstration : soit y un zero commun
    des f_i ; tout element de J_d s'annule en y ; or y est un point de
    P^{n_1} x ... x P^{n_m}, donc dans chaque facteur au moins une coordonnee
    est non nulle, donc pour tout d >= 0 au moins un monome de multidegre d
    est non nul en y, donc R_d ne s'annule pas entierement en y. D'ou
    J_d != R_d. Contradiction.

    Le critere est donc SUFFISANT et ne peut pas donner de faux positif. Un
    echec a un degre donne n'est en revanche pas concluant : on monte en
    degre, et au-dela de `n_degres` on renvoie `certifie = False` avec
    `concluant = False`.

    J_d est l'image de (+)_i R_{d - deg f_i} --(x f_i)--> R_d. Si un
    d - deg f_i a une composante negative, le terme ne contribue pas : on
    teste alors un sous-module de J_d, ce qui ne peut que rendre le critere
    plus difficile a satisfaire, jamais plus permissif.

    ----------------------------------------------------------------------
    Reserve sur GF(p)
    ----------------------------------------------------------------------
    Le rang est calcule modulo p. Un rang PLEIN mod p implique un rang plein
    en caracteristique zero pour un relevement des memes coefficients -- le
    rang ne peut que chuter par specialisation. Le certificat vaut donc pour
    un relevement de ce f-la, pas pour un f arbitraire. C'est la meme reserve
    que celle de l'en-tete de `sections.py` : un resultat NUL obtenu mod p
    est concluant, un resultat non nul ne l'est pas.

    Un seul essai qui certifie suffit : on cherche l'EXISTENCE d'un f
    equivariant surjectif, pas une propriete de tous.

    Renvoie {'certifie', 'concluant', 'degre', 'essais'}.
    """
    from cy_landscape.core.sections import _mult_matrix
    if len(c_charges) != 1:
        return {'certifie': False, 'concluant': False, 'degre': None,
                'motif': 'rank_C >= 2', 'essais': []}
    m = len(b_charges[0])
    c = list(c_charges[0])
    journal = []

    for _ in range(n_essais):
        v = (rng.randint(0, p, size=base.shape[0]) @ base) % p
        # coefficients de chaque f_i, sur la base monomiale complete
        fpol = {}
        for i in range(len(b_charges)):
            if (0, i) not in offsets:
                continue
            deg = degres[0][i]
            S, idx, free, piv, Mred = anneau.quotient(deg)
            coeffs = np.zeros(len(S), dtype=np.int64)
            for t, k in enumerate(free):
                coeffs[k] = v[offsets[(0, i)] + t]
            if coeffs.any():
                fpol[i] = (deg, S, coeffs)
        if not fpol:
            journal.append(('f identiquement nul', None, None))
            continue

        for d in _degres_a_essayer(anneau.amb, c, maxdim, n_degres,
                                   b_charges=b_charges):
            cible = anneau.dimY(d)
            if cible == 0:
                continue
            # La source doit pouvoir couvrir la cible. Sinon le degre est
            # INUTILISABLE : le rang est majore par dim source, l'echec est
            # arithmetique et ne dit rien de la surjectivite. Ne pas le
            # filtrer revenait a compter comme « non certifie » des degres ou
            # le critere n'avait aucune chance -- et a payer la rref pour rien.
            sources = []
            for i, (deg, S, coeffs) in fpol.items():
                s = [d[k] - deg[k] for k in range(m)]
                if any(x < 0 for x in s) or anneau.dimY(s) == 0:
                    continue
                sources.append((i, s, deg, S, coeffs))
            dim_src = sum(anneau.dimY(s) for _, s, _, _, _ in sources)
            if dim_src < cible:
                journal.append((tuple(d), 'source insuffisante',
                                f'{dim_src} < {cible}'))
                continue
            cols = []
            for i, s, deg, S, coeffs in sources:
                blk = _mult_matrix(anneau, s, deg, (S, coeffs), d)
                if blk.size:
                    cols.append(blk)
            if not cols:
                continue
            M = np.hstack(cols)
            rang, _ = rref_complet(M, p)
            rang = rang.shape[0]
            journal.append((tuple(d), rang, cible))
            if rang == cible:
                return {'certifie': True, 'concluant': True,
                        'degre': tuple(d), 'essais': journal}
    return {'certifie': False, 'concluant': False, 'degre': None,
            'essais': journal}
