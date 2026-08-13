"""
gamma_action.py -- Action de Gamma sur les sections, et equivariance de f.

--------------------------------------------------------------------------
Format des generateurs de Braun
--------------------------------------------------------------------------
Les entrees des matrices sont 0, 1, -1, ou des puissances de rt[n], ou
rt[n] designe exp(2 i pi / n). Exemple reel (CICY 7669, P2xP2xP2 avec trois
hypersurfaces de tridegre (1,1,1)) :

    "Z3 x Z3", { <matrice 9x9 diagonale en rt[3]^k>,
                 <matrice 9x9 de permutation des trois facteurs> },
               { <matrice 3x3 sur les polynomes>, <matrice 3x3> }

Le premier bloc agit sur les 9 coordonnees, le second sur les 3 polynomes
definissants.

--------------------------------------------------------------------------
Corps de travail
--------------------------------------------------------------------------
Plutot que de manipuler des racines de l'unite complexes, on choisit un
premier p tel que p = 1 mod n : GF(p) contient alors une racine primitive
n-ieme, et tout le calcul reste de l'algebre lineaire modulaire, dans la
continuite de sections.py. Aucune approximation flottante.

--------------------------------------------------------------------------
LIMITE CONNUE : le relevement n'est que PROJECTIF
--------------------------------------------------------------------------
Gamma agit sur X, mais son relevement aux fibres en droites n'est defini
qu'a une phase pres. Concretement, sur la CICY 7669, les deux generateurs
donnes par Braun sont une matrice de phases et une permutation cyclique des
coordonnees a l'interieur de chaque facteur. Ces deux matrices NE COMMUTENT
PAS : elles engendrent le groupe de Heisenberg d'ordre 27, et ce n'est que
sur la variete -- donc projectivement -- que le quotient est Z3 x Z3.

Consequence mesuree, par decomposition isotypique sur l'espace ambiant
(la somme des dimensions propres doit valoir la dimension totale) :

    degre (2,1,0) : 18 sur 18   -> coherent
    degre (1,1,2) : 54 sur 54   -> coherent
    degre (1,1,1) : 0 sur 27    -> AUCUN vecteur propre commun

Le degre (1,1,1) est precisement celui des polynomes definissants. Sur ces
degres le cocycle du relevement est non trivial : il n'existe pas de base de
vecteurs propres communs, et chercher f dans un espace propre commun -- ce
que fait `espace_f_equivariant` -- n'a pas de sens.

Traiter correctement ce cas demande de travailler avec l'extension centrale
(le groupe de Heisenberg) plutot qu'avec Gamma, et de choisir un relevement
compatible sur B et sur C. C'est le travail qui reste. En l'etat,
`espace_f_equivariant` n'est valide que sur les degres ou le cocycle est
trivial, et le controle de coherence ci-dessus permet de le detecter.

De plus, `Ring` construit les polynomes definissants avec des coefficients
ALEATOIRES, alors qu'ils doivent etre Gamma-covariants -- l'ideal n'est
sinon pas preserve et l'action ne descend pas au quotient. Le second bloc de
generateurs du fichier de Braun, celui qui agit sur les polynomes, sert
exactement a cela et n'est pas encore exploite.

--------------------------------------------------------------------------
Condition d'equivariance testee
--------------------------------------------------------------------------
Pour V = ker(f : B -> C) avec B = (+) O(b_i) et C = O(c) :

  - g envoie O(b) sur O(b o sigma), ou sigma est la permutation des facteurs.
    Il doit donc exister une permutation pi des facteurs de B telle que
    b_{pi(i)} = b_i o sigma. C'est la condition sur les CHARGES, deja testee
    par equivariance.py.
  - f = (f_i), f_i dans H^0(O(c - b_i)). L'equivariance impose
        g . f_{pi(i)} = chi(g) * f_i
    pour un caractere chi de Gamma. Chercher un tel f revient a chercher un
    vecteur propre commun aux operateurs induits par les generateurs, dans
    (+)_i H^0(Y, O(c - b_i)).

C'est cette derniere condition -- la seule qui porte sur les polynomes --
que ce module resout. Un espace propre non nul signifie qu'une monade
equivariante de ces charges existe ; il reste alors a verifier qu'un element
generique de cet espace donne bien un fibre (f surjective), ce que fait
`existe_monade_equivariante`.
"""
import re
from itertools import product as iprod

import numpy as np

from cy_landscape.core.sections import (basis_multi, rref_mod, Ring)


# ======================================================================
# Corps GF(p) contenant les racines de l'unite necessaires
# ======================================================================

def choisir_premier(ordres, minimum=10007):
    """Plus petit premier > minimum tel que p = 1 mod n pour tout n."""
    from math import gcd
    n = 1
    for o in ordres:
        n = n * o // gcd(n, o)
    p = minimum + (1 - minimum) % n
    while True:
        p += n
        if _est_premier(p):
            return p, n


def _est_premier(x):
    if x < 2:
        return False
    for d in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        if x % d == 0:
            return x == d
    d, r = x - 1, 0
    while d % 2 == 0:
        d //= 2; r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        v = pow(a, d, x)
        if v in (1, x - 1):
            continue
        for _ in range(r - 1):
            v = v * v % x
            if v == x - 1:
                break
        else:
            return False
    return True


def racine_primitive(p, n):
    """Racine primitive n-ieme de l'unite dans GF(p) (p = 1 mod n)."""
    assert (p - 1) % n == 0
    for g in range(2, p):
        r = pow(g, (p - 1) // n, p)
        if r != 1 and pow(r, n, p) == 1:
            # verifier la primitivite
            if all(pow(r, n // q, p) != 1 for q in _facteurs(n)):
                return r
    raise RuntimeError("racine primitive introuvable")


def _facteurs(n):
    f, d = set(), 2
    while d * d <= n:
        while n % d == 0:
            f.add(d); n //= d
        d += 1
    if n > 1:
        f.add(n)
    return f


# ======================================================================
# Lecture des entrees rt[n]^k
# ======================================================================

def parse_entree(txt, p, racines):
    """
    Convertit une entree Mathematica en element de GF(p).
    Gere 0, 1, -1, rt[n], rt[n]^k, et les signes.
    `racines` : {n: racine primitive n-ieme dans GF(p)}.
    """
    t = txt.strip().replace(' ', '')
    if not t:
        return 0
    signe = 1
    while t.startswith('-'):
        signe = -signe
        t = t[1:]
    if t in ('', '1'):
        return signe % p
    if t == '0':
        return 0
    m = re.fullmatch(r'rt\[(\d+)\](?:\^(\d+))?', t)
    if m:
        n = int(m.group(1))
        k = int(m.group(2) or 1)
        if n not in racines:
            raise ValueError(f"racine {n}-ieme non disponible")
        return (signe * pow(racines[n], k, p)) % p
    try:
        return (signe * int(t)) % p
    except ValueError:
        raise ValueError(f"entree non reconnue : {txt!r}")


def ordres_requis(txt):
    """Ensemble des n apparaissant dans les rt[n] d'un bloc."""
    return {int(x) for x in re.findall(r'rt\[(\d+)\]', txt)}


# ======================================================================
# Action sur les monomes
# ======================================================================

def action_monomes(M, ambient, degre, p):
    """
    Action de g sur H^0(A, O(degre)), pour M matrice MONOMIALE.

    La substitution x_j -> m_j x_tau(j) envoie un monome sur un multiple d'un
    autre monome, mais PAS forcement de meme multidegre : si g permute les
    facteurs projectifs, le degre part sur degre o sigma^-1. La version
    precedente cherchait l'image dans la base du meme degre, ce qui n'etait
    correct que pour un degre invariant -- une hypothese fausse des que
    sigma est non trivial.

    Renvoie (base_src, degre_img, base_img, perm, coeffs) :
        base_src[i]  ->  coeffs[i] * base_img[perm[i]]
    """
    tailles = [n + 1 for n in ambient]
    debuts, acc = [], 0
    for t in tailles:
        debuts.append(acc); acc += t
    total = acc

    tau = [None] * total
    mult = [0] * total
    for j in range(total):
        nz = [(k, M[j][k]) for k in range(total) if M[j][k] % p]
        if len(nz) != 1:
            return None
        tau[j], mult[j] = nz[0]

    # permutation des facteurs induite par tau
    sigma = [None] * len(ambient)
    for r in range(len(ambient)):
        cibles = set()
        for k in range(tailles[r]):
            j = debuts[r] + k
            for rr in range(len(ambient)):
                if debuts[rr] <= tau[j] < debuts[rr] + tailles[rr]:
                    cibles.add(rr)
        if len(cibles) != 1:
            return None
        sigma[r] = cibles.pop()
    if sorted(sigma) != list(range(len(ambient))):
        return None

    degre_img = [0] * len(ambient)
    for r in range(len(ambient)):
        degre_img[sigma[r]] = degre[r]

    base_src = basis_multi(ambient, list(degre))
    base_img = basis_multi(ambient, list(degre_img))
    idx_img = {m: i for i, m in enumerate(base_img)}

    perm = [0] * len(base_src)
    coeffs = [1] * len(base_src)
    for i, mono in enumerate(base_src):
        expo = [0] * total
        for r, per in enumerate(mono):
            for k, e in enumerate(per):
                expo[debuts[r] + k] = e
        cf = 1
        img = [0] * total
        for j in range(total):
            if expo[j]:
                cf = cf * pow(mult[j], expo[j], p) % p
                img[tau[j]] += expo[j]
        cible = tuple(tuple(img[debuts[r]:debuts[r] + tailles[r]])
                      for r in range(len(ambient)))
        k = idx_img.get(cible)
        if k is None:
            return None
        perm[i] = k
        coeffs[i] = cf
    return base_src, degre_img, base_img, perm, coeffs


def permutation_charges(charges, sigma):
    """
    Permutation pi telle que charges[pi[i]] = charges[i] o sigma^-1,
    c'est-a-dire l'appariement des facteurs induit par g. None si absente.
    """
    m = len(sigma)
    img = []
    for ch in charges:
        v = [0] * m
        for r in range(m):
            v[sigma[r]] = ch[r]
        img.append(tuple(v))
    reste = {i: tuple(ch) for i, ch in enumerate(charges)}
    pi = [None] * len(charges)
    for i, t in enumerate(img):
        for j, v in list(reste.items()):
            if v == t:
                pi[i] = j
                del reste[j]
                break
        else:
            return None
    return pi


def espace_f_equivariant(anneau, ambient, b_charges, c_charges,
                         generateurs, p, caracteres=None, maxdim=3000):
    """
    Dimension du sous-espace des f = (f_ji) equivariantes.

    f est la matrice rank_C x rank_B dont l'entree (j,i) vit dans
    H^0(Y, O(c_j - b_i)) -- nulle des qu'une composante du degre est
    negative. L'equivariance sous un generateur g s'ecrit

        g . f_{rho(j), pi(i)} = chi(g) * f_{j,i}

    ou pi et rho sont les permutations induites sur les facteurs de B et de C.

    `generateurs` : liste de (matrice_coordonnees, sigma).
    `caracteres`  : liste de listes chi(g) a essayer ; None -> caractere
                    trivial seulement. Un f equivariant peut n'exister que
                    pour un caractere non trivial, donc les balayer tous
                    est necessaire pour conclure a l'absence.

    Renvoie {caractere_index: dimension} ou None si hors domaine.
    """
    nB, nC = len(b_charges), len(c_charges)
    m = len(b_charges[0])

    degres = [[[c_charges[j][k] - b_charges[i][k] for k in range(m)]
               for i in range(nB)] for j in range(nC)]
    actif = [[all(x >= 0 for x in degres[j][i]) for i in range(nB)]
             for j in range(nC)]

    # coordonnees : une case par entree active
    cases = [(j, i) for j in range(nC) for i in range(nB) if actif[j][i]]
    if not cases:
        return None
    dims, offs, acc = {}, {}, 0
    for (j, i) in cases:
        d = anneau.dimY(degres[j][i])
        dims[(j, i)] = d
        offs[(j, i)] = acc
        acc += d
    N = acc
    if N == 0 or N > maxdim:
        return None

    # permutations induites
    paires = []
    for M, sigma in generateurs:
        pi = permutation_charges(b_charges, sigma)
        rho = permutation_charges(c_charges, sigma)
        if pi is None or rho is None:
            return None
        paires.append((M, sigma, pi, rho))

    if caracteres is None:
        caracteres = [[1] * len(generateurs)]

    resultats = {}
    for ic, chi in enumerate(caracteres):
        lignes = []
        for gi, (M, sigma, pi, rho) in enumerate(paires):
            for (j, i) in cases:
                jj, ii = rho[j], pi[i]
                if not actif[jj][ii]:
                    # f_{jj,ii} est nul : la contrainte impose f_{j,i} = 0
                    for t in range(dims[(j, i)]):
                        v = np.zeros(N, dtype=np.int64)
                        v[offs[(j, i)] + t] = 1
                        lignes.append(v)
                    continue
                act = action_monomes(M, ambient, degres[jj][ii], p)
                if act is None:
                    return None
                base_src, deg_img, base_img, perm, coeffs = act
                if list(deg_img) != list(degres[j][i]):
                    return None
                S_src, idx_src, free_src, _, _ = anneau.quotient(degres[jj][ii])
                S_img, idx_img, free_img, _, _ = anneau.quotient(degres[j][i])
                for col, pos in enumerate(free_src):
                    v = np.zeros(len(S_img), dtype=np.int64)
                    v[perm[pos]] = coeffs[pos] % p
                    red = anneau.reduce_vec(degres[j][i], v)
                    ligne = np.zeros(N, dtype=np.int64)
                    for t, val in enumerate(red):
                        if val:
                            ligne[offs[(j, i)] + t] = int(val) % p
                    ligne[offs[(jj, ii)] + col] = (
                        ligne[offs[(jj, ii)] + col] - chi[gi]) % p
                    lignes.append(ligne)
        if not lignes:
            resultats[ic] = N
            continue
        A = np.array(lignes, dtype=np.int64) % p
        rang, _ = rref_mod(A, p)
        resultats[ic] = N - rang
    return resultats


def caracteres_abeliens(ordres, p, racines):
    """Tous les caracteres d'un produit de cycliques, dans GF(p)."""
    listes = [[pow(racines[o], k, p) for k in range(o)] for o in ordres]
    return [list(t) for t in iprod(*listes)]
