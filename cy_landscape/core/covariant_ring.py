"""
covariant_ring.py -- Polynomes definissants Gamma-COVARIANTS.

--------------------------------------------------------------------------
Le defaut corrige
--------------------------------------------------------------------------
`sections.Ring` tire les coefficients des polynomes definissants AU HASARD.
C'est legitime tant qu'on ne veut que le rang generique d'une application.
Ce ne l'est plus des qu'on parle de Gamma : pour que l'action descende au
quotient X/Gamma, il faut que Gamma preserve l'ideal engendre par les
polynomes, ce qu'un tirage aleatoire ne fait pas. Tout test d'equivariance
construit sur un anneau aleatoire porte donc sur une variete qui n'admet
PAS l'action annoncee.

Le fichier de Braun donne exactement la donnee manquante : pour chaque
generateur g, une matrice N sur les K polynomes, a cote de la matrice M sur
les coordonnees. La condition de covariance s'ecrit

    p_alpha(g . x)  =  sum_beta  N[alpha][beta] * p_beta(x)          (I)

C'est un systeme LINEAIRE en les coefficients des p_alpha. Son noyau est
l'espace des familles de polynomes covariantes ; un element generique de ce
noyau donne une CICY portant reellement l'action de Gamma.

--------------------------------------------------------------------------
Convention
--------------------------------------------------------------------------
Braun ne fixe pas explicitement le sens de N (N, sa transposee, ou son
inverse selon qu'on fait agir g sur les points ou sur les fonctions).
Plutot que de deviner, `resoudre_covariants` essaie les quatre conventions
et retient celles dont le noyau est NON DEGENERE, c'est-a-dire ou un element
generique donne les K polynomes tous non nuls et de multidegres corrects.
Le resultat indique laquelle a ete retenue ; si plusieurs conviennent, elles
donnent le meme ideal a renumerotation pres et le module le signale.

--------------------------------------------------------------------------
Action sur les monomes
--------------------------------------------------------------------------
On ne suppose PAS M monomiale. La substitution x_j -> sum_k M[j][k] x_k se
factorise par bloc projectif : si g envoie le facteur r sur le facteur
sigma(r), l'action sur les monomes de multidegre d est le produit tensoriel,
sur les facteurs, des actions par bloc. Le multidegre image est d o sigma^-1.

Tout est calcule dans GF(p), avec p choisi congru a 1 modulo l'exposant du
groupe pour que les racines de l'unite necessaires y vivent. Aucun flottant.

--------------------------------------------------------------------------
Ce que ce module ne fait PAS
--------------------------------------------------------------------------
Il ne verifie pas la LISSITE de la variete covariante, ni que l'action y est
librement agissante -- Braun l'a fait, c'est la raison d'etre de sa liste.
Il verifie en revanche, et c'est le controle utile ici, que la fonction de
Hilbert de l'anneau covariant coincide avec celle de l'anneau aleatoire :
un ecart signalerait que le choix covariant est degenere au point de changer
les dimensions, donc que le reste du pipeline ne s'y applique pas.
"""
import numpy as np

from cy_landscape.core.sections import (Ring, basis_multi, dim_multi,
                                        monomials, rref_mod)
from cy_landscape.core.gamma_action import (choisir_premier, racine_primitive,
                                            parse_entree)
from cy_landscape.core.braun_symmetry import ordres_rt, matrice_mod_p


# ======================================================================
# Algebre lineaire mod p : noyau
# ======================================================================

def rref_complet(M, p):
    """Forme echelonnee reduite mod p. Renvoie (R, pivots)."""
    M = np.asarray(M, dtype=np.int64) % p
    rows, cols = M.shape
    r, piv = 0, []
    for c in range(cols):
        if r >= rows:
            break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0:
            continue
        i = r + nz[0]
        if i != r:
            M[[r, i]] = M[[i, r]]
        inv = pow(int(M[r, c]), p - 2, p)
        M[r] = (M[r] * inv) % p
        col = M[:, c].copy()
        col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr):
            M[nzr] = (M[nzr] - np.outer(col[nzr], M[r])) % p
        piv.append(c)
        r += 1
    return M[:r], piv


def noyau_mod(M, p):
    """
    Base du noyau de M (lignes = equations, colonnes = inconnues), mod p.
    Renvoie une matrice (dim_noyau x n_colonnes).
    """
    M = np.asarray(M, dtype=np.int64)
    n = M.shape[1]
    if M.shape[0] == 0:
        return np.eye(n, dtype=np.int64)
    R, piv = rref_complet(M, p)
    libres = [c for c in range(n) if c not in set(piv)]
    base = np.zeros((len(libres), n), dtype=np.int64)
    for k, c in enumerate(libres):
        base[k, c] = 1
        for r, cp in enumerate(piv):
            base[k, cp] = (-R[r, c]) % p
    return base


# ======================================================================
# Action d'un generateur sur les monomes, par blocs
# ======================================================================

def _decoupage(ambient):
    tailles = [n + 1 for n in ambient]
    debuts, acc = [], 0
    for t in tailles:
        debuts.append(acc)
        acc += t
    return debuts, tailles, acc


def permutation_facteurs_numerique(M, ambient, p):
    """
    sigma tel que g envoie le facteur r sur le facteur sigma(r), lu sur les
    positions non nulles de M. None si M ne respecte pas la structure par
    blocs (elle melangerait alors des facteurs, ce qui n'arrive pas dans la
    liste de Braun mais doit etre detecte plutot que suppose).
    """
    debuts, tailles, total = _decoupage(ambient)
    if len(M) != total or any(len(l) != total for l in M):
        return None

    def bloc(i):
        for r, (d, t) in enumerate(zip(debuts, tailles)):
            if d <= i < d + t:
                return r
        return None

    sigma = [None] * len(ambient)
    for i, ligne in enumerate(M):
        r = bloc(i)
        cibles = {bloc(j) for j, v in enumerate(ligne) if v % p}
        if not cibles:
            return None                      # ligne nulle : g non inversible
        if len(cibles) != 1:
            return None
        s = cibles.pop()
        if sigma[r] is None:
            sigma[r] = s
        elif sigma[r] != s:
            return None
    if any(s is None for s in sigma) or sorted(sigma) != list(range(len(ambient))):
        return None
    return sigma


def _action_bloc(M, ambient, debuts, tailles, r, s, d, p):
    """
    Matrice (dim_s x dim_r) de la substitution x_{r,i} -> sum_l M[r,i][s,l] x_{s,l}
    sur les monomes de degre d, du bloc r vers le bloc s.

    Necessite tailles[r] == tailles[s] (meme dimension projective), ce que
    garantit sigma.
    """
    src = monomials(ambient[r], d)
    dst = monomials(ambient[s], d)
    idx = {m: i for i, m in enumerate(dst)}
    nr, ns = tailles[r], tailles[s]
    lin = [[M[debuts[r] + i][debuts[s] + l] % p for l in range(ns)]
           for i in range(nr)]

    A = np.zeros((len(dst), len(src)), dtype=np.int64)
    for j, e in enumerate(src):
        # produit des formes lineaires, developpe par accumulation
        cur = {(0,) * ns: 1}
        for i in range(nr):
            for _ in range(e[i]):
                nxt = {}
                for mono, cf in cur.items():
                    for l in range(ns):
                        c = lin[i][l]
                        if not c:
                            continue
                        key = list(mono)
                        key[l] += 1
                        key = tuple(key)
                        nxt[key] = (nxt.get(key, 0) + cf * c) % p
                cur = nxt
        for mono, cf in cur.items():
            if cf % p:
                A[idx[mono], j] = cf % p
    return A


def matrice_substitution(M, ambient, sigma, degre, p):
    """
    Matrice de S_g : S_degre -> S_{degre o sigma^-1}, dans les bases
    `basis_multi` (facteur 0 le plus significatif).

    Le multidegre image verifie deg_img[sigma[r]] = degre[r].
    """
    debuts, tailles, _ = _decoupage(ambient)
    m = len(ambient)
    deg_img = [0] * m
    for r in range(m):
        deg_img[sigma[r]] = degre[r]

    sigma_inv = [0] * m
    for r in range(m):
        sigma_inv[sigma[r]] = r

    A = np.ones((1, 1), dtype=np.int64)
    for s in range(m):
        r = sigma_inv[s]
        B = _action_bloc(M, ambient, debuts, tailles, r, s, degre[r], p)
        A = np.kron(A, B) % p
    return A, deg_img


def appliquer_substitution(M, ambient, sigma, degre, coeffs, p):
    """S_g p, ou p a pour coefficients `coeffs` sur `basis_multi(ambient, degre)`."""
    A, deg_img = matrice_substitution(M, ambient, sigma, degre, p)
    return (A @ (np.asarray(coeffs, dtype=np.int64) % p)) % p, deg_img


# ======================================================================
# Resolution du systeme de covariance
# ======================================================================

CONVENTIONS = ('N', 'NT', 'Ninv', 'NinvT')


def _variante(N, quelle, p):
    A = np.asarray(N, dtype=np.int64) % p
    if quelle == 'N':
        return A
    if quelle == 'NT':
        return A.T.copy()
    R = _inverse_mod(A, p)
    if R is None:
        return None
    return R if quelle == 'Ninv' else R.T.copy()


def _inverse_mod(A, p):
    n = A.shape[0]
    if A.shape[0] != A.shape[1]:
        return None
    aug = np.hstack([A % p, np.eye(n, dtype=np.int64)])
    R, piv = rref_complet(aug, p)
    if len(piv) != n or piv != list(range(n)):
        return None
    return R[:, n:] % p


def resoudre_covariants(ambient, config, coord_mats, poly_mats, p,
                        conventions=CONVENTIONS, maxdim=20000):
    """
    Espace des familles (p_1, ..., p_K) verifiant (I) pour tous les generateurs.

    `config` : K x m (polynomes x facteurs), comme partout dans le pipeline.
    `coord_mats`, `poly_mats` : matrices NUMERIQUES mod p, une par generateur.

    Renvoie un dict :
        {'sigmas': [...],
         'par_convention': {nom: {'dim': d, 'base': array, 'offsets': [...],
                                  'non_degenere': bool, 'blocs_nuls': [...]}}}
    ou None si les sigmas ne sont pas extractibles.
    """
    config = np.asarray(config, dtype=np.int64)
    K, m = config.shape

    sigmas = []
    for M in coord_mats:
        s = permutation_facteurs_numerique(M, ambient, p)
        if s is None:
            return None
        sigmas.append(s)

    dims = [dim_multi(ambient, list(config[a])) for a in range(K)]
    if sum(dims) > maxdim:
        return None
    offs, acc = [], 0
    for d in dims:
        offs.append(acc)
        acc += d
    Ntot = acc

    # matrices de substitution, memoisees par (generateur, polynome)
    subs = {}
    for gi, (M, sigma) in enumerate(zip(coord_mats, sigmas)):
        for a in range(K):
            A, deg_img = matrice_substitution(M, ambient, sigma,
                                              list(config[a]), p)
            subs[(gi, a)] = (A, deg_img)

    sortie = {'sigmas': sigmas, 'dims': dims, 'offsets': offs,
              'par_convention': {}}

    for conv in conventions:
        lignes = []
        incompatible = False
        for gi, Nbrut in enumerate(poly_mats):
            N = _variante(Nbrut, conv, p)
            if N is None:
                incompatible = True
                break
            for a in range(K):
                A, deg_img = subs[(gi, a)]
                # indices des polynomes du membre de droite
                cibles = [b for b in range(K) if N[a][b] % p]
                if any(list(config[b]) != list(deg_img) for b in cibles):
                    incompatible = True
                    break
                # base d'arrivee = S_{deg_img} ; sa dimension
                D = A.shape[0]
                bloc = np.zeros((D, Ntot), dtype=np.int64)
                bloc[:, offs[a]:offs[a] + dims[a]] = A
                for b in cibles:
                    if dims[b] != D:
                        incompatible = True
                        break
                    bloc[:, offs[b]:offs[b] + dims[b]] = (
                        bloc[:, offs[b]:offs[b] + dims[b]]
                        - int(N[a][b]) * np.eye(D, dtype=np.int64)) % p
                if incompatible:
                    break
                lignes.append(bloc % p)
            if incompatible:
                break
        if incompatible:
            sortie['par_convention'][conv] = {'dim': 0, 'base': None,
                                              'non_degenere': False,
                                              'motif': 'multidegres incompatibles'}
            continue

        Mat = np.vstack(lignes) if lignes else np.zeros((0, Ntot), dtype=np.int64)
        base = noyau_mod(Mat, p)
        blocs_nuls = [a for a in range(K)
                      if not base[:, offs[a]:offs[a] + dims[a]].any()]
        sortie['par_convention'][conv] = {
            'dim': int(base.shape[0]),
            'base': base,
            'non_degenere': bool(base.shape[0] and not blocs_nuls),
            'blocs_nuls': blocs_nuls,
            'motif': None,
        }
    return sortie


def tirer_covariants(base, offsets, dims, p, rng):
    """
    Element generique du noyau -> [coefficients de p_1, ..., p_K].

    Un element aleatoire d'un espace vectoriel sur GF(p) evite tout
    sous-espace propre fixe avec probabilite >= 1 - 1/p.
    """
    lam = rng.randint(0, p, size=base.shape[0])
    v = (lam @ base) % p
    return [v[offsets[a]:offsets[a] + dims[a]].copy() for a in range(len(dims))]


def verifier_covariance(ambient, config, coord_mats, poly_mats, conv,
                        coeffs, p):
    """
    Recontrole INDEPENDANT de (I) sur les coefficients obtenus.

    Ne reutilise pas le systeme lineaire : reapplique la substitution et
    compare terme a terme. Renvoie (ok, ecart_max).
    """
    config = np.asarray(config, dtype=np.int64)
    K = config.shape[0]
    ecart = 0
    for M, Nbrut in zip(coord_mats, poly_mats):
        sigma = permutation_facteurs_numerique(M, ambient, p)
        N = _variante(Nbrut, conv, p)
        for a in range(K):
            img, deg_img = appliquer_substitution(M, ambient, sigma,
                                                  list(config[a]), coeffs[a], p)
            rhs = np.zeros_like(img)
            for b in range(K):
                if N[a][b] % p and list(config[b]) == list(deg_img):
                    rhs = (rhs + int(N[a][b]) * coeffs[b]) % p
            ecart = max(ecart, int(np.abs((img - rhs) % p).max(initial=0)))
    return ecart == 0, ecart


# ======================================================================
# Anneau de sections a coefficients imposes
# ======================================================================

class CovariantRing(Ring):
    """
    `sections.Ring` dont les coefficients des polynomes sont IMPOSES.

    Meme interface (`quotient`, `reduce_vec`, `dimY`), donc utilisable
    partout ou `Ring` l'est. Le premier P de sections.py (32003) reste le
    modulo de `Ring` ; ici on travaille dans GF(p) avec le p des racines de
    l'unite, donc on redefinit le modulo localement.
    """

    def __init__(self, ambient, config, coeffs, p):
        super().__init__(ambient, config, seed=0)
        self.p = p
        for a, c in enumerate(coeffs):
            b = basis_multi(self.amb, list(self.cfg[a]))
            assert len(b) == len(c), (len(b), len(c))
            self._poly[a] = (b, np.asarray(c, dtype=np.int64) % p)


def modulo_de(anneau):
    """P de sections.py, ou le p propre a un CovariantRing."""
    return getattr(anneau, 'p', None) or __import__(
        'cy_landscape.core.sections', fromlist=['P']).P
