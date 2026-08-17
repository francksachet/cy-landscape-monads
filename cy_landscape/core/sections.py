"""
sections.py -- Anneaux de sections et RANGS REELS des applications induites.

--------------------------------------------------------------------------
Ce que ce module debloque
--------------------------------------------------------------------------
Tout le reste du pipeline ne manipule que des DIMENSIONS de groupes de
cohomologie. C'est suffisant pour chi, et pour des bornes, mais pas pour
h^0(wedge^2 V) : cette quantite est le noyau d'une application explicite, et
un noyau ne se lit pas sur les dimensions de la source et du but.

Consequence mesuree jusqu'ici : les monades de rang 4 et 5 sortaient en
« non testable » faute de h^0(wedge^2 V), d'ou un catalogue a 402 E6 pour
16 SO(10) et zero SU(5).

Ce module descend d'un cran : il construit les matrices, et calcule leur
rang.

--------------------------------------------------------------------------
Methode
--------------------------------------------------------------------------
Pour V = ker(f : B -> C) avec C de rang 1, wedge^2 V est le noyau de

    wedge^2 B --> B (x) C,     b1 ^ b2 |-> b1 (x) f(b2) - b2 (x) f(b1)

d'ou, en prenant les sections globales :

    h^0(wedge^2 V) = dim ker( H^0(wedge^2 B) --> H^0(B (x) C) )

Les deux membres sont des sommes de fibres en droites. On represente
H^0(Y, O(a)) comme le quotient S_a / I_a, ou S_a est l'espace des monomes
de multidegre a sur l'espace ambiant et I_a le sous-espace engendre par les
polynomes definissants. L'application est alors une matrice explicite, dont
le rang se calcule par elimination de Gauss.

Les calculs se font modulo un grand nombre premier (32003) avec des
coefficients tires au hasard pour les polynomes definissants et pour f.
C'est la formalisation de « monade generique » : le rang obtenu est le rang
generique avec une probabilite tres proche de 1, et il ne peut qu'etre
inferieur ou egal au rang vrai -- donc h^0 calcule est une borne SUPERIEURE
du h^0 generique. Un h^0 nul obtenu ainsi est donc concluant.

--------------------------------------------------------------------------
Domaine de validite -- A RESPECTER
--------------------------------------------------------------------------
Le modele S_a / I_a ne represente H^0(Y, O(a)) que si la suite de Koszul ne
contribue pas au-dela du terme p = 1 en degre 0. `domaine_valide` verifie
deux conditions suffisantes : toutes les charges en jeu positives ou nulles,
et h^0 certifie par `koszul_cohomology_ex` pour chacune d'elles.

Validation, sur le domaine valide :
  h^0(V) explicite contre le calcul de Koszul independant : 174/177 (98 %).
  Les 3 ecarts sont des erreurs de Koszul, pas du module : dans ces cas la
  source est de dimension 10 pour un but de dimension 9, donc le noyau est
  necessairement non nul, alors que Koszul renvoyait 0.

  h^0(wedge^2 V) explicite, rang 4 : 38/38 a l'interieur des bornes
  rigoureuses de `monad_wedge`.

Hors du domaine valide, ce module ne doit pas etre utilise : le quotient
S_a / I_a n'est alors pas H^0(Y, O(a)).

--------------------------------------------------------------------------
Cout
--------------------------------------------------------------------------
38 ms par monade en moyenne (rang 4, max_charge 2), matrices de taille
mediane 43 x 122 et jusqu'a 888. Le cout croit vite avec les charges :
`maxdim` plafonne la taille et fait renvoyer None au-dela, auquel cas
l'appelant retombe sur les bornes rigoureuses.
"""
import numpy as np, itertools
from math import comb
P = 32003

def monomials(n, d):
    """Monomes de degre d en n+1 variables (P^n) -> tuples d'exposants."""
    if d < 0: return []
    def rec(k, rem):
        if k == 1: yield (rem,); return
        for e in range(rem+1):
            for r in rec(k-1, rem-e): yield (e,)+r
    return list(rec(n+1, d))

def basis_multi(ambient, a):
    """Base monomiale de H^0(A, O(a)) sur A = prod P^{n_r}."""
    if any(x < 0 for x in a): return []
    per = [monomials(ambient[r], a[r]) for r in range(len(ambient))]
    return list(itertools.product(*per))

def dim_multi(ambient, a):
    if any(x < 0 for x in a): return 0
    d = 1
    for r, x in enumerate(a): d *= comb(x + ambient[r], ambient[r])
    return d

def mult_index(base_map, mono_a, mono_b):
    key = tuple(tuple(np.add(mono_a[r], mono_b[r])) for r in range(len(mono_a)))
    return base_map.get(key)

def rref_mod_full(M, p=P):
    """
    Forme echelonnee REDUITE mod p ; renvoie (R, pivots) avec R de rang plein.

    `rref_mod` ne renvoyait que (rang, pivots) et travaillait sur une COPIE
    (`M = M % p` cree un nouveau tableau) : la matrice de l'appelant restait
    donc brute. `Ring.quotient` faisait ensuite `Mred = M[:rank]`, c'est-a-dire
    gardait les `rank` premieres lignes NON REDUITES, et `reduce_vec`
    soustrayait `v[c] * Mred[r]` en supposant un pivot egal a 1. Mesure :
    30/30 elements de l'ideal renvoyaient un reste non nul -- l'application
    n'etait pas la projection sur le quotient. Voir tests_regression.
    """
    M = np.asarray(M, dtype=np.int64) % p
    rows, cols = M.shape
    r = 0; piv = []
    for c in range(cols):
        if r >= rows: break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0: continue
        i = r + nz[0]
        if i != r: M[[r, i]] = M[[i, r]]
        inv = pow(int(M[r, c]), p-2, p)
        M[r] = (M[r] * inv) % p
        col = M[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr):
            M[nzr] = (M[nzr] - np.outer(col[nzr], M[r])) % p
        piv.append(c); r += 1
    return M[:r], piv


def rref_mod(M, p=P):
    """Reduction de Gauss mod p ; renvoie (rang, pivots)."""
    M = M % p
    rows, cols = M.shape
    r = 0; piv = []
    for c in range(cols):
        if r >= rows: break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0: continue
        i = r + nz[0]
        if i != r: M[[r, i]] = M[[i, r]]
        inv = pow(int(M[r, c]), p-2, p)
        M[r] = (M[r] * inv) % p
        col = M[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr):
            M[nzr] = (M[nzr] - np.outer(col[nzr], M[r])) % p
        piv.append(c); r += 1
    return r, piv

class Ring:
    """
    R_a = S_a / I_a pour Y = CICY, coefficients generiques mod p.

    `p` est un parametre et non plus la constante P : `covariant_ring` a
    besoin d'un premier congru a 1 modulo l'exposant de Gamma, pour que les
    racines de l'unite du relevement vivent dans le corps. Tant que p n'est
    pas donne, le comportement est strictement l'ancien (p = P = 32003).
    """
    def __init__(self, ambient, config, seed=0, p=P):
        self.amb = ambient; self.cfg = np.asarray(config)
        self.K = self.cfg.shape[0]
        self.p = p
        self.rng = np.random.RandomState(seed)
        self._poly = {}
        self._quot = {}
    def poly(self, alpha):
        """Coefficients generiques du polynome definissant alpha."""
        if alpha not in self._poly:
            b = basis_multi(self.amb, list(self.cfg[alpha]))
            self._poly[alpha] = (b, self.rng.randint(1, self.p, size=len(b)))
        return self._poly[alpha]
    def quotient(self, a):
        """(base de S_a, indices des monomes libres, matrice de reduction)."""
        key = tuple(a)
        if key in self._quot: return self._quot[key]
        p = self.p
        S = basis_multi(self.amb, list(a))
        idx = {m: i for i, m in enumerate(S)}
        rows = []
        for al in range(self.K):
            q = list(self.cfg[al])
            sub = basis_multi(self.amb, [a[r]-q[r] for r in range(len(a))])
            if not sub: continue
            pb, pc = self.poly(al)
            for s in sub:
                v = np.zeros(len(S), dtype=np.int64)
                for m, cf in zip(pb, pc):
                    j = idx.get(tuple(tuple(np.add(m[r], s[r])) for r in range(len(a))))
                    if j is not None: v[j] = (v[j] + cf) % p
                rows.append(v)
        if rows:
            Mred, piv = rref_mod_full(np.array(rows, dtype=np.int64), p)
        else:
            Mred = np.zeros((0, len(S)), dtype=np.int64); piv = []
        free = [i for i in range(len(S)) if i not in set(piv)]
        self._quot[key] = (S, idx, free, piv, Mred)
        return self._quot[key]
    def reduce_vec(self, a, v):
        """Reduit v (dans S_a) modulo I_a -> coordonnees sur les monomes libres."""
        S, idx, free, piv, Mred = self.quotient(a)
        p = self.p
        v = v.copy() % p
        for r, c in enumerate(piv):
            if v[c]:
                v = (v - v[c] * Mred[r]) % p
        return v[free]
    def dimY(self, a):
        S, idx, free, piv, Mred = self.quotient(a)
        return len(free)

def _mult_matrix(R, a_src, a_f, fpoly, a_dst):
    """Matrice de la multiplication R_{a_src} --(. f)--> R_{a_dst}."""
    S_src, idx_src, free_src, _, _ = R.quotient(a_src)
    S_dst, idx_dst, free_dst, _, _ = R.quotient(a_dst)
    fb, fc = fpoly
    cols = []
    for i in free_src:
        m = S_src[i]
        v = np.zeros(len(S_dst), dtype=np.int64)
        for fm, cf in zip(fb, fc):
            key = tuple(tuple(np.add(m[r], fm[r])) for r in range(len(a_src)))
            j = idx_dst.get(key)
            if j is not None: v[j] = (v[j] + cf) % P
        cols.append(R.reduce_vec(a_dst, v))
    if not cols: return np.zeros((len(free_dst), 0), dtype=np.int64)
    return np.array(cols, dtype=np.int64).T

def h0_V_explicit(R, b_charges, c_charges, maxdim=4000):
    """h0(V) = dim ker( H0(B) -> H0(C) ), rangs reels mod p."""
    m = len(b_charges[0]); amb = R.amb
    assert len(c_charges) == 1
    c = c_charges[0]
    dsrc = sum(R.dimY(b) for b in b_charges)
    ddst = R.dimY(c)
    if dsrc > maxdim or ddst > maxdim: return None, (dsrc, ddst)
    blocks = []
    for b in b_charges:
        deg = [c[k]-b[k] for k in range(m)]
        if any(x < 0 for x in deg):
            blocks.append(np.zeros((ddst, R.dimY(b)), dtype=np.int64)); continue
        fb = basis_multi(amb, deg)
        fc = R.rng.randint(1, P, size=len(fb))
        blocks.append(_mult_matrix(R, b, deg, (fb, fc), c))
    M = np.hstack(blocks) if blocks else np.zeros((ddst,0),dtype=np.int64)
    rank, _ = rref_mod(M.T.copy())
    return dsrc - rank, (dsrc, ddst)

def h0_wedge2_V_explicit(R, b_charges, c_charges, maxdim=4000):
    """h0(w2V) = dim ker( H0(w2B) -> H0(B(x)C) )."""
    m = len(b_charges[0]); amb = R.amb
    assert len(c_charges) == 1
    c = c_charges[0]; n = len(b_charges)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    dims_src = {(i,j): R.dimY([b_charges[i][k]+b_charges[j][k] for k in range(m)]) for i,j in pairs}
    dsrc = sum(dims_src.values())
    dst = [[b_charges[k][t]+c[t] for t in range(m)] for k in range(n)]
    ddst = sum(R.dimY(d) for d in dst)
    if dsrc > maxdim or ddst > maxdim: return None, (dsrc, ddst)
    fpolys = {}
    for j in range(n):
        deg = [c[t]-b_charges[j][t] for t in range(m)]
        fb = basis_multi(amb, deg) if all(x>=0 for x in deg) else []
        fpolys[j] = (deg, (fb, R.rng.randint(1,P,size=len(fb)) if fb else np.zeros(0,dtype=np.int64)))
    offs_d = {}; o=0
    for k in range(n): offs_d[k]=o; o+=R.dimY(dst[k])
    M = np.zeros((ddst, dsrc), dtype=np.int64); oc=0
    for (i,j) in pairs:
        src = [b_charges[i][t]+b_charges[j][t] for t in range(m)]
        w = dims_src[(i,j)]
        for (k,l,sgn) in ((i,j,1),(j,i,-1)):
            deg, fp = fpolys[l]
            if not fp[0]: continue
            blk = _mult_matrix(R, src, deg, fp, dst[k])
            if blk.size:
                M[offs_d[k]:offs_d[k]+blk.shape[0], oc:oc+w] = (M[offs_d[k]:offs_d[k]+blk.shape[0], oc:oc+w] + sgn*blk) % P
        oc += w
    rank, _ = rref_mod(M.T.copy())
    return dsrc - rank, (dsrc, ddst)


def domaine_valide(ambient, config, b_charges, c_charges, rank_c_max=1):
    """
    Conditions suffisantes pour que S_a / I_a represente bien H^0(Y, O(a))
    sur toutes les charges intervenant dans le calcul. Voir l'en-tete.

    `rank_c_max` plafonne le rang de C. La valeur par defaut 1 reproduit
    exactement l'ancien comportement, dont depend `hoppe_fast` : le chemin
    wedge^2 explicite (`h0_wedge2_V_explicit`) suppose rank_C = 1 et leverait
    une assertion au-dela. Passer `rank_c_max=None` autorise rank_C >= 2, ce
    qui n'a de sens que pour les fonctions effectivement generalisees --
    h^0(V) et la decomposition de H^1(V). Mesure ayant motive ce parametre :
    sur le scan « gros Gamma », les 26 candidats portant un groupe d'ordre
    compatible sont TOUS E6 a rank_C = 2, et 24 d'entre eux satisfont toutes
    les autres conditions. Ce n'etait donc pas la cohomologie qui les
    ecartait, mais cette seule ligne.
    """
    from cy_landscape.core.exact_cohomology import koszul_cohomology_ex
    if len(c_charges) < 1:
        return False
    if rank_c_max is not None and len(c_charges) > rank_c_max:
        return False
    m = len(c_charges[0])
    # Charges qui doivent PORTER des sections : on exige qu'elles soient
    # positives et que h^0 y soit certifie.
    ch = [list(x) for x in b_charges] + [list(x) for x in c_charges]
    ch += [[b_charges[i][k] + b_charges[j][k] for k in range(m)]
           for i in range(len(b_charges)) for j in range(i + 1, len(b_charges))]
    for c in c_charges:
        ch += [[c[k] + x[k] for k in range(m)] for x in b_charges]
    if any(any(v < 0 for v in x) for x in ch):
        return False

    # Les degres c_j - b_i sont d'une AUTRE nature : ce sont les cases de la
    # matrice f. Une case de degre negatif n'est pas un defaut du modele --
    # elle signifie simplement H^0(O(c_j - b_i)) = 0, donc une case
    # identiquement nulle. Toute la machinerie le traite deja ainsi :
    # `espace_f_equivariant` calcule `actif = all(x >= 0 ...)`,
    # `h0_V_generique` insere un bloc de zeros, et `decomposition_h1_V`
    # saute les cases absentes de `offsets`.
    #
    # Les exiger positives etait donc une condition PLUS STRICTE que ce que
    # le code consommateur demande, et elle avait un cout precis : les SEPT
    # candidats du catalogue portant un groupe Z4 -- cyclique, d'ordre 4,
    # donc la seule route vers le rang 4 exempte du cocycle du §5.27 --
    # etaient tous ecartes « hors domaine » a cause d'UNE case negative
    # chacun, les 36 autres charges etant certifiees.
    for c in c_charges:
        for x in b_charges:
            d = [c[k] - x[k] for k in range(m)]
            if any(v < 0 for v in d):
                continue          # case nulle : rien a certifier
            ch.append(d)

    for x in ch:
        r = koszul_cohomology_ex(ambient, config, x)
        if not r['certified_by_degree'][0]:
            return False
    return True


_RINGS = {}


def get_ring(ambient, config, seed=0):
    """Anneau de sections memoise par geometrie."""
    key = (tuple(ambient), np.asarray(config).tobytes(),
           np.asarray(config).shape, seed)
    if key not in _RINGS:
        _RINGS[key] = Ring(ambient, config, seed=seed)
    return _RINGS[key]
