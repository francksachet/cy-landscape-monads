"""
cech.py -- Base explicite des classes que le modele S_a / I_a ne voit pas.

--------------------------------------------------------------------------
D'ou vient ce module
--------------------------------------------------------------------------
Le §5.29 a montre que S_a / I_a SOUS-COMPTE H^0(Y, O(a)) : sur #6836,
la charge (0,0,0,0,1) donne 4 la ou h^0 vaut 8. Le §5.30 a localise le
manque : H^0(Y) recoit toute la diagonale q = p de la suite spectrale

    E_1^{-p, q} = H^q(A, /\\^p)  ==>  H^{q-p}(Y, O(a))

et le modele ne represente que le coin (0, 0). `analyse_modele` COMPTE ce
qui manque ; ce module en construit une BASE.

--------------------------------------------------------------------------
La forme des classes
--------------------------------------------------------------------------
Sur A = produit de P^{n_i}, la formule de Kunneth donne

    H^q(A, O(b)) = (+)_{q_1 + ... + q_m = q} (x)_i H^{q_i}(P^{n_i}, O(b_i))

et chaque facteur ne contribue qu'en deux endroits (Bott) :

    q_i = 0    si b_i >= 0        base : les monomes de degre b_i
    q_i = n_i  si b_i <= -n_i - 1 base : les monomes NEGATIFS, exposants
                                  tous <= -1 de somme b_i
    rien sinon

Un monome negatif x_0^{a_0} ... x_n^{a_n} avec tous les a_j <= -1 est la
notation usuelle de la classe de Cech correspondante sur le recouvrement
standard : pour n = 1 et b = -2, l'unique element est 1/(x_0 x_1), qui
engendre H^1(P^1, O(-2)).

Une classe manquante de H^0(Y, O(a)) est donc un couple (S, w) : un
sous-ensemble S d'equations, de cardinal p, et un element w de la base de
H^p(A, O(a - somme_{k dans S} d_k)).

--------------------------------------------------------------------------
Ce que ce module fait, et ne fait pas
--------------------------------------------------------------------------
Il ENUMERE. Il ne multiplie pas encore, et ne fait pas agir Gamma : ce sont
les deux etapes suivantes du §5.30, et ce sont elles qui rendront les
candidats Z4 calculables. Le controle disponible des maintenant est que le
cardinal de la base redonne exactement le `manquant` d'`analyse_modele`,
lui-meme valide contre Koszul sans contre-exemple.
"""
from itertools import combinations
from math import comb


def monomes_positifs(n, d):
    """
    Exposants des monomes de degre d en n+1 variables : base de
    H^0(P^n, O(d)) pour d >= 0. Liste vide sinon.
    """
    if d < 0:
        return []
    return list(_compositions(n + 1, d))


def monomes_negatifs(n, d):
    """
    Base de H^{n}(P^n, O(d)) pour d <= -n-1 : les exposants tous <= -1 de
    somme d. Liste vide sinon.

    Ecriture : a_j = -1 - c_j avec c_j >= 0 et somme(c_j) = -d - (n+1),
    ce qui donne bien C(-d-1, n) elements -- le cardinal de Bott.
    """
    if d > -n - 1:
        return []
    reste = -d - (n + 1)
    return [tuple(-1 - c for c in comp)
            for comp in _compositions(n + 1, reste)]


def _compositions(k, total):
    """Suites de k entiers >= 0 de somme `total`."""
    if k == 1:
        yield (total,)
        return
    for premier in range(total + 1):
        for suite in _compositions(k - 1, total - premier):
            yield (premier,) + suite


def base_hq_facteur(n, d, q):
    """Base de H^q(P^n, O(d)) -- vide hors des deux cas de Bott."""
    if q == 0:
        return monomes_positifs(n, d)
    if q == n:
        return monomes_negatifs(n, d)
    return []


def base_hq_ambiant(ambient, b, q):
    """
    Base de H^q(A, O(b)) par Kunneth.

    Un element est un tuple d'exposants, un par facteur projectif. Les
    facteurs ou q_i = n_i portent des exposants NEGATIFS : c'est ce qui
    distingue une classe de Cech d'une section ordinaire.
    """
    m = len(ambient)
    sorties = []
    # repartitions admissibles de q : chaque facteur contribue 0 ou n_i
    for choix in _repartitions(ambient, list(b), q):
        bases = [base_hq_facteur(ambient[i], b[i], choix[i])
                 for i in range(m)]
        if any(not x for x in bases):
            continue
        for combo in _produit(bases):
            sorties.append(tuple(combo))
    return sorties


def _repartitions(ambient, b, q, i=0, courant=None):
    """Suites (q_i) avec q_i dans {0, n_i}, somme q, et H^{q_i} non nul."""
    if courant is None:
        courant = []
    m = len(ambient)
    if i == m:
        if q == 0:
            yield list(courant)
        return
    n = ambient[i]
    if b[i] >= 0:                       # q_i = 0 possible
        yield from _repartitions(ambient, b, q, i + 1, courant + [0])
    if b[i] <= -n - 1 and q - n >= 0:   # q_i = n_i possible
        yield from _repartitions(ambient, b, q - n, i + 1, courant + [n])


def _produit(listes):
    if not listes:
        yield []
        return
    for tete in listes[0]:
        for suite in _produit(listes[1:]):
            yield [tete] + suite


def base_classes_manquantes(ambient, config, a, p_max=None):
    """
    Base des classes de H^0(Y, O(a)) que S_a / I_a ne represente pas.

    Rend une liste de (S, w) : S le sous-ensemble d'equations (tuple
    d'indices, de cardinal p >= 1), w un element de la base de
    H^p(A, O(a - somme_S d_k)) sous la forme d'un tuple d'exposants par
    facteur.

    Le cardinal doit valoir `analyse_modele(...)['manquant']`.
    """
    # int() explicite : `config` vient de numpy, et des np.int64 se
    # propageraient dans les exposants -- illisibles a l'affichage et
    # sournois a la comparaison (np.int64(0) == 0 est vrai, mais les
    # tuples ne se comparent pas toujours comme on croit).
    cfg = [[int(v) for v in r] for r in config]
    K = len(cfg)
    m = len(ambient)
    a = [int(v) for v in a]
    if p_max is None:
        p_max = K
    out = []
    for p in range(1, min(p_max, K) + 1):
        for S in combinations(range(K), p):
            b = [a[i] - sum(cfg[k][i] for k in S) for i in range(m)]
            for w in base_hq_ambiant(ambient, b, p):
                out.append((S, w))
    return out


def cardinal_hq(ambient, b, q):
    """Cardinal de `base_hq_ambiant`, sans rien construire (controle)."""
    total = 0
    for choix in _repartitions(ambient, list(b), q):
        d = 1
        for i, n in enumerate(ambient):
            if choix[i] == 0:
                d *= comb(b[i] + n, n) if b[i] >= 0 else 0
            else:
                d *= comb(-b[i] - 1, n) if b[i] <= -n - 1 else 0
        total += d
    return total


# ======================================================================
# Produit par une section ordinaire
# ======================================================================
#
# LA REGLE, ET POURQUOI ELLE EST SI SIMPLE
# ----------------------------------------
# Sur le recouvrement standard de P^n, H^n(P^n, O(d)) est le quotient des
# monomes de Laurent de degre d par ceux qui ont AU MOINS un exposant >= 0
# (ceux-la viennent des faces C^{n-1}). Une base est donc donnee par les
# monomes a exposants tous <= -1, et multiplier par un monome ordinaire
# revient a additionner les exposants PUIS a projeter :
#
#     le produit survit  <=>  tous les exposants restent <= -1
#
# Verifie a la main sur P^1 : 1/(x0 x1) fois x0 vaut zero (l'exposant de x0
# passe a 0, on retombe dans l'image de C^0), tandis que 1/(x0^2 x1) fois x0
# donne 1/(x0 x1). C'est tout -- il n'y a pas de terme correctif a ce
# niveau, la difficulte est ailleurs (voir la note en fin de module).

def produit_classe(w, monome):
    """
    Multiplie une classe de Kunneth `w` par un monome ordinaire.

    `w`     : tuple d'exposants par facteur (certains facteurs negatifs).
    `monome`: tuple d'exposants par facteur, tous >= 0.

    Rend le tuple produit, ou None si la classe MEURT -- c'est-a-dire si un
    facteur a exposants negatifs voit l'un d'eux atteindre 0.
    """
    sortie = []
    for exp, mon in zip(w, monome):
        negatif = any(v < 0 for v in exp)
        prod = tuple(int(a) + int(b) for a, b in zip(exp, mon))
        if negatif and any(v >= 0 for v in prod):
            return None
        sortie.append(prod)
    return tuple(sortie)


def matrice_produit(base_src, base_dst, polynome):
    """
    Matrice du produit par `polynome` de <base_src> vers <base_dst>.

    `polynome` : liste de (coefficient, monome). Les bases sont des listes
    de classes au format de `base_hq_ambiant`.

    Les classes qui meurent ne contribuent a rien -- c'est la projection du
    quotient de Cech, pas une troncature arbitraire.
    """
    index = {w: j for j, w in enumerate(base_dst)}
    M = [[0] * len(base_src) for _ in range(len(base_dst))]
    for j, w in enumerate(base_src):
        for coef, mon in polynome:
            r = produit_classe(w, mon)
            if r is None:
                continue
            i = index.get(r)
            if i is None:
                # Le produit est une classe valide du degre cible mais
                # absente de la base fournie : c'est une erreur d'appel,
                # pas un cas a ignorer en silence.
                raise KeyError(f"classe produite hors de la base cible : {r}")
            M[i][j] += coef
    return M


# ----------------------------------------------------------------------
# CE QUI N'EST PAS ENCORE FAIT, ET OU EST LA DIFFICULTE
# ----------------------------------------------------------------------
# Le produit ci-dessus est celui de la suite spectrale ASSOCIEE GRADUEE.
# Dans H^0(Y, O(a)), les sections venant de l'ambiant forment un
# SOUS-espace (le coin E_inf^{0,0}) et les classes de Cech n'en sont qu'un
# QUOTIENT. La matrice du produit dans la base (ordinaire, Cech) est donc
# triangulaire par blocs
#
#     [ ordinaire -> ordinaire      correction ]
#     [ 0                     Cech -> Cech     ]
#
# et le bloc `correction`, de Cech vers ordinaire, n'est pas donne par la
# regle monomiale : il vient de la differentielle de Koszul. C'est lui qui
# reste a construire, et c'est le seul endroit du chantier dont je ne
# puisse pas dire d'avance qu'il tombe juste.
