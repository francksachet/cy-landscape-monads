"""
symetrie_config.py -- Symetries de la matrice de configuration d'une CICY,
et forme canonique d'une monade sous ces symetries.

--------------------------------------------------------------------------
A quoi cela sert
--------------------------------------------------------------------------
Le generateur enumere du §5.23 produit la famille COMPLETE des sommes de
vecteurs unite. Sur une CICY dont la matrice de configuration est
symetrique, cette famille contient les images les unes des autres : #6947
a 12 monades survivantes, qui sont les 12 arrangements d'un meme motif de
multiplicites (3, 1, 1, 0) sur ses quatre facteurs P^1 -- et son groupe
d'automorphismes est justement d'ordre 24, de sorte que ces 12 forment UNE
orbite. Les evaluer separement, c'est refaire douze fois le meme calcul.

`equivariance_f.py` est l'etape qui coute : 14 945 candidats, environ
cinquante heures. Y replier les orbites divise ce cout par |Aut| local.

--------------------------------------------------------------------------
Ce que ce module NE prouve PAS
--------------------------------------------------------------------------
Une permutation des facteurs qui preserve la matrice de configuration
envoie la CICY sur une CICY de la MEME famille. Elle n'a aucune raison, en
revanche, de commuter avec l'action de Gamma lue chez Braun, qui est
attachee a des coordonnees precises. Rien ici ne demontre donc que deux
monades d'une meme orbite recoivent le meme verdict d'equivariance.

C'est pourquoi le repli n'est PAS applique par defaut, et pourquoi
`verifier_invariance` existe : elle confronte la prediction aux verdicts
deja calcules. Tant que cette confrontation n'a pas ete faite sur un lot
reel, la deduplication est une hypothese, pas un resultat -- et une
hypothese qui, si elle est fausse, fait disparaitre des candidats sans
laisser de trace. C'est exactement le mecanisme du §5.23.

|Aut| mesure : #6947 -> 24, #5 -> 12, #6715 -> 6, #6890 -> 1. Sur #6890 le
repli ne peut rien economiser : sa matrice n'a aucune symetrie.
"""
from itertools import permutations

import numpy as np


def automorphismes(ambient, config, plafond_m=9):
    """
    Permutations des facteurs projectifs qui laissent la configuration fixe.

    Une permutation p est retenue si :
      - elle preserve les dimensions des facteurs (P^n -> P^n) ;
      - permuter les COLONNES de la matrice par p redonne le meme
        MULTIENSEMBLE de lignes -- les equations peuvent etre reordonnees,
        elles ne sont pas numerotees.

    Renvoie la liste des permutations (tuples), identite comprise.

    `plafond_m` borne l'enumeration : m! permutations deviennent couteuses
    au-dela de 9 facteurs. Les 194 CICYs a quotient libre ont m <= 10 ; le
    cas m = 10 est traite par le decoupage en blocs de dimensions egales
    ci-dessous, qui ramene le nombre de candidats a prod(n_k!).
    """
    amb = list(ambient)
    cfg = np.asarray(config)
    m = len(amb)
    if cfg.shape[1] != m:
        raise ValueError(f"configuration {cfg.shape} incompatible avec m={m}")

    # Seuls les facteurs de MEME dimension peuvent s'echanger : on enumere
    # les permutations bloc par bloc plutot que les m! globales.
    blocs = {}
    for i, n in enumerate(amb):
        blocs.setdefault(n, []).append(i)
    cout = 1
    for idx in blocs.values():
        for k in range(2, len(idx) + 1):
            cout *= k
    if cout > 10 ** 6:
        # Ne pas engendrer silencieusement une liste ingerable : renvoyer
        # le groupe trivial revient a ne rien deduire, ce qui est correct
        # mais doit se voir.
        return [tuple(range(m))], False

    ref = sorted(map(tuple, cfg))
    autos = []
    listes = [sorted(v) for v in blocs.values()]
    for combo in _produit_permutations(listes):
        p = [0] * m
        for orig, cible in combo:
            p[orig] = cible
        pt = tuple(p)
        if sorted(map(tuple, cfg[:, list(pt)])) == ref:
            autos.append(pt)
    return autos, True


def _produit_permutations(listes):
    """Produit cartesien des permutations de chaque bloc d'indices."""
    if not listes:
        yield []
        return
    tete, reste = listes[0], listes[1:]
    for perm in permutations(tete):
        base = list(zip(tete, perm))
        for suite in _produit_permutations(reste):
            yield base + suite


def canonique(b_charges, c_charges, autos):
    """
    Representant canonique de (B, C) sous le groupe `autos`.

    Les charges sont des multiensembles : deux monades qui ne different que
    par l'ordre des b_i sont la meme. On trie donc avant de comparer, et on
    prend le minimum lexicographique sur toutes les permutations.
    """
    meilleur = None
    for p in autos:
        b = tuple(sorted(tuple(x[i] for i in p) for x in b_charges))
        c = tuple(sorted(tuple(x[i] for i in p) for x in c_charges))
        cle = (b, c)
        if meilleur is None or cle < meilleur:
            meilleur = cle
    return meilleur


def orbites(candidats, ambient, config):
    """
    Regroupe des candidats d'une meme CICY par orbite sous Aut(config).

    `candidats` : liste de dicts portant 'b_charges' et 'c_charges'.
    Renvoie (classes, autos, complet) ou `classes` associe a chaque forme
    canonique la liste des INDICES des candidats de l'orbite.
    """
    autos, complet = automorphismes(ambient, config)
    classes = {}
    for i, r in enumerate(candidats):
        cle = canonique(r['b_charges'], r['c_charges'], autos)
        classes.setdefault(cle, []).append(i)
    return classes, autos, complet


def verifier_invariance(lignes, ambient, config, cles_verdict=None):
    """
    Confronte l'hypothese de repli aux verdicts DEJA calcules.

    Pour chaque orbite, compare le verdict de tous ses membres. Renvoie
    un rapport : nombre d'orbites, nombre d'orbites non triviales, et
    surtout le nombre d'orbites DISCORDANTES -- celles dont les membres
    n'ont pas recu le meme verdict.

    Une seule discordance suffit a invalider le repli : cela signifierait
    qu'evaluer un representant et recopier son verdict sur ses freres
    inventerait des resultats.

    `lignes` : sorties de `equivariance_f.py` (un dict par (candidat,
    groupe, lambda)). Le verdict est compare a (groupe, lambda, survit,
    etat, n_gen_quotient) -- pas au degre temoin, qui est lui-meme permute
    par l'automorphisme et n'a donc aucune raison de coincider.
    """
    if cles_verdict is None:
        cles_verdict = ('groupe', 'lambda', 'survit', 'etat',
                        'n_gen_quotient')
    autos, complet = automorphismes(ambient, config)

    par_orbite = {}
    for L in lignes:
        cle = canonique(L['b_charges'], L['c_charges'], autos)
        v = tuple(str(L.get(k)) for k in cles_verdict)
        par_orbite.setdefault(cle, {}).setdefault(
            tuple(sorted(tuple(x) for x in L['b_charges'])), set()).add(v)

    n_orb = len(par_orbite)
    n_multi = sum(1 for v in par_orbite.values() if len(v) > 1)
    discordantes = []
    for cle, membres in par_orbite.items():
        if len(membres) < 2:
            continue
        jeux = list(membres.values())
        if any(j != jeux[0] for j in jeux[1:]):
            discordantes.append(cle)

    return {
        'ordre_aut': len(autos),
        'groupe_complet': complet,
        'orbites': n_orb,
        'orbites_non_triviales': n_multi,
        'discordantes': len(discordantes),
        'exemples_discordants': discordantes[:3],
        'lignes': len(lignes),
    }
