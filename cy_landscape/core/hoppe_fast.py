"""
hoppe_fast.py -- Critere de Hoppe avec sortie anticipee et cache.

--------------------------------------------------------------------------
PORTEE DU VERDICT -- a lire avant d'interpreter un `stable: True`
--------------------------------------------------------------------------
Le critere s'enonce : pour c1(V) = 0, V est stable si et seulement si
h0(wedge^p V) = 0 pour p = 1..rk-1. Cette EQUIVALENCE suppose Pic(X) de
rang 1. Sur une CICY a m > 1 facteurs projectifs, la stabilite de pente
depend de la classe de Kahler J choisie, et la condition ci-dessus reste
NECESSAIRE sans etre suffisante : elle ne voit pas la polarisation.

`stable: True` doit donc se lire « non elimine par Hoppe », et non
« stable ». C'est un filtre, pas un certificat.

Pour un fibre d'extension, ou la construction exhibe des sous-faisceaux
explicites, `extensions.pente_extension` teste ce que Hoppe ne voit pas et
elimine sur certificat -- 105 des 2 647 extensions declarees Hoppe-stables
d'un scan de controle. Pour une monade, aucun sous-faisceau ne se presente
a bon compte : la reserve reste entiere et non instruite.

Optimisations par rapport a stability_full.py :
1. Teste le rang 1 d'abord (le moins cher) -> sortie immediate si echec
2. Teste un seul H ample (le plus simple) avant les combinaisons
3. Utilise le cache de cohomologie (evite ~90% des recalculs)
4. Skip le rang 3 si rk(V) <= 3
"""
import numpy as np
from itertools import product as iprod

from cy_landscape.core.cache import cached_koszul
from cy_landscape.core.monads import (
    MonadBundle, compute_monad_cohomology, compute_monad_cohomology_ex)
from cy_landscape.core.monad_wedge import cohomology_wedge2_V


def _shift(charges_list, H):
    return [[ch[k] - H[k] for k in range(len(H))] for ch in charges_list]


def _monad_h0_twisted(ambient, config, monad, H):
    """H⁰(V(-H)) via monade decalee, avec cache."""
    shifted = MonadBundle(_shift(monad.b_charges, H), _shift(monad.c_charges, H))

    hB0 = 0
    for b in shifted.b_charges:
        h = cached_koszul(ambient, config, b)
        hB0 += h.get(0, 0)

    hC0 = 0
    for c in shifted.c_charges:
        h = cached_koszul(ambient, config, c)
        hC0 += h.get(0, 0)

    return max(0, hB0 - min(hB0, hC0))


def _wedge2_h0_twisted(ambient, config, monad, H):
    """H⁰(∧²V(-H)) via suite exacte, avec cache."""
    from itertools import combinations
    m = monad.m
    b = monad.b_charges
    c = monad.c_charges[0] if monad.rank_C == 1 else None
    if c is None:
        return 0

    # H⁰(∧²B(-H))
    h0_w2B = 0
    for i, j in combinations(range(monad.rank_B), 2):
        ch = [b[i][k] + b[j][k] - H[k] for k in range(m)]
        h = cached_koszul(ambient, config, ch)
        h0_w2B += h.get(0, 0)

    # H⁰(V⊗C(-H))
    bc_charges = [[b_i[k] + c[k] - H[k] for k in range(m)] for b_i in b]
    h0_BC = sum(cached_koszul(ambient, config, ch).get(0, 0) for ch in bc_charges)
    c2 = [2*c[k] - H[k] for k in range(m)]
    h0_CC = cached_koszul(ambient, config, c2).get(0, 0)
    rank_f = min(h0_BC, h0_CC)
    h0_VC = h0_BC - rank_f

    rank_g = min(h0_w2B, h0_VC)
    return h0_w2B - rank_g


def _wedge3_h0_twisted(ambient, config, monad, H):
    """H⁰(∧³V(-H)) = H³(V(H)) pour rk4 c1=0."""
    neg_H = [-H[k] for k in range(len(H))]
    shifted = MonadBundle(_shift(monad.b_charges, neg_H), _shift(monad.c_charges, neg_H))

    hB = {i: 0 for i in range(4)}
    for b in shifted.b_charges:
        h = cached_koszul(ambient, config, b)
        for i in range(4):
            hB[i] += h.get(i, 0)

    hC = {i: 0 for i in range(4)}
    for c in shifted.c_charges:
        h = cached_koszul(ambient, config, c)
        for i in range(4):
            hC[i] += h.get(i, 0)

    rank_f3 = min(hB[3], hC[3])
    return max(0, hB[3] - rank_f3)


def hoppe_fast(ambient, config, monad, max_H=1, D=None):
    """
    Critere de Hoppe avec sortie anticipee.

    Ordre des tests (du moins cher au plus cher) :
    0. H = 0  <-- L'ENONCE MEME DU CRITERE pour c1(V) = 0
    0 bis. TWISTS : tout H du polytope deg_J(H) >= 0, si `D` est fourni.
       C'est la phase qui voit les H a composantes de signes MELANGES,
       hors de portee de `max_H` : sur le catalogue scan_wilson2 elle
       elimine #7484 (rang 4 SO(10)) par H = (-2, 0, 1), la ou H = 0 et
       H = e_i ne trouvent rien. `D = hoppe_fast.vecteur_D(d_ijk, J)`.
    1. Rang 1, H = e_i (generateurs)
    2. Rang 2, H = e_i
    3. Rang 3, H = e_i (si rk >= 4)
    4. Combinaisons H = (a_1,...,a_m) si max_H > 1

    Sortie immediate au premier echec.

    ----------------------------------------------------------------------
    PHASE 0 -- le test qui manquait
    ----------------------------------------------------------------------
    Le critere de Hoppe, pour c1(V) = 0, s'enonce :

        V stable  <=>  h0(wedge^p V) = 0  pour p = 1..rk-1

    soit exactement H = 0. Les phases 1 et 2 ne testaient que H = e_i et
    ses combinaisons. Or pour H ample, H0(wedge^p V(-H)) est un SOUS-espace
    de H0(wedge^p V) : ces tests sont donc strictement plus FAIBLES que le
    critere et laissent passer des fibres non stables. C'est ce qui
    expliquait les candidats a h0(V) != 0 ou h3(V) != 0 des scans
    precedents -- impossibles pour un fibre stable de pente nulle.

    En depliant det V = O et la dualite de Serre (h0(F*) = h3(F) sur un
    CY3), wedge^(rk-1) V = V* et wedge^(rk-2) V = wedge^2 V*, d'ou :

        rk = 3 : h0(V), h3(V)
        rk = 4 : h0(V), h0(w2V), h3(V)
        rk = 5 : h0(V), h0(w2V), h3(w2V), h3(V)

    En rang 3 le critere se reduit donc entierement a h0(V) = h3(V) = 0.

    ATTENTION : `_wedge3_h0_twisted` implemente h0(w3V) = h3(V), identite
    valable en rang 4 SEULEMENT (voir sa docstring), alors que les phases
    1 et 2 l'appellent des que rk >= 4 -- donc aussi en rang 5, ou elle est
    fausse. La phase 0 n'y fait pas appel et utilise le tableau ci-dessus.
    """
    m = monad.m
    rk = monad.rank_V
    tests = 0

    # ---- Phase 0 : H = 0 -------------------------------------------------
    H0 = [0] * m
    tests += 1
    h0 = _monad_h0_twisted(ambient, config, monad, H0)
    if h0 > 0:
        return {"stable": False, "reason": f"H⁰(V) = {h0}", "tests": tests}

    # ---- Phase 0 bis : le polytope des twists ---------------------------
    # Purement ELIMINATOIRE : un h0(V(-H)) > 0 certifie avec deg_J(H) >= 0
    # exhibe un sous-faisceau de pente >= mu(V) = 0. L'absence de twist
    # destabilisant ne prouve rien -- le verdict `stable` de cette fonction
    # reste « non elimine » (voir l'entete du module).
    if D is not None:
        tw = hoppe_twists(ambient, config, monad, D)
        tests += tw['n_twists']
        if tw['instable']:
            return {"stable": False, "reason": tw['motif'], "tests": tests,
                    "twist_temoin": tw['temoin']}

    # h3(V) : on utilise les BORNES rigoureuses, sans supposer la stabilite
    # (ce serait circulaire ici). Une borne INFERIEURE strictement positive
    # prouve h3(V) > 0, donc la non-stabilite. Une borne inferieure nulle ne
    # permet pas de conclure et n'est pas traitee comme un succes.
    cohV = compute_monad_cohomology_ex(monad, ambient, config,
                                       stable_hypothesis=False)
    if cohV is not None:
        tests += 1
        lo3 = cohV['bounds'][3][0]
        if lo3 > 0:
            return {"stable": False, "reason": f"H³(V) >= {lo3}", "tests": tests}

    if rk >= 4:
        if monad.rank_C != 1:
            # wedge^2 V n'est pas calculable par ce chemin : on ne peut pas
            # conclure. Signale explicitement plutot que de faire passer.
            return {"stable": False,
                    "reason": "non testable (rank_C >= 2, wedge^2 V indisponible)",
                    "tests": tests}
        # RESERVE IMPORTANTE : `cohomology_wedge2_V` est defaillant --
        # sur 456 cas testes, 79 % des vecteurs de Betti retournes
        # contredisent le chi calcule par la fonction elle-meme. Les
        # verdicts de rang 4 et 5 issus de ce chemin sont donc AUSSI
        # fiables que ce module, c'est-a-dire peu. Seuls les verdicts de
        # rang 3 (qui ne font intervenir que h0(V) et h3(V)) sont
        # actuellement dignes de confiance.
        try:
            res = cohomology_wedge2_V(ambient, config,
                                      monad.b_charges, monad.c_charges,
                                      cohom_V=(dict(cohV) if cohV else None))
        except Exception:
            return {"stable": False,
                    "reason": "non testable (echec du calcul de wedge^2 V)",
                    "tests": tests}
        w2 = dict(res.get('wedge2V') or {})
        cert = dict(res.get('certified') or {})

        # --- RANGS REELS : on tente le calcul explicite de h0(wedge^2 V) ---
        # Les bornes de la suite exacte longue ne determinent presque jamais
        # h0(wedge^2 V), ce qui faisait sortir les rangs 4 et 5 en « non
        # testable ». `sections.py` construit la matrice de
        # H0(wedge^2 B) -> H0(B (x) C) et en calcule le rang : le noyau donne
        # h0(wedge^2 V) directement. Voir l'en-tete de sections.py pour le
        # domaine de validite et la validation (174/177 sur h0(V), 38/38
        # dans les bornes sur h0(wedge^2 V)).
        try:
            from cy_landscape.core.sections import (
                domaine_valide, get_ring, h0_wedge2_V_explicit)
            if domaine_valide(ambient, config, monad.b_charges, monad.c_charges):
                anneau = get_ring(ambient, config)
                val, _dims = h0_wedge2_V_explicit(
                    anneau, monad.b_charges, monad.c_charges, maxdim=1500)
                if val is not None:
                    w2[0] = int(val)
                    cert[0] = True
        except Exception:
            pass

        # Un degre absent veut dire INCONNU, pas zero. On ne conclut pas a
        # la stabilite sur une valeur non determinee : c'est exactement le
        # mecanisme par lequel la v1 laissait passer des fibres non stables.
        degres = [0] + ([3] if rk == 5 else [])
        for dgr in degres:
            tests += 1
            if not cert.get(dgr):
                return {"stable": False,
                        "reason": f"non testable (h^{dgr}(∧²V) non determine, "
                                  f"bornes {res.get('bounds', {}).get(dgr)})",
                        "tests": tests}
            if int(w2.get(dgr, 0)) > 0:
                return {"stable": False,
                        "reason": f"H^{dgr}(∧²V) = {int(w2.get(dgr, 0))}",
                        "tests": tests}

    # Phase 1 : generateurs e_i, rang par rang
    for p in range(1, rk):
        for i in range(m):
            H = [0] * m
            H[i] = 1
            tests += 1

            if p == 1:
                h0 = _monad_h0_twisted(ambient, config, monad, H)
            elif p == 2:
                h0 = _wedge2_h0_twisted(ambient, config, monad, H)
            elif p == 3 and rk >= 4:
                h0 = _wedge3_h0_twisted(ambient, config, monad, H)
            else:
                continue

            if h0 > 0:
                return {"stable": False,
                        "reason": f"H⁰(∧^{p}V(-e_{i+1})) = {h0}",
                        "tests": tests}

    # Phase 2 : combinaisons (si max_H > 1)
    if max_H > 1:
        for H in iprod(*[range(1, max_H + 1) for _ in range(m)]):
            H = list(H)
            if all(h <= 1 for h in H):
                continue  # Deja teste

            for p in range(1, rk):
                tests += 1
                if p == 1:
                    h0 = _monad_h0_twisted(ambient, config, monad, H)
                elif p == 2:
                    h0 = _wedge2_h0_twisted(ambient, config, monad, H)
                elif p == 3 and rk >= 4:
                    h0 = _wedge3_h0_twisted(ambient, config, monad, H)
                else:
                    continue

                if h0 > 0:
                    return {"stable": False,
                            "reason": f"H⁰(∧^{p}V(-{H})) = {h0}",
                            "tests": tests}

    return {"stable": True,
            "reason": f"Hoppe satisfait ({tests} tests, rangs 1..{rk-1})",
            "tests": tests}


# ======================================================================
# LE POLYTOPE DES TWISTS -- la forme SUFFISANTE du critere (§5.14)
# ======================================================================


def vecteur_D(d_ijk, J):
    """D_k(J) = sum_jk d_kjl J_j J_l. deg_J(O(a)) = sum_k a_k D_k(J)."""
    d = np.asarray(d_ijk, dtype=np.int64)
    Jv = np.asarray(J, dtype=np.int64)
    return np.einsum('ijk,j,k->i', d, Jv, Jv)


def charges_wedge_B(b_charges, p_ext):
    """Charges de wedge^p B : sommes sur les p-uplets croissants de b_i."""
    from itertools import combinations
    m = len(b_charges[0])
    return [[sum(b_charges[i][k] for i in idx) for k in range(m)]
            for idx in combinations(range(len(b_charges)), p_ext)]


def polytope_twists(b_charges, p_ext, D, plafond=200000):
    """
    TOUS les H a tester pour rendre le critere de Hoppe SUFFISANT au degre p.

    ----------------------------------------------------------------------
    Pourquoi cet ensemble, et pourquoi il est fini
    ----------------------------------------------------------------------
    V est mu_J-stable DES QUE h^0(wedge^p V(-H)) = 0 pour p = 1..rk-1 et
    tout O(H) de degre deg_J(H) >= 0. Ecrit ainsi l'ensemble des H parait
    infini. Il ne l'est pas :

      - wedge^p V(-H) est un sous-faisceau de wedge^p B(-H), donc
        h^0(wedge^p V(-H)) <= sum_I h^0(O(ch_I - H)). Dans le modele
        S/I employe ici, h^0(O(a)) = 0 des qu'une composante de a est
        negative. Il faut donc H <= ch_I pour AU MOINS un I : H est borne
        SUPERIEUREMENT par les charges de wedge^p B.

      - deg_J(H) = sum_k H_k D_k(J) >= 0 avec tous les D_k(J) > 0 et
        H_k <= hi_k borne alors H INFERIEUREMENT :
        H_k >= -(sum_{l != k} hi_l D_l) / D_k.

    Le polytope est donc compact, et petit en pratique : 110 twists pour
    #6890 et 143 pour #6947, tous degres p confondus.

    Renvoie None si un D_k(J) est nul -- la borne inferieure disparait et
    l'ensemble n'est plus fini par cet argument -- ou si le decompte depasse
    `plafond`. Dans les deux cas l'appelant ne doit PAS conclure a la
    stabilite : c'est un critere suffisant, il n'a de valeur que complet.
    """
    from itertools import product
    ch = charges_wedge_B(b_charges, p_ext)
    if not ch:
        return []
    m = len(ch[0])
    D = [int(x) for x in D]
    if any(x <= 0 for x in D):
        return None
    hi = [max(c[k] for c in ch) for k in range(m)]
    lo = []
    for k in range(m):
        reste = sum(hi[l] * D[l] for l in range(m) if l != k)
        lo.append(-(reste // D[k]) - 1)
    taille = 1
    for k in range(m):
        taille *= (hi[k] - lo[k] + 1)
        if taille > plafond:
            return None
    out = []
    for H in product(*[range(lo[k], hi[k] + 1) for k in range(m)]):
        if sum(H[k] * D[k] for k in range(m)) < 0:
            continue
        if any(all(H[k] <= c[k] for k in range(m)) for c in ch):
            out.append(list(H))
    return out


def borne_h0_V_twist(ambient, config, monad, H):
    """
    (borne_inferieure, certifie) pour h0(V(-H)), V = ker(B -> C).

    dim ker >= dim source - dim cible : la borne est INCONDITIONNELLE, et
    une valeur strictement positive PROUVE h0(V(-H)) > 0.

    Chaque h0 de fibre en droites doit etre CERTIFIE par
    `koszul_cohomology_ex` : sans cela la borne porterait sur des nombres
    faux dans ~30 % des cas (§4.2), et une elimination « demontree » ne le
    serait pas. `certifie = False` => l'appelant ne doit pas conclure.
    """
    from cy_landscape.core.exact_cohomology import koszul_cohomology_ex
    m = len(H)
    src = dst = 0
    cert = True
    for b in monad.b_charges:
        r = koszul_cohomology_ex(ambient, config, [b[k] - H[k] for k in range(m)])
        cert = cert and bool(r['certified_by_degree'][0])
        src += r[0]
    for c in monad.c_charges:
        r = koszul_cohomology_ex(ambient, config, [c[k] - H[k] for k in range(m)])
        cert = cert and bool(r['certified_by_degree'][0])
        dst += r[0]
    return max(0, src - dst), cert


def hoppe_twists(ambient, config, monad, D, plafond=200000):
    """
    Phase des TWISTS : cherche un H de degre >= 0 destabilisant V, par la
    seule borne de comptage. ELIMINATOIRE, jamais confirmatoire.

    Un h0(V(-H)) > 0 avec deg_J(H) >= 0 fournit un sous-faisceau O(H) c V
    de pente >= 0 = mu(V) : V n'est pas mu_J-stable. C'est ce que les
    phases H = 0 et H = e_i de `hoppe_fast` ne voient pas -- le polytope
    complet contient des H a composantes de signes MELANGES, hors de
    portee de `max_H`.

    Mesure sur le catalogue `scan_wilson2` : sur 77 entrees a rank_C = 1,
    une (#7484, rang 4 SO(10)) est eliminee par H = (-2, 0, 1), de degre 4,
    avec source 13 > cible 12 et tous les h0 certifies. Elle figurait comme
    Hoppe-stable. Aucune n'est eliminee par H = 0 ni par H = e_i.

    Renvoie {'instable': bool|None, 'temoin': H|None, 'n_twists': int,
             'non_certifies': int}. `instable = None` signifie qu'on n'a
    rien trouve mais que des twists n'ont pas pu etre certifies : on ne
    conclut pas a la stabilite pour autant -- ce test ne la prouve jamais.
    """
    Hs = polytope_twists(monad.b_charges, 1, D, plafond=plafond)
    if Hs is None:
        return {'instable': None, 'temoin': None, 'n_twists': 0,
                'non_certifies': 0, 'motif': 'polytope non borne ou hors plafond'}
    n_nc = 0
    for H in Hs:
        borne, cert = borne_h0_V_twist(ambient, config, monad, H)
        if not cert:
            n_nc += 1
            continue
        if borne > 0:
            return {'instable': True, 'temoin': list(H), 'n_twists': len(Hs),
                    'non_certifies': n_nc,
                    'motif': f'h0(V(-H)) >= {borne} avec H = {list(H)}, '
                             f'deg_J(H) >= 0'}
    return {'instable': (None if n_nc else False), 'temoin': None,
            'n_twists': len(Hs), 'non_certifies': n_nc,
            'motif': (f'{n_nc} twists non certifies sur {len(Hs)}' if n_nc
                      else f'{len(Hs)} twists testes, aucun destabilisant')}
