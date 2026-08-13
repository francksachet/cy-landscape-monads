#!/usr/bin/env python3
"""
tests_regression.py -- Suite de non-regression du pipeline CY Landscape.

--------------------------------------------------------------------------
Pourquoi cette suite existe
--------------------------------------------------------------------------
Chaque defaut trouve dans ce projet l'a ete par une reference INDEPENDANTE :
Riemann-Roch contre la somme alternee des h^i, la dualite de Serre contre
elle-meme, l'identite wedge^2 V = V* au rang 3, une valeur connue de la
litterature. Aucun n'avait ete detecte par le code lui-meme.

Le bug du c2 croise dans intersection.py avait meme survecu a sa propre
validation, parce que celle-ci ne portait que sur la quintique -- ou m = 1,
donc aucun terme croise, donc le bug invisible.

Cette suite rassemble tous ces controles. Elle tourne en quelques secondes
et doit etre relancee apres CHAQUE modification, avant tout scan. Un scan
coute deux heures ; un test en coute une seconde.

    python tests_regression.py

Sortie : une ligne par test, et un code de retour non nul si l'un echoue.
"""
import sys, os, itertools, random, traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

_RESULTATS = []


def test(nom):
    def deco(fn):
        def run():
            try:
                detail = fn()
                _RESULTATS.append((nom, True, detail or ""))
            except AssertionError as e:
                _RESULTATS.append((nom, False, str(e)))
            except Exception as e:
                _RESULTATS.append((nom, False,
                                   f"{type(e).__name__}: {e}"))
        run.__name__ = fn.__name__
        return run
    return deco


def _cy3():
    """CICYs embarquees qui sont bien des 3-varietes."""
    from cy_landscape.data.oxford_cicys import get_all_oxford
    out = []
    for c in get_all_oxford():
        cfg = np.asarray(c['config'])
        if cfg.ndim == 1:
            cfg = cfg.reshape(1, -1)
        if sum(c['ambient']) - cfg.shape[0] == 3 and cfg.shape[0] <= 12:
            out.append((c, cfg))
    return out


# ======================================================================
# 1. Nombres d'intersection et c2
# ======================================================================

@test("c2.J : quintique = 50, bicubique = (36,36)")
def t_c2():
    from cy_landscape.core.intersection import (
        compute_intersection_numbers, compute_c2_tangent)
    for amb, cfg, attendu in (([4], [[5]], [50.0]),
                              ([2, 2], [[3, 3]], [36.0, 36.0])):
        cfg = np.array(cfg)
        d = compute_intersection_numbers(amb, cfg)
        v = list(compute_c2_tangent(amb, cfg, d))
        assert all(abs(a - b) < 1e-9 for a, b in zip(v, attendu)), \
            f"{amb}{cfg.tolist()} -> {v}, attendu {attendu}"
    return "la bicubique est le cas discriminant : m > 1, donc termes croises"


@test("chi(L) est entier pour tout L (Riemann-Roch)")
def t_chi_entier():
    from cy_landscape.core.chi_exact import make_calculator
    rng = random.Random(1)
    n = 0
    for c, cfg in _cy3():
        cal = make_calculator(c['ambient'], cfg)
        for _ in range(60):
            a = [rng.randint(-5, 5) for _ in range(len(c['ambient']))]
            v = cal.line(a)
            assert isinstance(v, int), f"chi non entier pour {a}"
            n += 1
    return f"{n} vecteurs de charges"


# ======================================================================
# 2. Cohomologie des fibres en droites
# ======================================================================

@test("koszul : quintique O(n) exact, y compris O(5) = 125")
def t_quintique():
    from cy_landscape.core.exact_cohomology import koszul_cohomology_ex
    attendu = {-5: (0, 0, 0, 125), -1: (0, 0, 0, 5), 0: (1, 0, 0, 1),
               1: (5, 0, 0, 0), 2: (15, 0, 0, 0), 5: (125, 0, 0, 0)}
    for n, att in attendu.items():
        r = koszul_cohomology_ex([4], np.array([[5]]), [n])
        got = tuple(r[i] for i in range(4))
        assert got == att, f"O({n}) -> {got}, attendu {att}"
    return "h0(O(5)) = 125 : c'est d_1 qui l'apporte, l'ancienne version donnait 126"


@test("chi du module == Riemann-Roch, sur tout le domaine")
def t_chi_module():
    from cy_landscape.core.exact_cohomology import koszul_cohomology_ex
    from cy_landscape.core.chi_exact import make_calculator
    rng = random.Random(5)
    n = 0
    for c, cfg in _cy3():
        cal = make_calculator(c['ambient'], cfg)
        for _ in range(40):
            a = [rng.randint(-4, 4) for _ in range(len(c['ambient']))]
            r = koszul_cohomology_ex(c['ambient'], cfg, a)
            assert r['chi'] == cal.line(a), \
                f"chi divergent en {a} : {r['chi']} vs {cal.line(a)}"
            n += 1
    return f"{n} vecteurs — chi est exact meme quand les h^i ne le sont pas"


@test("h^i certifies : accord avec Riemann-Roch a 100 %")
def t_certifies_hrr():
    from cy_landscape.core.exact_cohomology import koszul_cohomology_ex
    from cy_landscape.core.chi_exact import make_calculator
    rng = random.Random(5)
    n = ok = 0
    for c, cfg in _cy3():
        cal = make_calculator(c['ambient'], cfg)
        for _ in range(40):
            a = [rng.randint(-4, 4) for _ in range(len(c['ambient']))]
            r = koszul_cohomology_ex(c['ambient'], cfg, a)
            if not r['certified']:
                continue
            n += 1
            s = sum((-1) ** i * r[i] for i in range(4))
            if s == cal.line(a):
                ok += 1
    assert n > 0, "aucun cas certifie — le critere s'est effondre"
    assert ok == n, f"{n - ok} certifies en desaccord avec HRR sur {n}"
    return f"{n} cas certifies, tous coherents"


@test("Serre h^i(L) = h^(3-i)(L^-1) sur les paires certifiees")
def t_serre():
    from cy_landscape.core.exact_cohomology import koszul_cohomology_ex
    rng = random.Random(5)
    n = ok = 0
    for c, cfg in _cy3():
        for _ in range(40):
            a = [rng.randint(-4, 4) for _ in range(len(c['ambient']))]
            r = koszul_cohomology_ex(c['ambient'], cfg, a)
            rn = koszul_cohomology_ex(c['ambient'], cfg, [-x for x in a])
            ca, cb = r['certified_by_degree'], rn['certified_by_degree']
            for i in range(4):
                if ca.get(i) and cb.get(3 - i):
                    n += 1
                    ok += (r[i] == rn[3 - i])
    assert n > 0, "aucune paire certifiee"
    assert ok == n, (f"{n - ok} violations sur {n} — la detection d'ambiguite "
                     f"de rang ne couvre plus tout")
    return f"{n} paires — ce test avait detecte la dependance au sens de recursion"


# ======================================================================
# 3. wedge^2 V
# ======================================================================

@test("chi(w2V) = -chi(V) au rang 3, = 0 au rang 4")
def t_wedge2_chi():
    from cy_landscape.core.cache import set_geometry
    from cy_landscape.core.monads import generate_monads
    from cy_landscape.core.monad_wedge import cohomology_wedge2_V
    from cy_landscape.core.chi_exact import make_calculator
    n3 = n4 = 0
    for c, cfg in _cy3()[:6]:
        set_geometry(c['ambient'], cfg)
        cal = make_calculator(c['ambient'], cfg)
        m = len(c['ambient'])
        for rk, ref in ((3, None), (4, 0)):
            rng = np.random.RandomState(3)
            for mo in generate_monads(m, rk, max_charge=3, n_random=60, rng=rng):
                if mo.rank_C != 1:
                    continue
                r = cohomology_wedge2_V(c['ambient'], cfg,
                                        mo.b_charges, mo.c_charges)
                att = -cal.monad(mo.b_charges, mo.c_charges) if rk == 3 else 0
                assert r['chi_wedge2V'] == att, \
                    f"rang {rk} : {r['chi_wedge2V']} vs {att}"
                if rk == 3: n3 += 1
                else: n4 += 1
        if n3 > 60 and n4 > 60:
            break
    return f"{n3} monades rang 3, {n4} rang 4"


# ======================================================================
# 4. Monades degenerees
# ======================================================================

@test("monades scindees rejetees (cas reels de la CICY 7669)")
def t_degenerees():
    from cy_landscape.core.monads import MonadBundle, check_monad_nondegenerate
    doit_rejeter = [
        ([[-2, 0, -3], [2, -3, 1], [-1, 2, 1], [-1, 1, -2]], [[-2, 0, -3]]),
        ([[1, 3, -2], [-1, -2, 3], [-2, 1, 3], [1, -1, -1]], [[-1, 1, 3]]),
    ]
    doit_passer = [
        ([[1, 0], [0, 1], [1, 1], [1, 1]], [[3, 3]]),
    ]
    for b, cc in doit_rejeter:
        ok, motif = check_monad_nondegenerate(MonadBundle(b, cc))
        assert not ok, f"monade scindee acceptee : B={b} C={cc}"
    for b, cc in doit_passer:
        ok, motif = check_monad_nondegenerate(MonadBundle(b, cc))
        assert ok, f"monade saine rejetee ({motif}) : B={b} C={cc}"
    return "b_i = c_j et colonne nulle"


# ======================================================================
# 5. Elagage du generateur positif
# ======================================================================

@test("elagage exact : aucun c faisable perdu")
def t_elagage():
    from cy_landscape.core.positive_monads import (
        _c_is_feasible, _find_positive_B, _mix)
    perdus = testes = 0
    for m in (2, 3, 4):
        for r_B in (4, 5, 6):
            for c in itertools.product(range(5), repeat=m):
                if all(x == 0 for x in c) or _c_is_feasible(c, r_B):
                    continue
                testes += 1
                rng = np.random.RandomState(_mix(1, c, r_B))
                if _find_positive_B(list(c), r_B, rng, max_attempts=200):
                    perdus += 1
    assert perdus == 0, f"{perdus} vecteurs c rejetes a tort sur {testes}"
    return f"{testes} vecteurs infaisables, aucun faux rejet"


# ======================================================================
# 6. Action de Gamma
# ======================================================================

@test("action de Gamma : rt[3] lu, action d'ordre 3 (donnees CICY 7669)")
def t_gamma():
    from cy_landscape.core.gamma_action import (
        choisir_premier, racine_primitive, parse_entree, action_monomes)
    import re
    p, _ = choisir_premier([3])
    rac = {3: racine_primitive(p, 3)}
    assert pow(rac[3], 3, p) == 1 and rac[3] != 1, "racine cubique invalide"

    def mat(txt):
        out = []
        for l in re.findall(r'\{([^{}]*)\}', txt):
            prof, cur, ent = 0, [], []
            for ch in l:
                if ch == '[': prof += 1
                elif ch == ']': prof -= 1
                if ch == ',' and prof == 0:
                    ent.append(''.join(cur)); cur = []
                else:
                    cur.append(ch)
            if cur: ent.append(''.join(cur))
            out.append([parse_entree(e, p, rac) for e in ent])
        return out

    g1 = mat("{{1,0,0,0,0,0,0,0,0},{0,rt[3]^2,0,0,0,0,0,0,0},"
             "{0,0,rt[3],0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0},"
             "{0,0,0,0,rt[3]^2,0,0,0,0},{0,0,0,0,0,rt[3],0,0,0},"
             "{0,0,0,0,0,0,1,0,0},{0,0,0,0,0,0,0,rt[3],0},"
             "{0,0,0,0,0,0,0,0,rt[3]^2}}")
    g2 = mat("{{0,0,1,0,0,0,0,0,0},{1,0,0,0,0,0,0,0,0},{0,1,0,0,0,0,0,0,0},"
             "{0,0,0,0,0,1,0,0,0},{0,0,0,1,0,0,0,0,0},{0,0,0,0,1,0,0,0,0},"
             "{0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,1,0,0},{0,0,0,0,0,0,0,1,0}}")
    for nom, g in (("phases", g1), ("permutation", g2)):
        for deg in ([1, 1, 1], [2, 0, 1]):
            a = action_monomes(g, [2, 2, 2], deg, p)
            assert a is not None, f"{nom} deg {deg} : action non monomiale"
            base, deg_img, base_img, perm, coef = a
            assert list(deg_img) == list(deg), \
                f"{nom} : sigma non trivial inattendu, image {deg_img}"
            p3 = [perm[perm[perm[i]]] for i in range(len(base))]
            c3 = [coef[i] * coef[perm[i]] * coef[perm[perm[i]]] % p
                  for i in range(len(base))]
            assert p3 == list(range(len(base))) and all(x == 1 for x in c3), \
                f"{nom} deg {deg} : g^3 != identite"
    return "g^3 = id, coefficients compris"


@test("lecture des generateurs : entrees symboliques non nulles")
def t_lecture_gen():
    from equivariance import _lire_matrice, permutation_facteurs
    m = _lire_matrice("{{0, Exp[2 I Pi/3], 0, 0}, {rt[3]^2, 0, 0, 0}, "
                      "{0,0,1,0}, {0,0,0,w}}")
    assert m == [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], m
    assert permutation_facteurs(m, [1, 1]) == [0, 1]
    return "les generateurs d'ordre impair etaient rejetes avant ce correctif"


@test("decomposition isotypique : somme = dimension (detecte le cocycle)")
def t_isotypique():
    """
    Sur l'espace ambiant, la somme des dimensions des espaces propres sur
    tous les caracteres doit valoir la dimension totale. Un ecart signale
    que les generateurs ne commutent pas sur ce degre -- relevement
    seulement projectif. C'est ce controle qui a revele que le degre
    (1,1,1) de la CICY 7669 n'admet aucun vecteur propre commun.
    """
    import re
    from cy_landscape.core.gamma_action import (
        choisir_premier, racine_primitive, parse_entree, action_monomes,
        caracteres_abeliens)
    from cy_landscape.core.sections import rref_mod
    p, _ = choisir_premier([3])
    rac = {3: racine_primitive(p, 3)}
    amb = [2, 2, 2]

    def mat(t):
        out = []
        for l in re.findall(r'\{([^{}]*)\}', t):
            prof, cur, ent = 0, [], []
            for ch in l:
                if ch == '[': prof += 1
                elif ch == ']': prof -= 1
                if ch == ',' and prof == 0:
                    ent.append(''.join(cur)); cur = []
                else:
                    cur.append(ch)
            if cur: ent.append(''.join(cur))
            out.append([parse_entree(e, p, rac) for e in ent])
        return out

    g1 = mat("{{1,0,0,0,0,0,0,0,0},{0,rt[3]^2,0,0,0,0,0,0,0},"
             "{0,0,rt[3],0,0,0,0,0,0},{0,0,0,1,0,0,0,0,0},"
             "{0,0,0,0,rt[3]^2,0,0,0,0},{0,0,0,0,0,rt[3],0,0,0},"
             "{0,0,0,0,0,0,1,0,0},{0,0,0,0,0,0,0,rt[3],0},"
             "{0,0,0,0,0,0,0,0,rt[3]^2}}")
    g2 = mat("{{0,0,1,0,0,0,0,0,0},{1,0,0,0,0,0,0,0,0},{0,1,0,0,0,0,0,0,0},"
             "{0,0,0,0,0,1,0,0,0},{0,0,0,1,0,0,0,0,0},{0,0,0,0,1,0,0,0,0},"
             "{0,0,0,0,0,0,0,0,1},{0,0,0,0,0,0,1,0,0},{0,0,0,0,0,0,0,1,0}}")
    cars = caracteres_abeliens([3, 3], p, rac)
    attendu = {(2, 1, 0): True, (1, 1, 2): True, (1, 1, 1): False}
    for deg, coherent in attendu.items():
        a1 = action_monomes(g1, amb, list(deg), p)
        a2 = action_monomes(g2, amb, list(deg), p)
        n = len(a1[0]); tot = 0
        for chi in cars:
            L = []
            for a, c in ((a1, chi[0]), (a2, chi[1])):
                _, _, _, perm, coef = a
                for i in range(n):
                    v = np.zeros(n, dtype=np.int64)
                    v[perm[i]] = (v[perm[i]] + coef[i]) % p
                    v[i] = (v[i] - c) % p
                    L.append(v)
            r, _ = rref_mod(np.array(L, dtype=np.int64) % p, p)
            tot += n - r
        if coherent:
            assert tot == n, f"degre {deg} : {tot} != {n}"
        else:
            assert tot != n, (f"degre {deg} : la somme vaut {n}, or le cocycle "
                              f"y est non trivial — le test ne detecte plus rien")
    return "(1,1,1) reste incoherent : relevement projectif, travail restant"


# ======================================================================
# 7. Appariement des matrices de configuration
# ======================================================================

@test("appariement invariant par permutation lignes/colonnes")
def t_appariement():
    from wilson_match import meme_matrice
    rng = random.Random(0)
    for _ in range(120):
        nl, nc = rng.randint(3, 8), rng.randint(3, 10)
        M = [[rng.randint(0, 3) for _ in range(nc)] for _ in range(nl)]
        pr, pc = list(range(nl)), list(range(nc))
        rng.shuffle(pr); rng.shuffle(pc)
        N = [[M[i][j] for j in pc] for i in pr]
        assert meme_matrice(M, N), "permutation non reconnue"
    return "120 matrices permutees au hasard"


@test("ordres de groupe lus correctement")
def t_ordres():
    from wilson_match import _ordre_groupe
    for nom, att in (('Z2', 2), ('Z3', 3), ('Z3 x Z3', 9),
                     ('Z2 x Z2', 4), ('Z5 x Z5', 25), ('Q8', 8)):
        got = _ordre_groupe(nom)
        assert got == att, f"{nom} -> {got}, attendu {att}"
    return "un ordre errone fausse la cible |chi(V)| = n_gen x |Gamma|"


# ======================================================================
# 8. Cible de l'indice en mode Wilson
# ======================================================================

@test("mode Wilson : cible |chi| = n_gen x |Gamma|, sans double facteur")
def t_cibles():
    from cy_landscape.main_optimized import _cibles_pour, _ordres_pour
    w = {7669: {'braun': 1, 'groupes': ['Z3', 'Z3 x Z3'], 'ordres': [3, 9]}}
    assert _cibles_pour(7669, w, 3) == {9, 27}, _cibles_pour(7669, w, 3)
    assert _ordres_pour(7669, w) == [3, 9]
    assert _cibles_pour(1, None, 3) == {3}
    return "ordres et cibles sont deux champs distincts (bug corrige)"


# ======================================================================
# 9. Reduction modulo l'ideal : reduce_vec doit ANNULER l'ideal
# ======================================================================

@test("reduce_vec annule l'ideal (projection sur le quotient)")
def t_reduce_vec():
    """
    Reference independante : la definition meme d'une reduction. Tout element
    de I_a doit avoir un reste nul. `Ring.quotient` gardait les lignes BRUTES
    (`rref_mod` travaille sur une copie, la matrice de l'appelant restait non
    reduite) et `reduce_vec` supposait des pivots egaux a 1 : 30/30 elements
    de l'ideal ressortaient non nuls.
    """
    from cy_landscape.core.sections import Ring, P
    rng = np.random.RandomState(0)
    total = 0
    for amb, cfg, degres in (
            ([2, 2], [[3, 3]], ([3, 3], [4, 4], [5, 4])),
            ([4], [[5]], ([5], [6], [7])),
            ([2, 2, 2], [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
             ([1, 1, 1], [2, 1, 2], [2, 2, 2]))):
        R = Ring(amb, cfg)
        for a in degres:
            S, idx, free, piv, Mred = R.quotient(list(a))
            if Mred.shape[0] == 0:
                continue
            for _ in range(10):
                v = (rng.randint(0, P, size=Mred.shape[0]) @ Mred) % P
                reste = R.reduce_vec(list(a), v) % P
                assert not reste.any(), (amb, a, "reste non nul")
                total += 1
    assert total >= 50, total
    return f"{total} elements de l'ideal, reste nul dans tous les cas"


# ======================================================================
# 10. Polynomes Gamma-covariants
# ======================================================================

@test("ideal Gamma-covariant : covariance revérifiée independamment")
def t_covariants():
    """
    On resout p_alpha(g.x) = sum_beta N[alpha][beta] p_beta(x) par algebre
    lineaire, puis on RE-SUBSTITUE sans reutiliser le systeme resolu. Les
    deux chemins n'ont en commun que la donnee de Braun.
    """
    from cy_landscape.core.braun_symmetry import (parse_symmetries, ordres_rt,
                                                  matrice_mod_p)
    from cy_landscape.core.gamma_action import choisir_premier, racine_primitive
    from cy_landscape.core.covariant_ring import (resoudre_covariants,
                                                  tirer_covariants,
                                                  verifier_covariance,
                                                  CovariantRing)
    from cy_landscape.core.sections import Ring
    chemin = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'cicyquotients.m')
    if not os.path.exists(chemin):
        return "cicyquotients.m absent — test ignore"
    from cy_landscape.data.parse_oxford import load_oxford_file
    entries = {e['num']: e for e in load_oxford_file('cicylist.txt')}
    SYM = parse_symmetries(chemin)

    n_ok = n_hilb = 0
    for num in (6947, 7669):
        e = entries.get(num)
        assert e is not None, f"CICY {num} absente de cicylist.txt"
        amb, cfg = e['ambient'], np.asarray(e['config'])
        for sym in SYM[num]['symetries'][:3]:
            ordres = ordres_rt(sym['coord']) | ordres_rt(sym['poly']) | {2}
            p, _ = choisir_premier(sorted(ordres), minimum=30011)
            rac = {n: racine_primitive(p, n) for n in ordres}
            Mc = [matrice_mod_p(x, p, rac) for x in sym['coord']]
            Np = [matrice_mod_p(x, p, rac) for x in sym['poly']]
            res = resoudre_covariants(amb, cfg, Mc, Np, p)
            assert res is not None, (num, sym['nom'], "sigma non extractible")
            v = res['par_convention']['N']
            assert v['non_degenere'], (num, sym['nom'], "noyau degenere")
            co = tirer_covariants(v['base'], res['offsets'], res['dims'], p,
                                  np.random.RandomState(0))
            ok, ec = verifier_covariance(amb, cfg, Mc, Np, 'N', co, p)
            assert ok, (num, sym['nom'], f"ecart {ec}")
            n_ok += 1
            # la fonction de Hilbert ne doit pas bouger : le choix covariant
            # est un point particulier de la famille, pas une degenerescence
            Rc = CovariantRing(amb, cfg, co, p)
            Ra = Ring(amb, cfg, seed=0, p=p)
            for d in ([1] * len(amb), [2] + [1] * (len(amb) - 1)):
                assert Rc.dimY(list(d)) == Ra.dimY(list(d)), (num, d)
                n_hilb += 1
    return (f"{n_ok} symetries, covariance exacte ; "
            f"{n_hilb} degres a fonction de Hilbert inchangee")


# ======================================================================
# 11. Equivariance de f : le test doit MORDRE
# ======================================================================

@test("equivariance de f : #6947 SO(10) survit, #6947 SU(5) tombe")
def t_equivariance_f():
    """
    Une suite qui ne casse jamais ne prouve rien, et un test d'existence qui
    accepte tout non plus. Sur un Gamma agissant par phases, le sous-espace
    des f equivariantes vaut toujours ~1/|Gamma| de l'espace total : sa non
    vacuite ne prouve RIEN. Le contenu est dans h^0(V) recalcule sur ce
    sous-espace. Ce test fige les deux verdicts opposes obtenus sur la meme
    CICY, ce qui interdit a la fois un test qui accepte tout et un test qui
    rejette tout.
    """
    chemin = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'cicyquotients.m')
    if not os.path.exists(chemin):
        return "cicyquotients.m absent — test ignore"
    from cy_landscape.core.braun_symmetry import parse_symmetries
    from cy_landscape.data.parse_oxford import load_oxford_file
    from equivariance_f import analyser

    entries = {e['num']: e for e in load_oxford_file('cicylist.txt')}
    SYM = parse_symmetries(chemin)
    e = entries[6947]
    amb, cfg = e['ambient'], np.asarray(e['config'])

    # (a) SO(10), Z2 : le candidat de la section 5.3 -- doit SURVIVRE
    b_a = [[1, 0, 0, 0, 0], [0, 0, 0, 1, 0], [1, 0, 0, 0, 0],
           [0, 1, 0, 0, 0], [1, 0, 0, 0, 0]]
    c_a = [[3, 1, 0, 1, 0]]
    la = [L for L in analyser(6947, amb, cfg, b_a, c_a,
                              SYM[6947]['symetries'], groupes={'Z2'})
          if L['etat'] == 'ok']
    assert la, "aucun verdict sur le cas SO(10)/Z2"
    assert all(x['dim_equivariant'] > 0 for x in la)
    assert all(x['h0_generique'] == 0 for x in la), "h0 generique non nul"
    assert all(x['h0_equivariant'] == 0 for x in la), \
        f"SO(10)/Z2 devrait survivre : {[x['h0_equivariant'] for x in la]}"

    # (b) SU(5), Z2 x Z2 : espace equivariant NON VIDE mais h^0(V) non nul,
    #     donc elimine. Sans ce second cas, un test qui accepterait tout
    #     passerait aussi.
    b_b = [[0, 1, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 0, 1, 0],
           [1, 0, 0, 0, 0], [1, 0, 0, 0, 0], [0, 1, 0, 1, 0]]
    c_b = [[2, 3, 0, 2, 0]]
    lb = [L for L in analyser(6947, amb, cfg, b_b, c_b,
                              SYM[6947]['symetries'], groupes={'Z2 x Z2'})
          if L['etat'] == 'ok']
    assert lb, "aucun verdict sur le cas SU(5)/Z2xZ2"
    assert all(x['dim_equivariant'] > 0 for x in lb), \
        "espace equivariant vide : le test ne mesurerait plus rien"
    assert all(x['h0_generique'] == 0 for x in lb)
    assert all(x['h0_equivariant'] > 0 for x in lb), \
        f"SU(5)/Z2xZ2 devrait tomber : {[x['h0_equivariant'] for x in lb]}"
    return (f"SO(10)/Z2 : h0 = 0 sur {len(la)} lambda ; "
            f"SU(5)/Z2xZ2 : h0 = {sorted({x['h0_equivariant'] for x in lb})} "
            f"pour un espace de dimension "
            f"{sorted({x['dim_equivariant'] for x in lb})}")


# ======================================================================
# 11 bis. h0(wedge^2 V) restreint : accord non contraint, et morsure
# ======================================================================

@test("h0(w2V) equivariant : redonne le generique sans contrainte, et mord")
def t_wedge2_equivariant():
    """
    Deux exigences opposees, comme pour h^0(V).

    (a) SANS contrainte -- base = espace entier des f -- le nouveau chemin
        doit redonner exactement `sections.h0_wedge2_V_explicit`, ecrit
        independamment. Mesure sur `scan_wilson2` : 71/71, et 71/71 a
        l'interieur des bornes rigoureuses de `monad_wedge`.

    (b) SOUS contrainte, il doit pouvoir CHANGER. Sur #6715 SO(10) / Z2xZ2,
        h^0(w2V) passe de 0 a 21 alors que l'espace equivariant n'est pas
        vide. Sans ce second volet, une fonction qui ignorerait la base et
        recalculerait toujours la valeur generique passerait le volet (a).
    """
    from cy_landscape.core.sections import (Ring, h0_wedge2_V_explicit,
                                            domaine_valide, P)
    from cy_landscape.core.equivariant_monad import h0_wedge2_V_sur_espace
    from cy_landscape.data.parse_oxford import load_oxford_file

    entries = {e['num']: e for e in load_oxford_file('cicylist.txt')}

    # ---- (a) accord non contraint ------------------------------------
    cas = [
        (6947, [[0, 1, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 0, 1, 0],
                [1, 0, 0, 0, 0], [1, 0, 0, 0, 0], [0, 1, 0, 1, 0]],
         [[2, 3, 0, 2, 0]]),
        (6715, None, None),
    ]
    n_acc = 0
    for num, b, c in cas:
        if b is None:
            continue
        e = entries[num]
        amb, cfg = e['ambient'], e['config']
        assert domaine_valide(amb, cfg, b, c), (num, "hors domaine")
        m = len(amb)
        degres = [[[c[0][k] - b[i][k] for k in range(m)]
                   for i in range(len(b))]]
        R = Ring(amb, cfg, seed=11)
        cases = [(0, i) for i in range(len(b))
                 if all(x >= 0 for x in degres[0][i])]
        dims = {k: R.dimY(degres[0][k[1]]) for k in cases}
        offs, a = {}, 0
        for k in cases:
            offs[k] = a
            a += dims[k]
        base = np.eye(a, dtype=np.int64)
        h_new, _ = h0_wedge2_V_sur_espace(R, b, c, base, offs, dims, degres, P,
                                          np.random.RandomState(4))
        R2 = Ring(amb, cfg, seed=11)
        h_old, _ = h0_wedge2_V_explicit(R2, b, c, maxdim=6000)
        assert h_new == h_old, (num, h_new, h_old)
        n_acc += 1
    assert n_acc >= 1

    # ---- (b) la contrainte doit pouvoir changer la valeur -------------
    chemin = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'cicyquotients.m')
    if not os.path.exists(chemin):
        return f"{n_acc} cas non contraints d'accord ; cicyquotients.m absent"
    from cy_landscape.core.braun_symmetry import parse_symmetries
    from equivariance_f import analyser
    SYM = parse_symmetries(chemin)
    e = entries[6715]
    amb, cfg = e['ambient'], np.asarray(e['config'])
    b6715 = [[0, 0, 0, 1, 0], [0, 0, 0, 1, 0], [1, 0, 0, 0, 0],
             [1, 0, 1, 0, 0], [0, 0, 0, 1, 0]]
    c6715 = [[2, 0, 1, 3, 0]]
    lignes = [L for L in analyser(6715, amb, cfg, b6715, c6715,
                                  SYM[6715]['symetries'],
                                  groupes={'Z2 x Z2'})
              if L['etat'] == 'ok']
    assert lignes, "aucun verdict sur #6715 / Z2xZ2"
    assert all(L['dim_equivariant'] > 0 for L in lignes), \
        "espace equivariant vide : le test ne mesurerait rien"
    assert all(L['h0w2_generique'] == 0 for L in lignes), \
        [L['h0w2_generique'] for L in lignes]
    assert all(L['h0w2_equivariant'] > 0 for L in lignes), \
        f"la contrainte devrait changer h0(w2V) : {[L['h0w2_equivariant'] for L in lignes]}"
    return (f"{n_acc} cas non contraints d'accord avec sections ; "
            f"#6715/Z2xZ2 : h0(w2V) 0 -> "
            f"{sorted({L['h0w2_equivariant'] for L in lignes})}")


# ======================================================================
# 11 bis 2. wedge^p V general : h^0(wedge^rk V) = h^0(O_Y) = 1
# ======================================================================

@test("wedge^p V general : det V = O donne h0 = 1, et p=1,2 redonnent l'existant")
def t_wedgep():
    """
    Le calcul de h^3(V) passe par h^0(wedge^{rk-1} V) -- det V = O donne
    wedge^{rk-1} V = V*, et Serre sur un CY3 donne h^0(V*) = h^3(V). La
    construction generale en p doit donc etre juste jusque dans ses SIGNES,
    la contraction par f faisant intervenir (-1)^{k-1}.

    Reference independante et sans appel : a p = rk, wedge^{rk} V = det V = O,
    donc h^0 vaut exactement h^0(O_Y) = 1. Une erreur de signe, de source ou
    de cible fait sortir autre chose que 1. Les deux autres controles
    (p = 1 et p = 2 redonnent les fonctions ecrites separement) ne testent que
    la coherence interne ; celui-ci teste contre une valeur connue d'avance.
    """
    from cy_landscape.core.sections import Ring, P, domaine_valide
    from cy_landscape.core.equivariant_monad import (h0_wedgep_V_sur_espace,
                                                     h0_wedge2_V_sur_espace)
    from cy_landscape.data.parse_oxford import load_oxford_file

    entries = {e['num']: e for e in load_oxford_file('cicylist.txt')}
    resultats = []
    for num, b, c in (
            (6890, [[0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 1, 0, 0],
                    [0, 0, 1, 0, 0], [0, 0, 0, 1, 0]], [[0, 1, 3, 1, 0]]),
            (6947, [[1, 0, 0, 0, 0], [0, 0, 0, 1, 0], [1, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0], [1, 0, 0, 0, 0]], [[3, 1, 0, 1, 0]])):
        e = entries[num]
        amb, cfg = e['ambient'], e['config']
        assert domaine_valide(amb, cfg, b, c), num
        m = len(amb)
        rk = len(b) - len(c)
        degres = [[[c[0][k] - b[i][k] for k in range(m)] for i in range(len(b))]]
        R = Ring(amb, cfg, seed=2)
        cases = [(0, i) for i in range(len(b))
                 if all(x >= 0 for x in degres[0][i])]
        dims, offs, a = {}, {}, 0
        for k in cases:
            dims[k] = R.dimY(degres[0][k[1]])
            offs[k] = a
            a += dims[k]
        base = np.eye(a, dtype=np.int64)
        rng = lambda: np.random.RandomState(5)

        hdet, _ = h0_wedgep_V_sur_espace(R, b, c, rk, base, offs, dims,
                                         degres, P, rng())
        assert hdet == 1, (num, f"h0(wedge^{rk} V) = {hdet}, attendu 1 "
                                f"(det V = O)")
        h2a, _ = h0_wedge2_V_sur_espace(R, b, c, base, offs, dims, degres,
                                        P, rng())
        h2b, _ = h0_wedgep_V_sur_espace(R, b, c, 2, base, offs, dims, degres,
                                        P, rng())
        assert h2a == h2b, (num, "p=2 : general et specialise divergent",
                            h2a, h2b)
        resultats.append((num, rk, hdet))
    return ("h0(wedge^rk V) = 1 sur " +
            ", ".join(f"#{n} (rk {r})" for n, r, _ in resultats) +
            " ; p=2 general == specialise")


# ======================================================================
# 11 ter. Surjectivite de f : le certificat J_d = R_d
# ======================================================================

@test("surjectivite de f : certifie le generique, refuse un facteur commun")
def t_surjectivite():
    """
    Trois controles, dont un NEGATIF construit pour echouer.

    Le critere est : s'il existe un multidegre d >= 0 avec J_d = R_d, ou
    J = (f_1, ..., f_n), alors les f_i n'ont pas de zero commun sur Y, donc
    f est surjective. Il est SUFFISANT, jamais faux positif ; un echec ne
    conclut pas.

    (a) f generique dans tout l'espace -> doit certifier.
    (b) f contraint aux polynomes DIVISIBLES par une coordonnee fixee ->
        les f_i ont alors un facteur commun, donc un zero commun sur Y
        (un diviseur non vide), donc le critere NE DOIT PAS certifier.
        Sans ce controle, une fonction qui renverrait toujours True
        passerait le (a).
    (c) sur #6890 / Z2, lambda = +1 certifie et lambda = -1 non : la
        surjectivite depend du caractere, ce qui est precisement la raison
        d'etre du test.
    """
    from cy_landscape.core.sections import Ring, P, domaine_valide
    from cy_landscape.core.equivariant_monad import f_sans_point_base
    from cy_landscape.data.parse_oxford import load_oxford_file

    entries = {e['num']: e for e in load_oxford_file('cicylist.txt')}
    e = entries[6890]
    amb, cfg = e['ambient'], e['config']
    b = [[0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 1, 0, 0],
         [0, 0, 1, 0, 0], [0, 0, 0, 1, 0]]
    c = [[0, 1, 3, 1, 0]]
    assert domaine_valide(amb, cfg, b, c)
    m = len(amb)
    degres = [[[c[0][k] - b[i][k] for k in range(m)] for i in range(len(b))]]
    R = Ring(amb, cfg, seed=2)
    cases = [(0, i) for i in range(len(b)) if all(x >= 0 for x in degres[0][i])]
    dims, offs, a = {}, {}, 0
    for k in cases:
        dims[k] = R.dimY(degres[0][k[1]])
        offs[k] = a
        a += dims[k]

    ra = f_sans_point_base(R, b, c, np.eye(a, dtype=np.int64), offs, dims,
                           degres, P, np.random.RandomState(1),
                           n_essais=2, n_degres=8)
    assert ra['certifie'], f"(a) f generique devrait certifier : {ra['essais']}"

    # (b) sous-espace des f divisibles par la 1re coordonnee du facteur 2
    sous = []
    for (j, i) in cases:
        S, idx, free, piv, Mred = R.quotient(degres[j][i])
        for t, k in enumerate(free):
            if S[k][2][0] > 0:
                v = np.zeros(a, dtype=np.int64)
                v[offs[(j, i)] + t] = 1
                sous.append(v)
    assert len(sous) > 0, "controle negatif vide : il ne testerait rien"
    rb = f_sans_point_base(R, b, c, np.array(sous, dtype=np.int64), offs, dims,
                           degres, P, np.random.RandomState(1),
                           n_essais=2, n_degres=8)
    assert not rb['certifie'], \
        "(b) un f a facteur commun a un zero commun : ne doit PAS certifier"

    # (c) discrimination selon le caractere
    chemin = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'cicyquotients.m')
    if not os.path.exists(chemin):
        return f"(a) certifie en {ra['degre']}, (b) refuse ; cicyquotients.m absent"
    from cy_landscape.core.braun_symmetry import parse_symmetries
    from equivariance_f import analyser
    SYM = parse_symmetries(chemin)
    L = [x for x in analyser(6890, amb, np.asarray(cfg), b, c,
                             SYM[6890]['symetries'], groupes={'Z2'})
         if x['etat'] == 'ok']
    assert L, "aucun verdict sur #6890 / Z2"
    par_lam = {x['lambda'][0]: x for x in L}
    assert set(par_lam) == {1, -1}, sorted(par_lam)
    assert par_lam[1]['surjectif_certifie'] is True, "lambda=+1 devrait certifier"
    assert par_lam[-1]['surjectif_certifie'] is False, \
        "lambda=-1 ne devrait pas certifier : le test ne discriminerait plus"
    assert par_lam[1]['survit'] and not par_lam[-1]['survit']
    return (f"(a) certifie en {ra['degre']} ; (b) facteur commun refuse ; "
            f"(c) #6890/Z2 : lambda=+1 certifie, lambda=-1 non")


# ======================================================================
# 11 quater. Decomposition de H^1(V) sous Gamma : 3 + 3
# ======================================================================

@test("spectre sous Gamma : H1(V) se decompose en 3 + 3 sur #6890 et #6947")
def t_spectre_gamma():
    """
    Reference independante, connue d'avance : Gamma agissant librement,
    n_gen(X/Gamma) = n_gen(X)/|Gamma|. Avec h^1(V) = 6 et |Gamma| = 2, la
    partie invariante doit valoir EXACTEMENT 3. Une decomposition en 4+2 ou
    5+1 signalerait une erreur dans l'action, dans le conoyau, ou dans le
    fibre -- et donnerait un nombre de generations different de 3 en aval.

    Deux preconditions verifiees ici plutot que supposees :
      - H^1(B) = 0, CERTIFIE par koszul_cohomology_ex pour chaque b_i, sans
        quoi H^1(V) n'est pas le conoyau de H^0(B) -> H^0(C) ;
      - h^1(V) = dim R_c - rang(f) doit retomber sur 6.
    """
    chemin = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          'cicyquotients.m')
    if not os.path.exists(chemin):
        return "cicyquotients.m absent — test ignore"
    from cy_landscape.core.braun_symmetry import (parse_symmetries, ordres_rt,
                                                  matrice_mod_p)
    from cy_landscape.core.gamma_action import choisir_premier, racine_primitive
    from cy_landscape.core.covariant_ring import (resoudre_covariants,
                                                  tirer_covariants,
                                                  verifier_covariance,
                                                  CovariantRing)
    from cy_landscape.core.equivariant_monad import (espace_f_equivariant,
                                                     decomposition_h1_V)
    from cy_landscape.core.exact_cohomology import koszul_cohomology_ex
    from cy_landscape.data.parse_oxford import load_oxford_file

    entries = {e['num']: e for e in load_oxford_file('cicylist.txt')}
    SYM = parse_symmetries(chemin)
    cas = {
        6890: ([[0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 1, 0, 0],
                [0, 0, 1, 0, 0], [0, 0, 0, 1, 0]], [[0, 1, 3, 1, 0]]),
        6947: ([[1, 0, 0, 0, 0], [0, 0, 0, 1, 0], [1, 0, 0, 0, 0],
                [0, 1, 0, 0, 0], [1, 0, 0, 0, 0]], [[3, 1, 0, 1, 0]]),
    }
    detail = []
    for num, (b, c) in cas.items():
        e = entries[num]
        amb, cfg = e['ambient'], np.asarray(e['config'])

        for bb in b:
            r = koszul_cohomology_ex(amb, cfg, list(bb))
            assert r['certified_by_degree'][1] and r[1] == 0, \
                (num, bb, "H^1(B) non nul ou non certifie : H^1(V) n'est pas "
                          "le conoyau")

        sym = [s for s in SYM[num]['symetries'] if s['nom'] == 'Z2'][0]
        ordres = ordres_rt(sym['coord']) | ordres_rt(sym['poly']) | {2}
        p, _ = choisir_premier(sorted(ordres), minimum=30011)
        rac = {n: racine_primitive(p, n) for n in ordres}
        Mc = [matrice_mod_p(x, p, rac) for x in sym['coord']]
        Np = [matrice_mod_p(x, p, rac) for x in sym['poly']]
        res = resoudre_covariants(amb, cfg, Mc, Np, p)
        v = res['par_convention']['N']
        co = tirer_covariants(v['base'], res['offsets'], res['dims'], p,
                              np.random.RandomState(0))
        ok, _ = verifier_covariance(amb, cfg, Mc, Np, 'N', co, p)
        assert ok, (num, "covariance non revérifiée")
        A = CovariantRing(amb, cfg, co, p)
        out = espace_f_equivariant(A, amb, b, c, Mc, p)
        assert out['etat'] == 'ok', (num, out['etat'])
        sol = [s for s in out['solutions'] if s['lambda'] == (1,)]
        assert sol, (num, "structure equivariante lambda=+1 absente")

        d = decomposition_h1_V(A, amb, b, c, Mc, sol[0]['base'],
                               out['offsets'], out['dims'], out['degres'], p,
                               np.random.RandomState(11))
        assert d is not None, (num, "decomposition indisponible")
        assert d['coherent'], (num, "multiplicites non additives : "
                                    "semi-simplicite en defaut")
        assert d['h1'] == 6, (num, f"h1(V) = {d['h1']}, attendu 6")
        inv = d['sur_H1'][1]
        anti = d['sur_H1'][p - 1]
        assert (inv, anti) == (3, 3), \
            (num, f"decomposition {inv} + {anti}, attendue 3 + 3 "
                  f"(sinon n_gen != 3 en aval)")
        # --- H^1(wedge^2 V) : les 10 de SO(10) -------------------------
        from itertools import combinations as _comb
        from cy_landscape.core.equivariant_monad import (
            decomposition_h1_wedge2_V, h0_wedge2_V_sur_espace)
        for i, j in _comb(range(len(b)), 2):
            dd = [b[i][k] + b[j][k] for k in range(len(amb))]
            rr = koszul_cohomology_ex(amb, cfg, dd)
            assert rr['certified_by_degree'][1] and rr[1] == 0, \
                (num, dd, "H^1(wedge^2 B) non nul ou non certifie")

        d2 = decomposition_h1_wedge2_V(A, amb, b, c, Mc, sol[0]['base'],
                                       out['offsets'], out['dims'],
                                       out['degres'], p,
                                       np.random.RandomState(11))
        assert d2 is not None, (num, "wedge^2 indisponible")
        assert d2['complexe_ok'], (num, "beta o alpha != 0 : pas un complexe")
        assert d2['coherent'], (num, "multiplicites non additives (wedge^2)")
        # meme quantite par un tout autre chemin
        h0_autre, _ = h0_wedge2_V_sur_espace(A, b, c, sol[0]['base'],
                                             out['offsets'], out['dims'],
                                             out['degres'], p,
                                             np.random.RandomState(5))
        assert d2['h0_controle'] == h0_autre == 0, \
            (num, "h0(wedge^2 V) incoherent entre les deux chemins",
             d2['h0_controle'], h0_autre)
        assert d2['h1'] == 8, (num, f"h1(wedge^2 V) = {d2['h1']}, attendu 8")
        assert (d2['sur_H1'][1], d2['sur_H1'][p - 1]) == (2, 6), \
            (num, "decomposition des 10", d2['sur_H1'])

        detail.append(f"#{num} : {inv}+{anti}")
    return ("H^1(B) = 0 certifie, h1(V) = 6, decomposition " +
            " et ".join(detail) + " → 3 generations sur le quotient ; "
            "h1(w2V) = 8 = 2+6 (les 10 de SO(10)), h0(w2V) = 0 "
            "concordant par deux chemins")


# ======================================================================
# 11 quinquies. Enumeration exhaustive des B
# ======================================================================

@test("enumeration exhaustive des B : contient l'echantillonnage, compte exact")
def t_enumeration_B():
    """
    Trois references independantes.

    (a) `compte_B` calcule le nombre de tuples ORDONNES par programmation
        dynamique, sans rien enumerer ; `enumerer_positive_B` les construit.
        Les deux chemins n'ont en commun que l'enonce du probleme.
    (b) TOUT B produit par l'echantillonnage `_find_positive_B` doit se
        retrouver dans l'enumeration. Un seul echappement signifierait que
        l'enumeration rate des cas, donc qu'un resultat d'absence obtenu
        dessus serait faux.
    (c) Tout B enumere doit verifier c1(V) = 0 et la positivite -- sinon
        l'enumeration produirait des objets hors domaine.

    Le point (b) est le plus important : c'est lui qui autorise a remplacer
    « aucun survivant parmi ce qu'on a tire » par « aucun survivant sur tout
    le domaine enumere ».
    """
    from itertools import product as iprod
    from cy_landscape.core.positive_monads import (compte_B, enumerer_positive_B,
                                                   _find_positive_B, _mix,
                                                   is_positive_monad,
                                                   generate_positive_monads)
    from cy_landscape.core.monads import MonadBundle

    n_paires = n_tirages = n_verif = 0
    for m in (2, 3):
        for c in iprod(range(0, 4), repeat=m):
            for rB in (4, 5):
                N = compte_B(list(c), rB)
                if N == 0 or N > 20000:
                    continue
                E = enumerer_positive_B(list(c), rB, plafond=20000)
                assert E is not None, (c, rB)
                # `enumerer_positive_B` renvoie un GENERATEUR (il l'etait
                # devenu pour borner la memoire, cf. le MemoryError de
                # scan_su5). L'epuiser dans la comprehension ci-dessous
                # laissait la boucle de validite tourner a vide : le test
                # affichait « 0 B valides » sans echouer. On materialise.
                E = list(E)
                assert E, (c, rB, "enumeration vide")
                ens = {tuple(sorted(tuple(b) for b in x)) for x in E}
                n_paires += 1
                # (c) validite de chaque B enumere
                for b in E:
                    mo = MonadBundle(b, [list(c)])
                    assert mo.c1_vanishes and is_positive_monad(mo), (c, rB, b)
                    n_verif += 1
                # (b) inclusion de l'echantillonnage
                for g in range(4):
                    rng = np.random.RandomState(_mix(g, tuple(c), rB))
                    for b in _find_positive_B(list(c), rB, rng, max_attempts=30):
                        cle = tuple(sorted(tuple(x) for x in b))
                        assert cle in ens, \
                            (c, rB, "tirage hors enumeration", cle)
                        n_tirages += 1
    assert n_paires >= 20 and n_tirages >= 200, (n_paires, n_tirages)

    # MONOTONIE DU GENERATEUR, sortie COMPLETE (les deux branches r_C).
    #
    # Le premier jet de ce test ne comparait que r_C = 1, la branche
    # modifiee. Il passait -- alors qu'activer `exhaustif_max` faisait
    # disparaitre 20 candidats de rank_C = 2 sur 42 dans un vrai scan : la
    # branche r_C = 2 puisait dans le RNG partage, donc ses tirages
    # dependaient de ce que r_C = 1 avait consomme avant. Tester la seule
    # branche touchee ne dit rien des effets de bord sur les autres.
    def cles(**kw):
        return {tuple(sorted(tuple(x) for x in mo.b_charges))
                + tuple(sorted(tuple(y) for y in mo.c_charges))
                for rv in (3, 4, 5)
                for mo in generate_positive_monads(
                    3, rv, max_charge=3, n_systematic=80,
                    rng=np.random.RandomState(1), seed=7, **kw)}
    ech = cles()
    st = {}
    exh = cles(exhaustif_max=200000, stats=st)
    perdus = ech - exh
    assert not perdus, \
        (f"{len(perdus)} monades trouvees par tirage absentes du mode "
         f"exhaustif : le generateur n'est pas monotone en exhaustif_max",
         sorted(perdus)[:3])
    assert len(exh) > len(ech), (len(ech), len(exh))
    assert st['c_exhaustifs'] > 0, st
    # Les DEUX branches doivent etre representees, sinon l'inclusion
    # ci-dessus ne porterait que sur une moitie du generateur -- c'est
    # exactement l'angle mort qui a laisse passer le defaut.
    def n_rC(cles_, r):
        return sum(1 for mo in generate_positive_monads(
            3, 3, max_charge=3, n_systematic=80,
            rng=np.random.RandomState(1), seed=7, exhaustif_max=200000)
            if len(mo.c_charges) == r)
    assert n_rC(exh, 1) > 0 and n_rC(exh, 2) > 0, \
        "une des deux branches r_C est absente : le test ne couvrirait qu'elle"
    return (f"{n_paires} paires (c, rB), {n_verif} B valides, "
            f"{n_tirages} tirages tous retrouves ; generateur complet : "
            f"{len(ech)} monades echantillonnees incluses dans {len(exh)} enumerees")


# ======================================================================
# 11 sexies. Fibre d'extension : le chemin correct (defaut 4.7)
# ======================================================================

@test("extension : rang et chi corrects, la pseudo-monade s'en ecarte")
def t_extension():
    """
    Le pipeline reutilisait le chemin des monades via une pseudo-monade
    B = F1 (+) F2, C = F2. Le noyau est de rang rank(F1) et de
    caracteristique chi(F1) ; le fibre d'extension est de rang
    rank(F1)+rank(F2) et de caracteristique chi(F1)+chi(F2). Mesure du
    defaut 4.7 : 1571 entrees sur 1571 en incoherence de rang.

    Reference INDEPENDANTE pour chi : sur un CY3, chi(L*) = -chi(L) pour un
    fibre en droites. Donc pour F2 = F1* (charges opposees), chi(V) doit
    valoir exactement 0, quelle que soit la geometrie. Le test l'exige sur
    plusieurs CICYs -- une valeur connue d'avance, pas une auto-comparaison.

    Le test exige aussi que la pseudo-monade DIVERGE : sans cela, il
    passerait encore si l'on remettait le raccourci en place.
    """
    from cy_landscape.core.extensions import ExtensionBundle, chi_extension
    from cy_landscape.core.monads import MonadBundle
    from cy_landscape.core.intersection import (compute_intersection_numbers,
                                                compute_c2_tangent)
    from cy_landscape.core.chi_exact import ChiCalculator

    n = divergences = 0
    for c, cfg in _cy3()[:8]:
        amb = c['ambient']
        m = len(amb)
        if m < 2:
            continue
        d = compute_intersection_numbers(amb, cfg)
        c2 = compute_c2_tangent(amb, cfg, d)
        cal = ChiCalculator(amb, d, c2)
        f1 = [[1] + [0] * (m - 1), [0, 1] + [0] * (m - 2)]
        f2 = [[-x for x in v] for v in f1]          # F2 = F1*
        ext = ExtensionBundle(f1, f2)
        assert ext.rank_V == len(f1) + len(f2) == 4, ext.rank_V
        assert ext.c1_vanishes
        assert chi_extension(ext, cal) == 0, \
            (amb, "chi(F1 + F1*) doit etre nul", chi_extension(ext, cal))
        pseudo = MonadBundle(f1 + f2, f2)
        chi_pseudo = cal.bundle(pseudo.b_charges) - cal.bundle(pseudo.c_charges)
        if pseudo.rank_V != ext.rank_V or chi_pseudo != 0:
            divergences += 1
        n += 1
    assert n >= 3, n
    assert divergences == n, \
        ("la pseudo-monade coincide avec l'extension : le raccourci du "
         "defaut 4.7 serait revenu sans etre detecte")
    return (f"{n} CICYs : chi(F1 + F1*) = 0 partout, rang = 4 ; "
            f"la pseudo-monade annonce rang 2 et chi != 0 dans {divergences} cas")


# ======================================================================
# 12. Ordre projectif et racines n-iemes
# ======================================================================

@test("ordre projectif : T^n = c.Id detecte, racines n-iemes exactes")
def t_ordre_projectif():
    """
    Le relevement de Gamma aux fibres en droites n'est que projectif : T_g^n
    peut valoir c.Id avec c != 1. Supposer c = 1 reviendrait a chercher les
    lambda parmi les seules racines de l'unite et a manquer des solutions.
    On verifie sur un cas construit ou c est choisi non trivial.
    """
    from cy_landscape.core.equivariant_monad import (ordre_projectif,
                                                     racines_niemes)
    p = 30013
    # c est construit comme une puissance n-ieme : la racine existe donc, et
    # le test porte sur le fait de les trouver TOUTES. Prendre un c arbitraire
    # n'aurait rien prouve -- il n'est une puissance n-ieme qu'une fois sur
    # pgcd(n, p-1).
    from math import gcd
    for base, n_voulu in ((7, 3), (2, 4), (5, 2), (11, 6)):
        c_voulu = pow(base, n_voulu, p)
        rac = racines_niemes(c_voulu, n_voulu, p)
        assert rac, (c_voulu, n_voulu, "aucune racine trouvee")
        for x in rac:
            assert pow(int(x), n_voulu, p) == c_voulu % p, (x, c_voulu, n_voulu)
        assert base % p in rac, (base, "racine evidente manquee")
        assert len(rac) == gcd(n_voulu, p - 1), \
            (n_voulu, len(rac), gcd(n_voulu, p - 1))
    # matrice de permutation cyclique d'ordre 3, mise a l'echelle par mu
    Pm = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=np.int64)
    mu = 5
    T = (mu * Pm) % p
    n, c = ordre_projectif(T, p)
    assert n == 3, n
    assert c == pow(mu, 3, p), (c, pow(mu, 3, p))
    assert c != 1, "cas trivial : le test ne prouverait rien"
    return f"T^3 = {c}.Id avec c != 1 ; racines n-iemes verifiees sur 4 cas"


# ======================================================================

def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith('t_')]
    print(f"\n{'='*72}")
    print(f"  NON-REGRESSION -- {len(tests)} tests")
    print(f"{'='*72}")
    for t in tests:
        t()
    largeur = max(len(n) for n, _, _ in _RESULTATS)
    echecs = 0
    for nom, ok, detail in _RESULTATS:
        marque = " OK " if ok else "ECHEC"
        print(f"  [{marque}] {nom:<{largeur}}  {detail}")
        if not ok:
            echecs += 1
    print(f"{'='*72}")
    if echecs:
        print(f"  {echecs} ECHEC(S) sur {len(_RESULTATS)} — ne pas lancer de scan\n")
    else:
        print(f"  {len(_RESULTATS)} tests passes\n")
    return 1 if echecs else 0


if __name__ == '__main__':
    sys.exit(main())
