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
# 11 septies. Enumeration des extensions : monotonie en max_charge
# ======================================================================

@test("extensions enumerees : monotones en max_charge, la ou le tirage perd 97 %")
def t_enumeration_extensions():
    """
    Le defaut mesure (§5.11) : `generate_extensions` echantillonnait, et
    passer de max_charge 2 a 3 FAISAIT PERDRE 216 extensions sur 222, soit
    97 %. Un generateur non monotone en son propre parametre de portee
    interdit tout enonce d'ABSENCE -- « aucun fibre stable a max_charge 3 »
    ne dit rien tant qu'on ignore si ce domaine contient celui de
    max_charge 2.

    ----------------------------------------------------------------------
    Quatre references, dont une NEGATIVE
    ----------------------------------------------------------------------
    (a) COMPTAGE. `compte_extensions` compte les tuples ordonnes par
        convolution, sans construire un seul tuple ; on le confronte a un
        comptage par enumeration explicite ecrit ici, dans le test. Les
        deux chemins n'ont en commun que l'enonce du domaine.

    (b) MONOTONIE, la propriete visee. E(q) doit etre inclus dans E(q+1),
        composant par composant du domaine. C'est vrai PAR CONSTRUCTION
        pour une enumeration de boites emboitees -- le test verifie que
        l'implementation realise bien cette construction.

    (c) INCLUSION DU TIRAGE. Tout ce que l'ancien chemin produisait doit
        se retrouver dans l'enumeration, sinon remplacer l'un par l'autre
        PERDRAIT des candidats. C'est ce point qui autorise a passer de
        « aucun survivant parmi ce qu'on a tire » a « aucun survivant sur
        le domaine ».

    (d) CONTROLE NEGATIF. Le tirage, lui, doit etre VU non monotone. Sans
        ce volet, un test qui accepterait n'importe quel generateur --
        y compris l'ancien -- passerait. Les deux exigences ensemble
        n'admettent qu'un generateur qui enumere reellement.

    Un cinquieme controle, trivial mais indispensable : l'enumeration ne
    doit pas etre vide et doit croitre strictement. L'ensemble vide est
    monotone et inclus dans tout ; il passerait (a), (b) et (c).
    """
    from itertools import product as iprod
    from cy_landscape.core.extensions import (
        compte_extensions, enumerer_extensions, _splits,
        _echantillonner_extensions)

    def cle(e):
        return (tuple(sorted(tuple(v) for v in e.f1_charges)),
                tuple(sorted(tuple(v) for v in e.f2_charges)))

    def brut_ordonnes(m, rv, q):
        """Comptage par enumeration explicite : reference independante."""
        n = 0
        boite = list(iprod(range(-q, q + 1), repeat=m))
        for _rk1, _rk2 in _splits(rv):
            for libres in iprod(boite, repeat=rv - 1):
                last = [-sum(v[k] for v in libres) for k in range(m)]
                if all(abs(x) <= q for x in last):
                    n += 1
        return n

    def enum_valide(m, rv, q):
        """Enumere en VALIDANT chaque objet : c1(V) = 0, rang, boite."""
        out = set()
        for e in enumerer_extensions(m, rv, q, plafond=10 ** 8):
            assert e.c1_vanishes, (m, rv, q, "c1(V) != 0", e)
            assert e.rank_V == rv, (m, rv, q, "rang", e.rank_V)
            assert all(abs(x) <= q for v in e.f1_charges + e.f2_charges
                       for x in v), (m, rv, q, "hors de la boite", e)
            out.add(cle(e))
        return out

    # (a) le compteur par convolution contre un comptage par enumeration
    n_comptages = 0
    for m in (1, 2):
        for rv in (3, 4, 5):
            for q in ((1, 2) if rv <= 4 else (1,)):
                a = compte_extensions(m, rv, q)
                b = brut_ordonnes(m, rv, q)
                assert a == b, (m, rv, q, "compte_extensions", a, "brut", b)
                assert a > 0, (m, rv, q, "domaine vide")
                n_comptages += 1

    # (b) monotonie en max_charge, et croissance stricte
    n_paires = n_objets = 0
    for m in (1, 2):
        for rv in (3, 4, 5):
            qs = (1, 2, 3) if rv <= 4 else (1, 2)
            E = {q: enum_valide(m, rv, q) for q in qs}
            for q in qs:
                assert E[q], (m, rv, q, "enumeration vide")
                n_objets += len(E[q])
            for q1, q2 in zip(qs, qs[1:]):
                perdus = E[q1] - E[q2]
                assert not perdus, (
                    f"m={m} rk={rv} : {len(perdus)} extensions de max_charge "
                    f"{q1} absentes de max_charge {q2} -- l'enumeration n'est "
                    f"pas monotone", sorted(perdus)[:2])
                assert len(E[q2]) > len(E[q1]), (m, rv, q1, q2,
                                                 len(E[q1]), len(E[q2]))
                n_paires += 1

    # (c) tout tirage doit se retrouver dans l'enumeration
    # (d) et le tirage, lui, doit PERDRE des extensions en passant de 2 a 3
    n_tirages = perdus_tirage = total_tirage = 0
    # (m, rank_V, q_bas, q_haut) -- m = 3 est teste a charge plus basse,
    # l'enumeration y coutant (2q+1)^3 par vecteur.
    for m, rv, qa, qb in ((2, 3, 2, 3), (2, 4, 2, 3), (3, 3, 1, 2)):
        S = {}
        for q in (qa, qb):
            S[q] = {cle(e) for e in
                    _echantillonner_extensions(m, rv, q, n_random=400)}
            E = enum_valide(m, rv, q)
            assert S[q], (m, rv, q, "tirage vide")
            hors = S[q] - E
            assert not hors, (
                f"m={m} rk={rv} q={q} : {len(hors)} tirages hors de "
                f"l'enumeration -- l'enumeration raterait des candidats",
                sorted(hors)[:2])
            n_tirages += len(S[q])
        perdus_tirage += len(S[qa] - S[qb])
        total_tirage += len(S[qa])
    assert perdus_tirage > 0.5 * total_tirage, (
        perdus_tirage, total_tirage,
        "le tirage ne perd plus rien entre max_charge 2 et 3 : le controle "
        "negatif ne mord plus, ce test passerait pour un generateur qui "
        "echantillonne")

    return (f"{n_comptages} comptages convolution == enumeration ; "
            f"{n_paires} paires (q, q+1) sans perte ; {n_tirages} tirages "
            f"tous retrouves ; le tirage, lui, perd {perdus_tirage}/"
            f"{total_tirage} extensions d'un cran de max_charge au suivant")


# ======================================================================
# 11 octies. Pente : les sous-faisceaux destabilisants d'une extension
# ======================================================================

@test("pente : certificat exact, et un echec de recherche n'elimine jamais")
def t_pente_extension():
    """
    Le critere de Hoppe « c1(V) = 0 => V stable <=> h0(w^p V) = 0 » suppose
    Pic(X) de rang 1. Sur une CICY a m > 1 il reste NECESSAIRE sans etre
    suffisant. Une extension, elle, exhibe ses sous-faisceaux -- F1, ses
    sous-sommes, les preimages des sous-sommes de F2 -- donc la pente les
    teste directement.

    ----------------------------------------------------------------------
    Le piege que ce test fige, et qui est la raison principale de son
    existence
    ----------------------------------------------------------------------
    La premiere version du critere concluait « instable » quand la
    recherche d'un temoin J echouait sur une grille [1,4]^m. Elle annoncait
    635 extensions destabilisees sur 2 647. Ce chiffre ne mesurait que la
    grille :

        J_max        3      6      12     24
        sans temoin  1748   1299   1042   925

    Aucune saturation -- exactement le faux lieu de base du §5.4, ou la
    mesure manquante etait la dimension de la source. Le volet (d) exige
    donc qu'un budget insuffisant rende `None` et jamais `False`, sur un
    cas REEL ou un temoin existe mais hors de portee d'une petite grille.

    Quatre references :

    (a) VALEUR CONNUE D'AVANCE. Le degre au point J = v vaut
        sum_ijk d_ijk v_i v_j v_k, que Riemann-Roch relie a chi(O(v)) par
        24*chi = 4*cube + 2*(c2.v). On confronte l'einsum a
        `ChiCalculator`, qui n'emprunte pas le meme chemin.

    (b) CONTROLE POSITIF. Un sous-faisceau a c1 <= 0 non nul doit trouver
        un temoin. Sans ce volet, une fonction refusant tout passerait.

    (c) CONTROLE NEGATIF CONSTRUIT. Un sous-faisceau a c1 >= 0 non nul doit
        etre certifie instable, sur toute CICY. Sans ce volet, une fonction
        acceptant tout passerait.

    (d) NON-ELIMINATION SUR ECHEC DE RECHERCHE. Voir ci-dessus.
    """
    from cy_landscape.core.extensions import (
        ExtensionBundle, pente_extension, certificat_instabilite,
        _sous_faisceaux, degre, ContextePente)
    from cy_landscape.core.intersection import (compute_intersection_numbers,
                                                compute_c2_tangent)
    from cy_landscape.core.chi_exact import ChiCalculator

    # (a) degre au point J = v contre Riemann-Roch
    n_rr = 0
    for c, cfg in _cy3()[:6]:
        amb = c['ambient']
        m = len(amb)
        d = compute_intersection_numbers(amb, cfg)
        c2 = compute_c2_tangent(amb, cfg, d)
        cal = ChiCalculator(amb, d, c2)
        for v in ([1] + [0] * (m - 1), [1] * m, [2] + [-1] * (m - 1)):
            cube = degre(d, v, v)
            c2v = sum(int(round(2 * float(x))) * y for x, y in zip(c2, v))
            assert 4 * cube + c2v == 24 * cal.line(v), \
                (amb, v, cube, cal.line(v))
            n_rr += 1
    assert n_rr >= 15, n_rr

    # (b) et (c) : deux verdicts OPPOSES construits sur les memes CICYs
    n_pos = n_neg = 0
    for c, cfg in _cy3()[:6]:
        amb = c['ambient']
        m = len(amb)
        if m < 2:
            continue
        d = compute_intersection_numbers(amb, cfg)
        ctx = ContextePente(d, m)

        # (c) F1 = O(a) avec a >= 0 non nul : deg_J(a) = sum_i a_i D_i(J)
        #     avec D_i(J) >= 0, donc >= 0 pour TOUTE classe de Kahler.
        a_pos = [1] + [0] * (m - 1)
        a_neg = [-x for x in a_pos]
        ext_bad = ExtensionBundle([a_pos], [a_neg])
        r_bad = pente_extension(ext_bad, ctx=ctx)
        assert r_bad['stable_possible'] is False, (amb, r_bad['etat'])
        assert r_bad['certificat'] is not None, amb
        n_neg += 1

        # (b) le meme, retourne : F1 = O(-a). Aucun certificat ne doit
        #     s'appliquer, et un temoin doit exister.
        ext_ok = ExtensionBundle([a_neg], [a_pos])
        assert certificat_instabilite(_sous_faisceaux(ext_ok)) is None, amb
        r_ok = pente_extension(ext_ok, ctx=ctx)
        assert r_ok['stable_possible'] is True, (amb, r_ok['etat'])
        assert r_ok['temoin'] is not None, amb
        # le temoin doit REELLEMENT verifier la condition
        for v in _sous_faisceaux(ext_ok):
            assert degre(d, v, r_ok['temoin']) < 0, (amb, v, r_ok['temoin'])
        n_pos += 1
    assert n_pos >= 3 and n_neg >= 3, (n_pos, n_neg)

    # (d) LE PIEGE. Cas reel : CICY 14, un temoin existe (J = [1,22,20]),
    #     hors de portee d'une grille [1,2]^3. Le verdict a budget
    #     insuffisant doit etre None -- surtout pas False.
    par_num = {c['num']: (c, cfg) for c, cfg in _cy3()}
    assert 14 in par_num, "CICY 14 absente des donnees"
    c14, cfg14 = par_num[14]
    d14 = compute_intersection_numbers(c14['ambient'], cfg14)
    ext14 = ExtensionBundle([[-2, -1, 2]], [[1, -1, 0], [1, 2, -2]])
    etroit = pente_extension(ext14, d14, J_max=2)
    large = pente_extension(ext14, d14, J_max=24)
    assert large['stable_possible'] is True, large['etat']
    assert etroit['stable_possible'] is None, (
        "un echec de recherche de temoin a ete converti en elimination : "
        "c'est l'erreur du §5.4, et elle annoncait 635 fibres destabilises "
        "qui n'existaient pas", etroit['etat'])
    assert etroit['certificat'] is None and large['certificat'] is None

    return (f"{n_rr} degres au point J=v == Riemann-Roch ; "
            f"{n_neg} certificats exacts sur c1 >= 0, {n_pos} temoins sur "
            f"c1 <= 0 ; sur #14 un temoin hors grille rend `None`, jamais "
            f"`False` (temoin large : {large['temoin']})")


# ======================================================================
# 11 nonies. Hoppe sous sa forme SUFFISANTE : les twists
# ======================================================================

@test("Hoppe suffisant : 110 twists sur #6890, tous a source non vide")
def t_hoppe_suffisant():
    """
    « c1(V) = 0 => V stable <=> h0(wedge^p V) = 0 » est une EQUIVALENCE
    seulement si Pic(X) est de rang 1. #6890 et #6947 ont h11 = 5 : leur
    verdict `stable` ne valait que « non elimine ». La forme suffisante
    ajoute les torsions -- h0(wedge^p V(-H)) = 0 pour tout H de degre
    deg_J(H) >= 0 -- et cet ensemble est fini parce que H est borne au
    dessus par les charges de wedge^p B et en dessous par deg_J(H) >= 0
    avec tous les D_k(J) > 0.

    ----------------------------------------------------------------------
    Le risque propre a ce test : un vert qui ne prouve rien
    ----------------------------------------------------------------------
    Un critere SUFFISANT verifie sur un ensemble vide est vrai sans rien
    demontrer. Si les twists avaient tous une source vide, h0 vaudrait 0
    sans qu'aucun rang ne soit calcule, et « stable » serait un artefact du
    decoupage -- la meme faute que les exotiques structurellement nuls du
    §4.8. Le volet (b) l'interdit explicitement.

    Quatre references :

    (a) VALEUR CONNUE D'AVANCE, negative. On construit une monade dont la
        premiere colonne de f est nulle : c - b_1 a une composante negative,
        donc f_1 = 0, donc O(b_1) est dans le noyau et
        h0(V) >= dim H0(O(b_1)). Le calcul doit rendre EXACTEMENT cette
        valeur, et les deux criteres doivent refuser.

    (b) ANTI-VACUITE. Sur #6890, les 110 twists doivent avoir une source
        NON VIDE -- le verdict porte sur 110 calculs de rang reels.

    (c) STRUCTURE DU POLYTOPE. H = 0 doit y figurer a tout p, sans quoi le
        critere suffisant ne CONTIENDRAIT pas l'ancien. Et un D a
        composante nulle doit rendre None : la borne inferieure disparait,
        l'ensemble n'est plus fini, on ne conclut pas.

    (d) CONTENANCE. Tout ce que Hoppe nu elimine, la forme suffisante doit
        l'eliminer aussi.
    """
    from cy_landscape.core.sections import Ring, P
    from cy_landscape.core.equivariant_monad import (hoppe_sur_espace,
                                                     hoppe_suffisant_sur_espace)
    from cy_landscape.core.hoppe_fast import vecteur_D, polytope_twists
    from cy_landscape.core.intersection import compute_intersection_numbers
    from cy_landscape.data.parse_oxford import load_oxford_file

    entries = {e['num']: e for e in load_oxford_file('cicylist.txt')}
    e = entries[6890]
    amb, cfg = e['ambient'], e['config']
    m = len(amb)
    D = [int(x) for x in vecteur_D(compute_intersection_numbers(amb, cfg),
                                   [1] * m)]
    assert all(x > 0 for x in D), ("D_k(1..1) doit etre > 0 : c'est lui qui "
                                   "borne le polytope", D)

    def monte(b, c):
        R = Ring(amb, cfg, seed=2)
        degres = [[[c[0][k] - b[i][k] for k in range(m)]
                   for i in range(len(b))]]
        cases = [(0, i) for i in range(len(b))
                 if all(x >= 0 for x in degres[0][i])]
        dims, offs, a = {}, {}, 0
        for k in cases:
            dims[k] = R.dimY(degres[0][k[1]])
            offs[k] = a
            a += dims[k]
        return R, np.eye(a, dtype=np.int64), offs, dims, degres

    # (c) structure du polytope
    B = [[0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 1, 0, 0],
         [0, 0, 1, 0, 0], [0, 0, 0, 1, 0]]
    C = [[0, 1, 3, 1, 0]]
    for q in (1, 2, 3):
        Hs = polytope_twists(B, q, D)
        assert Hs, (q, "polytope vide")
        assert [0] * m in Hs, (q, "H = 0 absent du polytope : le critere "
                                  "suffisant ne contiendrait pas l'ancien")
    D0 = list(D)
    D0[0] = 0
    assert polytope_twists(B, 1, D0) is None,         "un D_k nul doit rendre None : le polytope n'est plus borne"

    # (b) et positif : #6890 passe, et sur des sources non vides
    R, base, offs, dims, degres = monte(B, C)
    suf = hoppe_suffisant_sur_espace(R, B, C, base, offs, dims, degres, P,
                                     np.random.RandomState(5), D)
    assert suf['stable'] is True, (suf['stable'], suf['motif'])
    n_tw = suf['n_twists']
    assert n_tw >= 100, n_tw
    assert suf['sources_non_vides'] == n_tw, (
        f"{n_tw - suf['sources_non_vides']} twists a source VIDE : le "
        f"verdict porterait sur un ensemble sans contenu", suf['detail'])
    # LE TWIST DOIT AGIR. Si le parametre etait ignore, chaque twist
    # recalculerait H = 0 et toutes les sources seraient egales : le test
    # resterait vert alors que le critere ne serait plus que l'ancien.
    # Sur #6890 les h0 valent 0 partout, donc c'est la SEULE trace visible.
    agit = [q for q, v in suf['detail'].items()
            if v['source_H0'] is not None and v['source_max'] > v['source_H0']]
    assert agit, (
        "aucun twist n'agrandit la source par rapport a H = 0 : le "
        "parametre `twist` est probablement ignore, et le critere retombe "
        "sur l'ancien sans que rien ne le signale",
        {q: (v['source_H0'], v['source_max']) for q, v in suf['detail'].items()})

    # (a) valeur connue d'avance, negative : f_1 = 0 met O(b_1) dans le noyau
    Bd = [[0, 2, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 1, 0, 0],
          [0, 0, 1, 0, 0], [0, -1, 0, 1, 0]]
    assert [sum(x[k] for x in Bd) for k in range(m)] == C[0], "c1(V) != 0"
    R2, base2, offs2, dims2, degres2 = monte(Bd, C)
    attendu = R2.dimY(Bd[0])
    assert attendu > 0, attendu
    nu = hoppe_sur_espace(R2, Bd, C, base2, offs2, dims2, degres2, P,
                          np.random.RandomState(5))
    su = hoppe_suffisant_sur_espace(R2, Bd, C, base2, offs2, dims2, degres2,
                                    P, np.random.RandomState(5), D)
    assert nu['stable'] is False, nu
    assert f'= {attendu}' in nu['motif'], (nu['motif'], attendu)
    # (d) contenance
    assert su['stable'] is False, (
        "Hoppe nu elimine mais la forme suffisante non : elle doit le "
        "CONTENIR (H = 0 fait partie du polytope)", su['motif'])

    return (f"#6890 : {n_tw} twists de degre >= 0, tous a source non vide "
            f"(max {max(v['source_max'] for v in suf['detail'].values())}, "
            f"contre {suf['detail'][3]['source_H0']} a H = 0), "
            f"tous a h0 nul -> stable pour J = (1,1,1,1,1) ; controle "
            f"negatif construit : h0(V) = {attendu} = dim H0(O(b_1))")


# ======================================================================
# 11 decies. Phase des twists : elimine #7484, epargne #6890
# ======================================================================

@test("twists : #7484 elimine par H=(-2,0,1), #6890 epargne -- deux verdicts")
def t_hoppe_twists():
    """
    `hoppe_fast` ne testait que H = 0 et H = e_i. Le polytope
    deg_J(H) >= 0 contient des H a composantes de signes MELANGES, hors de
    portee de `max_H` quel qu'il soit, et c'est la que se trouvait un faux
    positif du catalogue.

    ----------------------------------------------------------------------
    Reference : un comptage de dimensions, pas une valeur du code
    ----------------------------------------------------------------------
    dim ker >= dim source - dim cible. Sur #7484 avec H = (-2, 0, 1) :

        h0(O(b_i - H)) = 6 + 0 + 3 + 4 + 0 = 13     (tous certifies)
        h0(O(c   - H)) = 12                         (certifie)

    donc h0(V(-H)) >= 1, et deg_J(H) = 4 > 0 : O(H) est un sous-faisceau de
    pente strictement positive dans un fibre de pente nulle. #7484 n'est
    PAS stable, alors qu'il figure au catalogue `scan_wilson2` comme
    Hoppe-stable. Le test fige les cinq h0 un par un : une derive de la
    cohomologie de Koszul le ferait tomber avant le verdict.

    ----------------------------------------------------------------------
    Deux verdicts opposes
    ----------------------------------------------------------------------
    #6890, demontre stable au §5.14, ne doit PAS etre elimine. Un test qui
    n'exigerait que l'elimination passerait pour une phase qui rejette
    tout ; un test qui n'exigerait que la survie passerait pour une phase
    qui n'elimine jamais. Les deux ensemble n'admettent qu'une phase qui
    discrimine.

    La CERTIFICATION est exigee : sans elle la borne porterait sur des h0
    faux dans ~30 % des cas (§4.2) et l'elimination ne serait pas demontree.
    """
    from cy_landscape.core.hoppe_fast import (hoppe_twists, vecteur_D,
                                              borne_h0_V_twist, hoppe_fast)
    from cy_landscape.core.monads import MonadBundle
    from cy_landscape.core.cache import set_geometry
    from cy_landscape.core.intersection import compute_intersection_numbers
    from cy_landscape.core.exact_cohomology import koszul_cohomology_ex
    from cy_landscape.data.parse_oxford import load_oxford_file

    entries = {e['num']: e for e in load_oxford_file('cicylist.txt')}

    # --- le cas qui doit tomber, avec ses h0 figes un par un -------------
    e = entries[7484]
    amb, cfg = e['ambient'], e['config']
    m = len(amb)
    set_geometry(amb, cfg)
    B = [[0, 1, 1], [1, 0, -1], [0, 0, 1], [1, 0, 1], [1, 0, -1]]
    C = [[3, 1, 1]]
    H = [-2, 0, 1]
    assert [sum(x[k] for x in B) for k in range(m)] == C[0], "c1(V) != 0"
    D = [int(x) for x in vecteur_D(compute_intersection_numbers(amb, cfg),
                                   [1] * m)]
    assert sum(H[k] * D[k] for k in range(m)) > 0, (H, D, "deg_J(H) <= 0")

    attendus = [6, 0, 3, 4, 0]
    src = 0
    for bb, att in zip(B, attendus):
        r = koszul_cohomology_ex(amb, cfg, [bb[k] - H[k] for k in range(m)])
        assert r['certified_by_degree'][0], (bb, "h0 non certifie")
        assert r[0] == att, (bb, r[0], att)
        src += r[0]
    rc = koszul_cohomology_ex(amb, cfg, [C[0][k] - H[k] for k in range(m)])
    assert rc['certified_by_degree'][0] and rc[0] == 12, rc[0]
    assert src == 13 and src - rc[0] == 1, (src, rc[0])

    mo = MonadBundle(B, C)
    borne, cert = borne_h0_V_twist(amb, cfg, mo, H)
    assert cert and borne == 1, (borne, cert)
    tw = hoppe_twists(amb, cfg, mo, D)
    assert tw['instable'] is True, tw
    assert tw['temoin'] == H, (tw['temoin'], H)
    # et hoppe_fast doit relayer le verdict quand on lui donne D
    hf = hoppe_fast(amb, cfg, mo, max_H=1, D=D)
    assert hf['stable'] is False, hf
    # ... alors qu'il le laissait passer sans D : c'est la mesure du gain
    hf0 = hoppe_fast(amb, cfg, mo, max_H=1)
    assert hf0['stable'] is not False, (
        "#7484 est deja elimine sans la phase des twists : le cas de test "
        "ne mesure plus rien", hf0)

    # --- le cas qui doit survivre ---------------------------------------
    e2 = entries[6890]
    amb2, cfg2 = e2['ambient'], e2['config']
    set_geometry(amb2, cfg2)
    B2 = [[0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 1, 0, 0],
          [0, 0, 1, 0, 0], [0, 0, 0, 1, 0]]
    C2 = [[0, 1, 3, 1, 0]]
    D2 = [int(x) for x in vecteur_D(compute_intersection_numbers(amb2, cfg2),
                                    [1] * len(amb2))]
    tw2 = hoppe_twists(amb2, cfg2, MonadBundle(B2, C2), D2)
    assert tw2['instable'] is not True, (
        "#6890 elimine par la phase des twists, alors qu'il est demontre "
        "stable au §5.14 : la phase rejette trop", tw2)
    assert tw2['n_twists'] >= 10, tw2

    # --- la certification doit bloquer la conclusion --------------------
    # Sur #21, 5 twists sur 45 ont un h0 non certifie. Le verdict doit
    # etre None, PAS False : « aucun twist destabilisant trouve » et
    # « aucun twist destabilisant sur des nombres verifies » sont deux
    # enonces differents, et le second seul vaut quelque chose. Sans le
    # suivi de la certification, cette entree passerait pour un blanc-seing.
    e3 = entries[21]
    amb3, cfg3 = e3['ambient'], e3['config']
    set_geometry(amb3, cfg3)
    B3 = [[0, 0, 2, 1, 0], [0, 0, 0, 0, 1], [0, 0, 0, 1, 0],
          [0, 0, 1, 1, 0], [0, 0, 1, 0, 0], [0, 0, 0, 0, 1]]
    C3 = [[0, 0, 4, 3, 2]]
    D3 = [int(x) for x in vecteur_D(compute_intersection_numbers(amb3, cfg3),
                                    [1] * len(amb3))]
    tw3 = hoppe_twists(amb3, cfg3, MonadBundle(B3, C3), D3)
    assert tw3['non_certifies'] > 0, (
        "#21 n'a plus de twist non certifie : le cas de test ne mesure "
        "plus la garde de certification", tw3)
    assert tw3['instable'] is None, (
        "des twists non certifies et pourtant un verdict tranche : la "
        "borne porterait sur des h0 faux dans ~30 % des cas (§4.2)", tw3)

    return (f"#7484 : source {src} > cible {rc[0]} a H = {H} (deg_J = "
            f"{sum(H[k]*D[k] for k in range(m))}), tous h0 certifies -> "
            f"elimine ; #6890 : {tw2['n_twists']} twists, aucun "
            f"destabilisant -> epargne ; #21 : {tw3['non_certifies']}/"
            f"{tw3['n_twists']} non certifies -> aucun verdict")


# ======================================================================
# 11 undecies. Multidegres du certificat de surjectivite
# ======================================================================

@test("multidegres : #21 passe de 4 degres impossibles a 4 viables")
def t_degres_surjectivite():
    """
    Le certificat de surjectivite cherche un multidegre d avec J_d = R_d.
    Le rang etant majore par la source, un d ou dim source < dim cible ne
    peut PAS aboutir : l'echec y est arithmetique et ne dit rien.

    ----------------------------------------------------------------------
    Ce que l'ancienne selection fabriquait
    ----------------------------------------------------------------------
    Elle engendrait des marches LONGUES sur un seul axe et des pas COURTS
    mixtes, jamais les deux a la fois, puis gardait les `n_degres` moins
    chers -- exactement le mauvais sens, la marge s'ameliorant avec la
    taille. Sur #21 (rang 5), les quatre degres retenus ont une marge
    predite de -22, -21, -23 et -22 : le critere ne pouvait aboutir sur
    AUCUN. D'ou 449 couples `indetermine : surjectivite de f non
    certifiee`, dont 420 de rang 5, et le constat du §5.4 « hors de portee
    au rang 5 » -- qui mesurait la liste des directions, pas la geometrie.

    Quatre volets :

    (a) CONTROLE NEGATIF, sur des nombres connus. L'ancienne selection sur
        #21 doit donner des marges toutes STRICTEMENT negatives.

    (b) LA CORRECTION MORD. La nouvelle doit donner des degres viables, et
        la verification EXACTE par `dimY` -- et non par l'estimation dans
        l'ambiant -- doit confirmer source > cible.

    (c) NON-REGRESSION. Les degres qui certifient #6890 et #6947 au §5.4
        doivent rester presents. Sans ce volet, on « corrigerait » #21 en
        cassant les deux seuls candidats du projet.

    (d) LE ZERO FALSY. Ces degres certifiants ont une marge EXACTEMENT
        NULLE (source = cible = 24). Un test ecrit `marge or -1` les rend
        falsy et les ecarte -- c'est le bug trouve en ecrivant ce test. Le
        volet exige donc explicitement qu'une marge nulle soit RETENUE.
    """
    from cy_landscape.core.equivariant_monad import (_degres_a_essayer,
                                                     _marge_predite)
    from cy_landscape.core.sections import Ring
    from cy_landscape.data.parse_oxford import load_oxford_file

    entries = {e['num']: e for e in load_oxford_file('cicylist.txt')}

    # --- (a) et (b) : #21, rang 5 ---------------------------------------
    e21 = entries[21]
    amb21 = e21['ambient']
    B21 = [[0, 0, 2, 1, 0], [0, 0, 0, 0, 1], [0, 0, 0, 1, 0],
           [0, 0, 1, 1, 0], [0, 0, 1, 0, 0], [0, 0, 0, 0, 1]]
    C21 = [0, 0, 4, 3, 2]
    assert [sum(x[k] for x in B21) for k in range(len(amb21))] == C21, "c1(V) != 0"

    anc = _degres_a_essayer(amb21, C21, 6000, 4)
    marges_anc = [_marge_predite(amb21, B21, C21, d) for d in anc]
    assert anc and all(mg is not None and mg < 0 for mg in marges_anc), (
        "l'ancienne selection sur #21 n'est plus impossible : le cas de "
        "test ne mesure plus la correction", list(zip(anc, marges_anc)))

    nou = _degres_a_essayer(amb21, C21, 6000, 4, b_charges=B21)
    marges_nou = [_marge_predite(amb21, B21, C21, d) for d in nou]
    assert nou and all(mg is not None and mg >= 0 for mg in marges_nou), (
        "la nouvelle selection propose encore des degres impossibles",
        list(zip(nou, marges_nou)))

    # verification EXACTE : la marge predite l'est dans l'ambiant, le test
    # reel porte sur R = S/I. Les deux peuvent differer -- on ne se fie pas
    # a l'estimation pour conclure.
    R21 = Ring(amb21, e21['config'], seed=2)
    d0 = nou[0]
    src = sum(R21.dimY([d0[k] - C21[k] + B21[i][k] for k in range(len(amb21))])
              for i in range(len(B21)))
    cible = R21.dimY(d0)
    assert src > cible > 0, (d0, src, cible,
                             "source insuffisante meme sur le degre retenu")

    # --- (c) et (d) : les deux candidats du §2 --------------------------
    ref = {6890: ([[0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 1, 0, 0],
                   [0, 0, 1, 0, 0], [0, 0, 0, 1, 0]],
                  [0, 1, 3, 1, 0], [0, 1, 5, 1, 0]),
           6947: ([[1, 0, 0, 0, 0], [0, 0, 0, 1, 0], [1, 0, 0, 0, 0],
                   [0, 1, 0, 0, 0], [1, 0, 0, 0, 0]],
                  [3, 1, 0, 1, 0], [5, 1, 0, 1, 0])}
    n_zero = 0
    for num, (b, c, dref) in ref.items():
        amb = entries[num]['ambient']
        mg = _marge_predite(amb, b, c, dref)
        # (d) : c'est bien le cas limite, marge nulle
        assert mg == 0, (num, dref, mg, "le degre certifiant du §5.4 n'a "
                                        "plus une marge nulle : le volet "
                                        "du zero falsy ne mesure plus rien")
        n_zero += 1
        lst = _degres_a_essayer(amb, c, 6000, 4, b_charges=b)
        assert dref in lst, (
            f"#{num} : le degre certifiant {dref} du §5.4 a disparu de la "
            f"selection. Une marge EXACTEMENT nulle est falsy en Python : "
            f"un `marge or -1` l'ecarte, et on casse les deux candidats du "
            f"projet en croyant reparer le rang 5", lst)

    return (f"#21 : ancien {marges_anc} (tous impossibles) -> nouveau "
            f"{marges_nou} ; verification exacte source {src} > cible "
            f"{cible} ; les {n_zero} degres certifiants a marge NULLE du "
            f"§5.4 sont conserves")


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
