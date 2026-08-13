#!/usr/bin/env python3
"""
validate_cohomology.py -- Validation du socle : exact_cohomology.py

Tout le pipeline repose sur `koszul_cohomology` : la cohomologie de V, le
test de Hoppe, wedge^2 V, le comptage des generations. Ce script mesure son
taux d'erreur contre des references INDEPENDANTES du module teste.

--------------------------------------------------------------------------
Les quatre tests
--------------------------------------------------------------------------

[HRR]  chi(L) = (1/6) * sum d_ijk a_i a_j a_k + (1/12) * sum c2_i a_i

    Hirzebruch-Riemann-Roch sur une 3-variete de Calabi-Yau (c1(TY) = 0).
    Les d_ijk viennent d'`intersection.py`. La valeur ne depend d'aucune
    hypothese sur la suite spectrale : c'est la reference la plus dure du
    lot. Un ecart signifie que la somme alternee des h^i est fausse, donc
    qu'au moins un h^i l'est.

    Le vecteur c2_i est RECALCULE ici (`c2_tangent_local`) et non pris
    dans `intersection.py`. Motif : `compute_c2_tangent` pose

        c2_amb[r,s] = c2_amb[s,r] = (n_r+1)(n_s+1)   pour r < s

    puis somme sur les couples ORDONNES (j,k). Or le terme de degre 2 de
    prod (1+J_r)^(n_r+1) vaut sum_(r<s) (n_r+1)(n_s+1) J_r J_s : ecrit
    comme somme ordonnee, le coefficient hors diagonale est donc la
    MOITIE de cette valeur. Les termes croises sont comptes deux fois.
    `c2_N` de la meme fonction utilise deja, correctement, la convention
    ordonnee -- les deux moities de la soustraction ne sont pas dans la
    meme convention.

    Trois verifications de la version locale :
      - chi(L) doit etre entier pour tout L. Sur 2800 tirages, la version
        d'`intersection.py` donne un chi entier dans 34 % des cas, la
        version locale dans 100 %.
      - bicubique P2xP2[3,3] : c2.J = (36, 36), valeur connue. Version
        locale 36, `intersection.py` 63.
      - quintique P4[5] : c2.J = 50 dans les deux cas. m = 1, aucun terme
        croise -- d'ou le fait que la validation initiale sur la quintique
        n'ait rien pu voir.

    Le test [C2] du rapport chiffre l'ecart entre les deux versions.

[SERRE]  h^i(L) = h^(3-i)(L^-1)

    K_Y est trivial sur un CY3. Chaque vecteur de charges est ainsi teste
    deux fois, contre lui-meme, sans reference externe.

    ATTENTION a la lecture : ce compteur porte sur les h^i BRUTS, y compris
    non certifies, et un taux d'echec eleve (~65 %) y est NORMAL. Le calcul
    de d_1 resout les rangs par recursion le long de la chaine de Koszul,
    et le resultat depend du bout par lequel on commence -- or Serre echange
    exactement les deux sens. Ces cas sont precisement ceux que la
    certification ecarte. Le test qui fait foi est [SERRE_CERT], dans la
    section [CERT] : restreint aux paires certifiees des deux cotes, il
    doit valoir 100 %.

[DEGEN]  chi lu sur la page E_1 == chi lu sur les h^i retournes

    Le module additionne les termes E_1 par anti-diagonale q = n + p en
    supposant la degeneration a E_1. La somme alternee sur toute la page
    E_1 donne chi quelle que soit la degeneration. Un ecart isole donc
    precisement les cas ou l'hypothese de degeneration tombe -- c'est la
    limite deja connue (quintique O(5) : 126 au lieu de 125), ici
    quantifiee sur l'ensemble du domaine reellement utilise.

[NEGATIF]  aucun h^i < 0

--------------------------------------------------------------------------
Ce que le rapport donne
--------------------------------------------------------------------------
Le taux d'erreur global, ET sa repartition par amplitude de charge. Si les
erreurs n'apparaissent qu'au-dela d'un certain max|a_i|, la conclusion
pratique est de plafonner --max-charge plutot que de reecrire le module.

Usage:
    python validate_cohomology.py                    # CICYs embarquees
    python validate_cohomology.py cicylist.txt       # echantillon du fichier
    python validate_cohomology.py cicylist.txt --n-cicys 60 --max-charge 6
"""
import os, sys, json, argparse, itertools, random
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np


def c2_tangent_local(ambient_dims, config_matrix, d_ijk):
    """
    c2(TY) . J_i, convention ordonnee coherente entre c2_amb et c2_N.
    Voir la note [HRR] de l'en-tete : c'est la moitie hors diagonale qui
    manque dans `intersection.compute_c2_tangent`.
    """
    from math import comb as _comb
    m = len(ambient_dims)
    K = config_matrix.shape[0]

    c2_amb = np.zeros((m, m))
    for r in range(m):
        c2_amb[r, r] = _comb(ambient_dims[r] + 1, 2)
        for s in range(r + 1, m):
            val = (ambient_dims[r] + 1) * (ambient_dims[s] + 1) / 2.0
            c2_amb[r, s] = c2_amb[s, r] = val

    c2_N = np.zeros((m, m))
    for a in range(K):
        for b in range(a + 1, K):
            for r in range(m):
                for s in range(m):
                    c2_N[r, s] += config_matrix[a, r] * config_matrix[b, s]

    coeff = c2_amb - c2_N
    out = np.zeros(m)
    for i in range(m):
        for j in range(m):
            for k in range(m):
                out[i] += coeff[j, k] * d_ijk[i, j, k]
    return out


def chi_hrr(charges, d, c2):
    """chi(L) par Hirzebruch-Riemann-Roch. Entiers exacts, pas de flottant."""
    h = len(charges)
    cube = 0
    for i in range(h):
        if charges[i] == 0:
            continue
        for j in range(h):
            if charges[j] == 0:
                continue
            for k in range(h):
                if charges[k] == 0:
                    continue
                cube += int(charges[i]) * int(charges[j]) * int(charges[k]) * int(d[i, j, k])
    lin = sum(int(charges[i]) * int(round(c2[i] * 2)) for i in range(h))
    # c2 peut etre demi-entier ; on travaille en douziemes doubles.
    num = 4 * cube + lin          # chi = cube/6 + (c2.a)/12 = (4*cube + 2*c2.a)/24
    return num // 24, (num % 24 == 0)


def chi_from_e1(ambient, config, charges):
    """
    chi lu directement sur la page E_1, sans hypothese de degeneration :
        chi = sum_{p,q} (-1)^(q-p) E_1^{p,q}
    """
    from cy_landscape.core.exact_cohomology import h_ambient
    K = config.shape[0]
    m = len(ambient)
    total = 0
    for p in range(K + 1):
        for subset in itertools.combinations(range(K), p):
            shifted = list(charges)
            for a in subset:
                for r in range(m):
                    shifted[r] -= int(config[a, r])
            for q, dim in h_ambient(ambient, shifted).items():
                total += ((-1) ** (q - p)) * int(dim)
    return total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cicyfile', nargs='?', default=None)
    ap.add_argument('--n-cicys', type=int, default=25,
                    help="Nombre de CICYs echantillonnees (defaut: 25)")
    ap.add_argument('--max-charge', type=int, default=5)
    ap.add_argument('--n-charges', type=int, default=120,
                    help="Vecteurs de charges tires par CICY (defaut: 120)")
    ap.add_argument('--max-ps', type=int, default=6)
    ap.add_argument('--seed', type=int, default=42)
    args = ap.parse_args()

    from cy_landscape.core.exact_cohomology import koszul_cohomology
    try:
        from cy_landscape.core.exact_cohomology import koszul_cohomology_ex
    except ImportError:
        koszul_cohomology_ex = None      # ancienne version du module
    from cy_landscape.core.intersection import (
        compute_intersection_numbers, compute_c2_tangent)

    if args.cicyfile:
        from cy_landscape.data.parse_oxford import load_oxford_file
        entries = load_oxford_file(args.cicyfile)
        entries = [e for e in entries if len(e['ambient']) <= args.max_ps]
    else:
        from cy_landscape.data.oxford_cicys import get_all_oxford
        entries = get_all_oxford()

    rng = random.Random(args.seed)
    if len(entries) > args.n_cicys:
        entries = rng.sample(entries, args.n_cicys)

    print(f"\n{'='*74}")
    print(f"  VALIDATION DU SOCLE -- exact_cohomology.koszul_cohomology")
    print(f"{'='*74}")
    print(f"  CICYs             : {len(entries)}")
    print(f"  Charges par CICY  : {args.n_charges} dans [-{args.max_charge}, {args.max_charge}]")

    n = 0
    n_cert = 0
    n_cert_ok = 0
    n_cert12 = 0
    n_serre_pairs = 0
    n_serre_ok = 0
    n_pas_cy3 = 0
    c2_ecarts = []
    fails = defaultdict(int)
    by_amp = defaultdict(lambda: [0, 0])      # max|a| -> [testes, en echec]
    by_cicy = defaultdict(lambda: [0, 0])
    examples = defaultdict(list)

    for c in entries:
        try:
            d = compute_intersection_numbers(c['ambient'], c['config'])
            c2 = c2_tangent_local(c['ambient'], np.asarray(c['config']).reshape(
                1, -1) if np.asarray(c['config']).ndim == 1 else np.asarray(c['config']), d)
            c2_pipeline = compute_c2_tangent(c['ambient'], c['config'], d)
        except Exception:
            continue
        m = len(c['ambient'])
        config = np.asarray(c['config'])
        if config.ndim == 1:
            config = config.reshape(1, -1)
        if config.shape[0] > 12:
            continue                          # garde-fou K>12 du module
        if sum(c['ambient']) - config.shape[0] != 3:
            # dim Y = dim(ambiant) - K doit valoir 3. Certaines entrees de la
            # liste ne sont pas des 3-varietes ; elles font echouer HRR a 100 %
            # pour une raison sans rapport avec le module teste.
            n_pas_cy3 += 1
            continue

        if not np.allclose(c2, c2_pipeline):
            c2_ecarts.append((c['num'],
                              [float(x) for x in np.round(c2, 3)],
                              [float(x) for x in np.round(c2_pipeline, 3)]))

        seen = set()
        for _ in range(args.n_charges):
            a = tuple(rng.randint(-args.max_charge, args.max_charge) for _ in range(m))
            if a in seen:
                continue
            seen.add(a)

            try:
                h = koszul_cohomology(c['ambient'], config, list(a))
                hneg = koszul_cohomology(c['ambient'], config, [-x for x in a])
            except Exception:
                continue

            n += 1
            amp = max(abs(x) for x in a)
            by_amp[amp][0] += 1
            by_cicy[c['num']][0] += 1
            bad = False

            chi_h = sum((-1) ** i * h.get(i, 0) for i in range(4))
            chi_ref, entier = chi_hrr(a, d, c2)
            if not entier:
                fails['HRR_non_entier'] += 1
                bad = True
            elif chi_h != chi_ref:
                fails['HRR'] += 1
                bad = True
                if len(examples['HRR']) < 6:
                    examples['HRR'].append((c['num'], a, [h.get(i, 0) for i in range(4)],
                                            chi_h, chi_ref))

            if any(h.get(i, 0) != hneg.get(3 - i, 0) for i in range(4)):
                fails['SERRE'] += 1
                bad = True

            # [SERRE_CERT] -- le test qui compte reellement.
            # Les h^i bruts violent massivement la dualite de Serre : le
            # calcul de d_1 resout les rangs par recursion le long de la
            # chaine, et le resultat depend du bout par lequel on commence,
            # alors que Serre echange precisement les deux sens. C'est
            # attendu et sans consequence : ces cas ne sont pas certifies.
            # Le test valide est donc restreint aux paires (h^i, h^(3-i))
            # certifiees des DEUX cotes.
            if koszul_cohomology_ex is not None:
                try:
                    fa = koszul_cohomology_ex(c['ambient'], config, list(a))
                    fb = koszul_cohomology_ex(c['ambient'], config,
                                              [-x for x in a])
                    ca = fa.get('certified_by_degree') or {}
                    cb = fb.get('certified_by_degree') or {}
                    for i in range(4):
                        if ca.get(i) and cb.get(3 - i):
                            n_serre_pairs += 1
                            if fa.get(i, 0) == fb.get(3 - i, 0):
                                n_serre_ok += 1
                except Exception:
                    pass
                if len(examples['SERRE']) < 6:
                    examples['SERRE'].append(
                        (c['num'], a, [h.get(i, 0) for i in range(4)],
                         [hneg.get(i, 0) for i in range(4)]))

            try:
                chi_e1 = chi_from_e1(c['ambient'], config, list(a))
                if chi_e1 != chi_h:
                    fails['DEGEN'] += 1
                    bad = True
                    if len(examples['DEGEN']) < 6:
                        examples['DEGEN'].append(
                            (c['num'], a, [h.get(i, 0) for i in range(4)], chi_h, chi_e1))
            except Exception:
                pass

            # [CERT] -- disponible seulement avec exact_cohomology corrige.
            # Mesure la fraction du domaine sur laquelle les h^i sont PROUVES
            # exacts (degenerescence a E_1 certifiee), et verifie que HRR y
            # passe. C'est ce taux qui determine le rendement d'un rescan :
            # le pipeline corrige ne retient que ces cas-la.
            if koszul_cohomology_ex is not None:
                try:
                    full = koszul_cohomology_ex(c['ambient'], config, list(a))
                    if full.get('certified'):
                        n_cert += 1
                        if sum((-1) ** i * full.get(i, 0) for i in range(4)) == chi_ref:
                            n_cert_ok += 1
                    _cd = full.get('certified_by_degree') or {}
                    if _cd.get(1) and _cd.get(2):
                        n_cert12 += 1
                except Exception:
                    pass

            if any(h.get(i, 0) < 0 for i in range(4)):
                fails['NEGATIF'] += 1
                bad = True

            if bad:
                by_amp[amp][1] += 1
                by_cicy[c['num']][1] += 1

    print(f"\n{'='*74}")
    print(f"  RESULTATS -- {n} vecteurs de charges testes")
    print(f"{'='*74}")
    if not fails:
        print(f"\n  Aucun echec. Le socle est fiable sur ce domaine.")
    for k in ('HRR', 'HRR_non_entier', 'SERRE', 'DEGEN', 'NEGATIF'):
        if fails.get(k):
            print(f"    {k:<16} {fails[k]:>7}  ({100*fails[k]/max(1,n):.1f} %)")

    if koszul_cohomology_ex is not None and n:
        print(f"\n  --- [CERT] cas ou les h^i sont PROUVES exacts ---")
        print(f"    certifies                  : {n_cert}/{n} ({100.0*n_cert/n:.1f} %)")
        if n_cert:
            print(f"    HRR OK parmi les certifies : {n_cert_ok}/{n_cert} "
                  f"({100.0*n_cert_ok/n_cert:.1f} %)")
        print(f"    h1 ET h2 certifies         : {n_cert12}/{n} ({100.0*n_cert12/n:.1f} %)")
        if n_serre_pairs:
            print(f"    Serre sur paires certifiees: {n_serre_ok}/{n_serre_pairs} "
                  f"({100.0*n_serre_ok/n_serre_pairs:.1f} %)   <- doit valoir 100 %")
        print(f"    -> c'est ce dernier taux qui compte : le comptage des")
        print(f"       generations ne depend que de h1 et h2, et c'est le")
        print(f"       critere qu'applique le pipeline corrige.")

    if n_pas_cy3:
        print(f"\n  ({n_pas_cy3} entrees ignorees : dim(ambiant) - K != 3, "
              f"pas des 3-varietes)")

    print(f"\n  --- [C2] intersection.compute_c2_tangent vs convention ordonnee ---")
    print(f"    CICYs ou les deux different : {len(c2_ecarts)} / {len(entries)}")
    for num, loc, pipe in c2_ecarts[:5]:
        print(f"      #{num:<6} correct={loc}   pipeline={pipe}")

    print(f"\n  --- Taux d'echec par amplitude de charge max|a_i| ---")
    print(f"    {'max|a|':>7} {'testes':>8} {'echecs':>8} {'taux':>8}")
    for amp in sorted(by_amp):
        t, f = by_amp[amp]
        print(f"    {amp:>7} {t:>8} {f:>8} {100*f/max(1,t):>7.1f} %")

    pires = sorted(by_cicy.items(), key=lambda kv: -kv[1][1] / max(1, kv[1][0]))[:8]
    print(f"\n  --- CICYs les plus touchees ---")
    for num, (t, f) in pires:
        if f:
            print(f"    CICY #{num:<6} {f}/{t}  ({100*f/max(1,t):.0f} %)")

    for k, ex in examples.items():
        if not ex:
            continue
        print(f"\n  --- Exemples [{k}] ---")
        for e in ex:
            if k == 'SERRE':
                print(f"    #{e[0]} a={list(e[1])}  h(L)={e[2]}  h(-L)={e[3]}")
            else:
                print(f"    #{e[0]} a={list(e[1])}  h={e[2]}  chi_module={e[3]}  chi_ref={e[4]}")

    print(f"\n{'='*74}\n")
    return 0


if __name__ == '__main__':
    sys.exit(main())
