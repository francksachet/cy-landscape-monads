#!/usr/bin/env python3
"""
equivariance.py -- Test NECESSAIRE d'equivariance des fibres sous Gamma.

--------------------------------------------------------------------------
Le probleme
--------------------------------------------------------------------------
Avoir |chi(V)| = n_gen * |Gamma| sur une CICY a symetrie librement agissante
est necessaire pour descendre sur le quotient, mais pas suffisant. Il faut
encore que l'action de Gamma sur X se releve a V -- que V soit
Gamma-EQUIVARIANT. Sinon le fibre ne descend pas et le modele n'existe pas.

--------------------------------------------------------------------------
La condition testee ici
--------------------------------------------------------------------------
Gamma agit sur l'espace ambiant en melangeant les coordonnees. Sur les
CICYs, cette action permute les facteurs projectifs : un generateur envoie
le facteur P^{n_r} sur un facteur P^{n_s} de meme dimension.

Si sigma est la permutation induite sur les facteurs, alors g* O(b) = O(b o
sigma), c'est-a-dire le fibre en droites dont le vecteur de charges est b
permute. Pour que la somme B = O(b_1) (+) ... (+) O(b_k) soit preservee par
Gamma, il faut donc que l'action de sigma sur les vecteurs permute
l'ENSEMBLE des b_i -- et de meme pour C.

Condition testee : pour chaque generateur, {b o sigma} = {b} comme
multiensembles, et pareil pour les c.

C'est NECESSAIRE, pas suffisant : meme si les ensembles de charges sont
preserves, il reste a verifier que l'application f : B -> C peut etre
choisie equivariante, ce qui est une condition sur les polynomes et non
plus sur les charges. Un fibre qui echoue ici ne descend certainement pas ;
un fibre qui passe reste un candidat.

--------------------------------------------------------------------------
Extraction de la permutation
--------------------------------------------------------------------------
Le fichier de Braun donne, pour chaque symetrie, des generateurs agissant
sur les coordonnees sous forme de matrices. Les coordonnees sont
concatenees facteur par facteur : le facteur r occupe n_r + 1 colonnes
consecutives. Une matrice qui envoie le bloc du facteur r sur celui du
facteur s definit sigma(r) = s. Les matrices diagonales par blocs (phases
seules) donnent l'identite -- l'action est alors interne a chaque facteur
et toute charge est preservee.

Usage:
    python equivariance.py CicyQuotients.m cicylist.txt scan_wilson
"""
import re, sys, json, os, argparse
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from wilson_match import (parse_braun, parse_cicylist, apparier,
                          _ordre_groupe)


# ======================================================================
# Permutation des facteurs induite par une matrice de coordonnees
# ======================================================================

def permutation_facteurs(matrice, dims):
    """
    dims : [n_0, n_1, ...] dimensions des facteurs projectifs.
    Les coordonnees occupent des blocs consecutifs de taille n_r + 1.

    Renvoie sigma sous forme de liste (sigma[r] = s), ou None si la matrice
    ne respecte pas la structure par blocs (auquel cas on ne conclut pas).
    """
    tailles = [n + 1 for n in dims]
    total = sum(tailles)
    if len(matrice) != total or any(len(l) != total for l in matrice):
        return None

    debuts = []
    acc = 0
    for t in tailles:
        debuts.append(acc)
        acc += t

    def bloc_de(i):
        for r, (d, t) in enumerate(zip(debuts, tailles)):
            if d <= i < d + t:
                return r
        return None

    sigma = [None] * len(dims)
    for i, ligne in enumerate(matrice):
        r = bloc_de(i)
        cibles = {bloc_de(j) for j, v in enumerate(ligne) if v != 0}
        if not cibles:
            continue
        if len(cibles) != 1:
            return None                 # melange plusieurs facteurs
        s = cibles.pop()
        if sigma[r] is None:
            sigma[r] = s
        elif sigma[r] != s:
            return None
    if any(s is None for s in sigma):
        return None
    if sorted(sigma) != list(range(len(dims))):
        return None                     # pas une permutation
    return sigma


def _premier_bloc(txt):
    """Premier groupe d'accolades equilibrees de txt, ou None."""
    debut = txt.find('{')
    if debut < 0:
        return None
    prof = 0
    for i in range(debut, len(txt)):
        if txt[i] == '{':
            prof += 1
        elif txt[i] == '}':
            prof -= 1
            if prof == 0:
                return txt[debut:i + 1]
    return None


def _generateurs_coordonnees(suite):
    """
    Liste des matrices de generateurs agissant sur les coordonnees.

    Format Braun, apres le nom du groupe :
        {{matrice_1}, {matrice_2}, ...}, {{poly_1}, ...}, H11 -> ...
    Le premier bloc equilibre est donc la LISTE des generateurs. Un groupe
    comme Z3 x Z3 en a deux ; la version precedente concatenait leurs lignes
    en une seule matrice, dont la taille ne correspondait plus au nombre de
    coordonnees -- d'ou son rejet. C'est ce qui faisait disparaitre tous les
    groupes a plusieurs generateurs, Z3 x Z3 en tete.

    On decoupe donc la liste en ses elements de profondeur 1.
    """
    bloc = _premier_bloc(suite)
    if bloc is None:
        return []
    interieur = bloc[1:-1]
    mats = []
    prof = 0
    debut = None
    for i, ch in enumerate(interieur):
        if ch == '{':
            if prof == 0:
                debut = i
            prof += 1
        elif ch == '}':
            prof -= 1
            if prof == 0 and debut is not None:
                mats.append(interieur[debut:i + 1])
                debut = None
    return mats


def _lire_matrice(txt):
    """
    Matrice de generateur, sous forme d'un motif ZERO / NON ZERO.

    Les generateurs d'ordre 3, 5, ... ne peuvent pas avoir d'entrees
    entieres : Braun les ecrit avec des racines de l'unite, sous des formes
    variees (Exp[2 I Pi/3], w, Zeta, ...). L'ancienne lecture ne
    retenait que les entiers, ce qui coupait ces lignes en morceaux de la
    mauvaise taille et faisait rejeter le generateur -- d'ou les 1006 puis
    151 « generateurs non exploitables », concentres precisement sur les
    groupes d'ordre impair.

    Pour extraire la PERMUTATION DES FACTEURS, seule la position des
    entrees non nulles compte ; leur valeur est sans importance. On decoupe
    donc chaque ligne sur les virgules de premier niveau et on marque
    chaque entree 0 ou 1.
    """
    lignes = re.findall(r'\{([^{}]*)\}', txt)
    m = []
    for l in lignes:
        entrees = _decoupe_virgules(l)
        if not entrees:
            continue
        m.append([0 if _est_nul(e) else 1 for e in entrees])
    return m


def _decoupe_virgules(txt):
    """Decoupe sur les virgules hors crochets (Exp[2 I Pi/3] reste entier)."""
    out, prof, cur = [], 0, []
    for ch in txt:
        if ch in '[(':
            prof += 1
        elif ch in '])':
            prof -= 1
        if ch == ',' and prof == 0:
            out.append(''.join(cur)); cur = []
        else:
            cur.append(ch)
    if cur:
        out.append(''.join(cur))
    return [x.strip() for x in out if x.strip()]


def _est_nul(entree):
    e = entree.strip()
    if not e:
        return True
    try:
        return int(e) == 0
    except ValueError:
        return False        # expression symbolique -> non nulle


# ======================================================================
# Test sur un fibre
# ======================================================================

def charges_preservees(charges, sigma):
    """{b o sigma} == {b} comme multiensembles."""
    ref = Counter(tuple(b) for b in charges)
    img = Counter(tuple(b[sigma[r]] for r in range(len(sigma))) for b in charges)
    return ref == img


def teste_fibre(b_charges, c_charges, sigmas):
    """
    Renvoie la liste des SYMETRIES sous lesquelles B et C sont preserves.

    Une symetrie peut avoir plusieurs generateurs (Z3 x Z3 en a deux). Pour
    que le groupe agisse sur le fibre, il faut que TOUS ses generateurs
    preservent les ensembles de charges -- pas seulement l'un d'eux. La
    version precedente acceptait des qu'un generateur convenait, ce qui
    laissait passer des fibres brises par le second : c'est pour cela
    qu'elle n'eliminait presque personne.

    Les generateurs d'une meme symetrie sont reperes par le nom de base,
    le script les ayant nommes "Nom", "Nom#2", "Nom#3"...
    """
    par_symetrie = {}
    for cle, sigma in sigmas:
        par_symetrie.setdefault(cle, []).append(sigma)

    ok = []
    for cle, gens in par_symetrie.items():
        if all(charges_preservees(b_charges, s) and
               charges_preservees(c_charges, s) for s in gens):
            ok.append(cle[1])
    return ok


def _sortie_tolerante():
    """
    Empeche un plantage d'encodage sur une console Windows.

    Les etiquettes de groupe de jauge contiennent des indices Unicode
    ("E\u2086"), que la console cp1252 ne sait pas encoder. `print` levait
    alors UnicodeEncodeError, le script mourait sans rien ecrire, et la
    commande suivante de l'enchainement tournait sur un fichier absent.
    Cas reel : scan_exh2, ou `equivariance.py` est mort ainsi et a fait
    echouer `equivariance_f.py` derriere lui.
    """
    import sys as _sys
    for flux in (_sys.stdout, _sys.stderr):
        try:
            flux.reconfigure(errors='replace')
        except (AttributeError, ValueError):
            pass


def main():
    _sortie_tolerante()
    ap = argparse.ArgumentParser()
    ap.add_argument('braun_m')
    ap.add_argument('cicylist')
    ap.add_argument('output_dir', help="dossier du scan (results.jsonl)")
    ap.add_argument('--input', default='results.jsonl')
    ap.add_argument('--n-gen', type=int, default=3,
                    help="Nombre de generations voulu sur le quotient.")
    ap.add_argument('--top', type=int, default=30)
    args = ap.parse_args()

    braun = parse_braun(args.braun_m)
    cicy = parse_cicylist(args.cicylist)
    corresp, _amb, _abs = apparier(braun, cicy)
    inv = {v: k for k, v in corresp.items()}
    par_num = {b['num']: b for b in braun}
    dims_par_cicy = {e['num']: [sum(r) - 1 for r in e['conf']] for e in cicy}

    # permutations disponibles par CICY (numerotation utilisateur)
    txt = open(args.braun_m, encoding='utf-8', errors='replace').read()
    sigmas_par_cicy = {}
    non_lus = 0
    for num_u, num_b in inv.items():
        dims = dims_par_cicy.get(num_u)
        if not dims:
            continue
        i = txt.find(f'Num -> {num_b},')
        if i < 0:
            continue
        j = txt.find('Num ->', i + 6)
        bloc = txt[i:j if j > 0 else len(txt)]
        k = bloc.find('Symmetries')
        if k < 0:
            continue
        # Une entree de symetrie commence par son nom entre guillemets ;
        # les matrices de generateurs sur les coordonnees suivent
        # immediatement. On lit la premiere accolade equilibree apres le nom.
        sig = []
        for idx_sym, m in enumerate(re.finditer(r'"([^"]+)"', bloc)):
            suite = bloc[m.end():]
            for k, brut_txt in enumerate(_generateurs_coordonnees(suite)):
                s = permutation_facteurs(_lire_matrice(brut_txt), dims)
                if s is None:
                    non_lus += 1
                    continue
                # Cle = (rang de l'entree de symetrie, nom). Surtout pas une
                # concatenation avec '#' : les noms de Braun en contiennent
                # deja ("Z2 x Z2#2", "Z4 x Z2$") pour distinguer des actions
                # inequivalentes du meme groupe abstrait. Les separer par '#'
                # fusionnait deux symetries distinctes en une seule, et le
                # test exigeait alors les generateurs des deux a la fois.
                sig.append(((idx_sym, m.group(1)), s))
        if sig:
            sigmas_par_cicy[num_u] = sig

    src = os.path.join(args.output_dir, args.input)
    rs = [json.loads(l) for l in open(src, encoding='utf-8') if l.strip()]

    print(f"\n{'='*74}")
    print(f"  TEST NECESSAIRE D'EQUIVARIANCE -- {len(rs)} candidats")
    print(f"{'='*74}")
    print(f"  CICYs avec permutation extraite : {len(sigmas_par_cicy)}")
    if non_lus:
        print(f"  Generateurs non exploitables    : {non_lus} "
              f"(structure par blocs non reconnue)")

    retenus, rejetes, inconnus = [], [], []
    for r in rs:
        sig = sigmas_par_cicy.get(r['cicy'])
        if not sig:
            inconnus.append(r); continue
        # Une permutation identite ne teste rien : toute charge la verifie.
        # On distingue donc les fibres qui passent un test NON TRIVIAL de ceux
        # dont la CICY n'offre que des permutations triviales (action interne
        # a chaque facteur, par phases). Sans cette distinction, un « passe »
        # peut n'etre qu'une absence de contrainte.
        # Une symetrie n'est « non triviale » que si au moins un de ses
        # generateurs permute reellement des facteurs. On garde alors TOUS
        # ses generateurs, sans quoi on testerait un sous-groupe.
        cles_nt = {cle for cle, sg in sig if sg != list(range(len(sg)))}
        non_triv = [(cle, sg) for cle, sg in sig if cle in cles_nt]
        ok = teste_fibre(r['b_charges'], r['c_charges'], sig)
        ok_nt = teste_fibre(r['b_charges'], r['c_charges'], non_triv)
        if ok:
            # Un groupe n'est utilisable que si son ORDRE ramene l'indice
            # amont a n_gen : |Gamma| = n_gen_amont / n_gen. Un fibre
            # d'indice 60 compatible avec un Z2 donnerait 30 generations,
            # pas 3 -- la compatibilite des charges ne suffit pas.
            amont = r.get('n_gen_amont') or abs(r['cohomology'][1] - r['cohomology'][2])
            besoin = amont / args.n_gen if args.n_gen else None
            utiles = [g for g in sorted(set(ok))
                      if besoin and _ordre_groupe(g) == besoin]
            r2 = dict(r)
            r2['equivariant_possible'] = sorted(set(ok))
            r2['groupes_utiles'] = utiles
            r2['test_non_trivial'] = bool(non_triv)
            r2['equivariant_non_trivial'] = sorted(set(ok_nt))
            retenus.append(r2)
        else:
            rejetes.append(r)

    n_nt = sum(1 for r in retenus if r.get('test_non_trivial'))
    print(f"\n  Passent le test necessaire : {len(retenus)}")
    print(f"    dont via une permutation NON TRIVIALE : {n_nt}")
    print(f"    dont sur une CICY sans permutation de facteurs : "
          f"{len(retenus) - n_nt}  (aucune contrainte : non concluant)")
    print(f"  Elimines                   : {len(rejetes)}")
    print(f"  Indetermines               : {len(inconnus)}")

    if retenus:
        propres = [r for r in retenus
                   if min(r['cohomology'][1], r['cohomology'][2]) == 0]
        avec_groupe = [r for r in retenus if r.get('groupes_utiles')]
        print(f"\n  Dont un groupe d'ORDRE compatible avec l'indice : "
              f"{len(avec_groupe)}")
        print(f"  Dont n_anti = 0 : {len(propres)}")
        cibles = [r for r in avec_groupe
                  if min(r['cohomology'][1], r['cohomology'][2]) == 0]
        print(f"  Les deux a la fois : {len(cibles)}   <<< candidats reels")
        print(f"\n    {'#CICY':>6} {'jauge':>7} {'rk':>2} {'n_anti':>6} "
              f"{'cohomologie':>16}  groupes compatibles")
        pool = sorted(cibles or propres or retenus,
                      key=lambda x: (min(x['cohomology'][1], x['cohomology'][2]),
                                     -max(x.get('ordres_gamma') or [0])))
        for r in pool[:args.top]:
            h = r['cohomology']
            print(f"    {r['cicy']:>6} {r.get('gauge',''):>7} "
                  f"{r.get('rank_V',''):>2} {min(h[1],h[2]):>6} "
                  f"{str(h):>16}  "
                  f"{', '.join(r.get('groupes_utiles') or r['equivariant_possible'])[:36]}")

    out = os.path.join(args.output_dir, 'results_equivariant.jsonl')
    with open(out, 'w', encoding='utf-8') as f:
        for r in retenus:
            f.write(json.dumps(r) + '\n')
    print(f"\n  Ecrit : {out}")
    print(f"{'='*74}\n")


if __name__ == '__main__':
    sys.exit(main())
