"""
braun_symmetry.py -- Lecture COMPLETE des symetries de CicyQuotients.m.

--------------------------------------------------------------------------
Pourquoi ce module
--------------------------------------------------------------------------
Jusqu'ici deux lectures partielles coexistaient :

  - `wilson_match.parse_braun` ne retenait que les NOMS des groupes ;
  - `equivariance._generateurs_coordonnees` ne lisait que le PREMIER bloc de
    generateurs (l'action sur les coordonnees), et seulement son motif
    zero / non-zero, ce qui suffit pour extraire la permutation des facteurs
    mais pas pour calculer quoi que ce soit.

Le SECOND bloc -- l'action de Gamma sur les polynomes definissants -- n'etait
lu nulle part. C'est pourtant lui qui dit comment choisir les coefficients
des polynomes pour que l'ideal soit Gamma-invariant. Sans lui,
`sections.Ring` tire des coefficients ALEATOIRES, l'ideal n'est pas preserve,
l'action ne descend pas au quotient, et tout test d'equivariance construit
par-dessus porte sur une geometrie qui n'est pas celle annoncee.

Ce module lit la structure entiere, sans l'interpreter :

    {{gap_id_1, gap_id_2}, "Nom", {gen_coord, ...}, {gen_poly, ...},
     H11 -> a, H21 -> b}

--------------------------------------------------------------------------
Analyse syntaxique
--------------------------------------------------------------------------
On n'utilise PAS d'expression reguliere pour la structure : les entrees
contiennent des `rt[4]`, `-rt[3]^2`, `Exp[2 I Pi/3]`, et les noms de groupes
contiennent eux-memes `#` et `$`. On tokenise donc les accolades en
respectant l'imbrication, et les crochets et parentheses sont traites comme
opaques -- une virgule a l'interieur de `Exp[2 I Pi/3]` ne separe rien.

Les feuilles sont rendues telles quelles, en chaines. Leur conversion en
elements de GF(p) est le travail de `gamma_action.parse_entree`.
"""
import re


# ======================================================================
# Tokenisation Mathematica
# ======================================================================

def split_top(s):
    """
    Decoupe `s` sur les virgules de PROFONDEUR ZERO.

    Les trois familles de delimiteurs comptent dans la profondeur : `{}` pour
    les listes, `[]` pour les appels de fonction (`rt[4]`, `Exp[...]`), `()`
    pour les groupements. Une virgule interne a l'un d'eux ne separe rien.
    """
    out, cur, prof = [], [], 0
    for ch in s:
        if ch in '{[(':
            prof += 1
        elif ch in '}])':
            prof -= 1
        if ch == ',' and prof == 0:
            out.append(''.join(cur))
            cur = []
        else:
            cur.append(ch)
    if cur:
        out.append(''.join(cur))
    return [x.strip() for x in out if x.strip()]


def parse_liste(s):
    """
    Liste Mathematica -> listes Python imbriquees, feuilles en chaines.

    `{{1, 0}, {0, rt[4]}}` -> [['1', '0'], ['0', 'rt[4]']]
    """
    s = s.strip()
    if s.startswith('{') and s.endswith('}'):
        return [parse_liste(x) for x in split_top(s[1:-1])]
    return s


def bloc_equilibre(txt, depart=0, ouvrant='{', fermant='}'):
    """
    Premier bloc d'accolades equilibrees a partir de `depart`.

    Renvoie (sous_chaine, index_juste_apres) ou (None, -1).
    """
    i = txt.find(ouvrant, depart)
    if i < 0:
        return None, -1
    prof = 0
    for j in range(i, len(txt)):
        if txt[j] == ouvrant:
            prof += 1
        elif txt[j] == fermant:
            prof -= 1
            if prof == 0:
                return txt[i:j + 1], j + 1
    return None, -1


# ======================================================================
# Lecture des entrees
# ======================================================================

def _blocs_num(txt):
    """Decoupe le fichier en blocs `Num -> n, ...`."""
    i = txt.find('CicyQuotients=')
    if i >= 0:
        txt = txt[i + len('CicyQuotients='):]
    positions = [(int(m.group(1)), m.start())
                 for m in re.finditer(r'Num\s*->\s*(\d+)', txt)]
    for k, (num, deb) in enumerate(positions):
        fin = positions[k + 1][1] if k + 1 < len(positions) else len(txt)
        yield num, txt[deb:fin]


def _entier(bloc, motif):
    m = re.search(motif, bloc)
    return int(m.group(1)) if m else None


def parse_symmetries(path):
    """
    {num_braun: {'conf': matrice, 'h11':, 'h21':, 'symetries': [...]}}

    Chaque symetrie :
        {'gap': (ordre, indice), 'nom': str,
         'coord': [matrice, ...],   matrices sur les COORDONNEES (chaines)
         'poly':  [matrice, ...],   matrices sur les POLYNOMES  (chaines)
         'h11':, 'h21':}

    `coord` et `poly` ont la meme longueur : un generateur du groupe donne
    une matrice de chaque cote. Si ce n'est pas le cas, l'entree est marquee
    `'suspect': True` plutot que rejetee silencieusement.
    """
    txt = open(path, encoding='utf-8', errors='replace').read()
    out = {}
    for num, bloc in _blocs_num(txt):
        conf_txt, _ = bloc_equilibre(bloc, bloc.find('Conf ->'))
        conf = None
        if conf_txt:
            brut = parse_liste(conf_txt)
            try:
                conf = [[int(x) for x in ligne] for ligne in brut]
            except (TypeError, ValueError):
                conf = None

        k = bloc.find('Symmetries')
        symetries = []
        if k >= 0:
            sym_txt, _ = bloc_equilibre(bloc, k)
            if sym_txt:
                for entree in split_top(sym_txt[1:-1]):
                    s = _parse_une_symetrie(entree)
                    if s:
                        symetries.append(s)

        out[num] = {
            'conf': conf,
            'h11': _entier(bloc, r'H11\s*->\s*(-?\d+)'),
            'h21': _entier(bloc, r'H21\s*->\s*(-?\d+)'),
            'symetries': symetries,
        }
    return out


def _parse_une_symetrie(entree):
    """
    `{{4, 2}, "Z2 x Z2", {coord...}, {poly...}, H11 -> 3, H21 -> 11}`
    """
    entree = entree.strip()
    if not (entree.startswith('{') and entree.endswith('}')):
        return None
    champs = split_top(entree[1:-1])
    if len(champs) < 4:
        return None

    gap = parse_liste(champs[0])
    try:
        gap = tuple(int(x) for x in gap)
    except (TypeError, ValueError):
        gap = None

    nom = champs[1].strip()
    if nom.startswith('"') and nom.endswith('"'):
        nom = nom[1:-1]

    coord = parse_liste(champs[2])
    poly = parse_liste(champs[3])

    h11 = h21 = None
    for c in champs[4:]:
        m = re.search(r'H11\s*->\s*(-?\d+)', c)
        if m:
            h11 = int(m.group(1))
        m = re.search(r'H21\s*->\s*(-?\d+)', c)
        if m:
            h21 = int(m.group(1))

    return {
        'gap': gap, 'nom': nom,
        'coord': coord, 'poly': poly,
        'h11': h11, 'h21': h21,
        'suspect': len(coord) != len(poly),
    }


# ======================================================================
# Conversion en matrices numeriques sur GF(p)
# ======================================================================

def ordres_rt(objet):
    """Ensemble des n apparaissant dans les `rt[n]` d'une structure lue."""
    if isinstance(objet, str):
        return {int(x) for x in re.findall(r'rt\[(\d+)\]', objet)}
    out = set()
    for x in objet:
        out |= ordres_rt(x)
    return out


def matrice_mod_p(mat, p, racines):
    """Matrice de chaines -> matrice d'entiers de GF(p)."""
    from cy_landscape.core.gamma_action import parse_entree
    return [[parse_entree(x, p, racines) for x in ligne] for ligne in mat]
