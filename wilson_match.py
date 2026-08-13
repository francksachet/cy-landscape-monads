#!/usr/bin/env python3
"""
wilson_match.py -- Croisement de la liste de Braun (195 quotients libres)
avec cicylist.txt et avec un catalogue de resultats.

--------------------------------------------------------------------------
A quoi ca sert
--------------------------------------------------------------------------
Une ligne de Wilson exige un quotient par une symetrie librement agissante.
Sans une telle symetrie, un fibre ne peut pas servir a briser E6 vers le
Modele Standard, quelle que soit la qualite de son spectre.

Or seules 195 des 7890 CICYs en possedent une (Braun, JHEP 1104 (2011) 005 ;
donnees de la page CicyQuotients d'Oxford). Ce script determine lesquelles
de TES CICYs -- celles qui portent des candidats -- sont dans cette liste.

--------------------------------------------------------------------------
Pourquoi on ne peut pas croiser sur les numeros
--------------------------------------------------------------------------
Les deux fichiers ne numerotent pas de la meme facon. Dans le fichier de
Braun, l'entree `Num -> 4` a une matrice 8x5 et (h11, h21) = (15, 15) ;
dans cicylist.txt, l'entree `Num : 1` a une matrice 7x6 et les memes
nombres de Hodge. Ce sont deux varietes differentes portant le meme couple
(15, 15) -- il y en a beaucoup.

Le croisement se fait donc sur la MATRICE DE CONFIGURATION, qui identifie
reellement la variete. Difficulte : elle n'est definie qu'a permutation
pres des lignes (facteurs projectifs) et des colonnes (polynomes). On
compare donc des formes canoniques.

Strategie : filtrer d'abord par invariants bon marche (dimensions, sommes
de lignes et de colonnes triees, multiensemble des lignes), puis verifier
les rares survivants par recherche exhaustive de permutation.
"""
import re, sys, json, os, argparse
from itertools import permutations
from collections import defaultdict


# ======================================================================
# Lecture du fichier Mathematica
# ======================================================================

def parse_braun(path):
    """
    Extrait [{num, conf, h11, h21, symetries:[noms]}] du fichier .m.

    Le fichier est du Mathematica : on ne l'interprete pas, on repere les
    motifs `Num -> n`, `Conf -> {{...}}` et les noms de groupes entre
    guillemets a l'interieur du bloc Symmetries.
    """
    txt = open(path, encoding='utf-8', errors='replace').read()
    # On ne garde que la liste elle-meme, apres l'affectation
    i = txt.find('CicyQuotients=')
    if i >= 0:
        txt = txt[i + len('CicyQuotients='):]

    entrees = []
    for m in re.finditer(r'Num\s*->\s*(\d+)', txt):
        debut = m.start()
        fin = txt.find('Num ->', m.end())
        if fin == -1:
            fin = len(txt)
        bloc = txt[debut:fin]

        conf = _extraire_conf(bloc)
        if conf is None:
            continue
        h11 = _premier_int(bloc, r'H11\s*->\s*(-?\d+)')
        h21 = _premier_int(bloc, r'H21\s*->\s*(-?\d+)')
        noms = re.findall(r'"([A-Za-z0-9_\\\[\]rt\.\*x ]+)"', bloc)
        noms = [n for n in noms if not n.startswith('\\')]
        entrees.append({'num': int(m.group(1)), 'conf': conf,
                        'h11': h11, 'h21': h21, 'symetries': noms})
    return entrees


def _premier_int(bloc, motif):
    m = re.search(motif, bloc)
    return int(m.group(1)) if m else None


def _extraire_conf(bloc):
    """Lit la matrice qui suit `Conf ->` en equilibrant les accolades."""
    m = re.search(r'Conf\s*->\s*\{', bloc)
    if not m:
        return None
    i = m.end() - 1
    prof = 0
    for j in range(i, len(bloc)):
        if bloc[j] == '{':
            prof += 1
        elif bloc[j] == '}':
            prof -= 1
            if prof == 0:
                brut = bloc[i:j + 1]
                break
    else:
        return None
    lignes = re.findall(r'\{([^{}]*)\}', brut)
    conf = []
    for l in lignes:
        vals = [int(x) for x in re.findall(r'-?\d+', l)]
        if vals:
            conf.append(vals)
    return conf or None


# ======================================================================
# Lecture de cicylist.txt
# ======================================================================

def parse_cicylist(path):
    """
    Format Oxford texte : blocs `Num : n`, `NumPs`, `NumPol`, `H11`, `H21`,
    puis NumPs lignes `{...}` formant la matrice de configuration.
    """
    entrees = []
    cur = None
    lignes_mat = []
    for ligne in open(path, encoding='utf-8', errors='replace'):
        s = ligne.strip()
        if not s:
            continue
        m = re.match(r'Num\s*:\s*(\d+)', s)
        if m:
            if cur is not None and lignes_mat:
                cur['conf'] = lignes_mat
                entrees.append(cur)
            cur = {'num': int(m.group(1))}
            lignes_mat = []
            continue
        if cur is None:
            continue
        for cle, motif in (('numps', r'NumPs\s*:\s*(-?\d+)'),
                           ('numpol', r'NumPol\s*:\s*(-?\d+)'),
                           ('h11', r'H11\s*:\s*(-?\d+)'),
                           ('h21', r'H21\s*:\s*(-?\d+)')):
            m2 = re.match(motif, s)
            if m2:
                cur[cle] = int(m2.group(1))
        if s.startswith('{') and not re.match(r'(C2|Redun)\s*:', s):
            vals = [int(x) for x in re.findall(r'-?\d+', s)]
            if vals and len(lignes_mat) < cur.get('numps', 10 ** 6):
                lignes_mat.append(vals)
    if cur is not None and lignes_mat:
        cur['conf'] = lignes_mat
        entrees.append(cur)
    return entrees


# ======================================================================
# Forme canonique d'une matrice de configuration
# ======================================================================

def invariants(conf):
    """Invariants bon marche sous permutation des lignes et des colonnes."""
    nl = len(conf)
    nc = len(conf[0]) if nl else 0
    sl = tuple(sorted(sum(r) for r in conf))
    sc = tuple(sorted(sum(conf[i][j] for i in range(nl)) for j in range(nc)))
    lignes = tuple(sorted(tuple(sorted(r)) for r in conf))
    return (nl, nc, sl, sc, lignes)


def _transposee(conf):
    return [list(c) for c in zip(*conf)]


def _raffiner(conf):
    """
    Raffinement iteratif facon Weisfeiler-Leman : on attribue a chaque ligne
    et a chaque colonne une couleur, initialisee par ses sommes, puis
    raffinee par le multiensemble des (valeur, couleur de l'autre indice).
    Renvoie (couleurs_lignes, couleurs_colonnes), invariantes par permutation.
    """
    nl, nc = len(conf), len(conf[0])
    cl = [sum(r) for r in conf]
    cc = [sum(conf[i][j] for i in range(nl)) for j in range(nc)]
    for _ in range(nl + nc + 2):
        nl_new = [ (cl[i], tuple(sorted((conf[i][j], cc[j]) for j in range(nc))))
                   for i in range(nl) ]
        nc_new = [ (cc[j], tuple(sorted((conf[i][j], cl[i]) for i in range(nl))))
                   for j in range(nc) ]
        rl = {v: k for k, v in enumerate(sorted(set(nl_new)))}
        rc = {v: k for k, v in enumerate(sorted(set(nc_new)))}
        nl2 = [rl[v] for v in nl_new]
        nc2 = [rc[v] for v in nc_new]
        if nl2 == cl and nc2 == cc:
            break
        cl, cc = nl2, nc2
    return cl, cc


def meme_matrice(a, b, max_total=200000):
    """
    Egalite a permutation pres des lignes ET des colonnes.

    Le raffinement de couleurs partitionne lignes et colonnes en classes ;
    seules les permutations preservant les classes peuvent convenir, ce qui
    fait tomber le cout de nc! a un produit de factorielles de tailles de
    classes. C'est ce qui manquait a la version precedente, qui se rabattait
    sur les invariants des que nc > 8 et declarait alors des correspondances
    ambigues (jusqu'a 78 candidates pour une seule entree de Braun).
    """
    if len(a) != len(b) or len(a[0]) != len(b[0]):
        return False
    cla, cca = _raffiner(a)
    clb, ccb = _raffiner(b)
    if sorted(cla) != sorted(clb) or sorted(cca) != sorted(ccb):
        return False

    nc = len(a[0])
    # classes de colonnes de b, par couleur
    from collections import defaultdict as _dd
    cls_b = _dd(list)
    for j, c in enumerate(ccb):
        cls_b[c].append(j)
    cls_a = _dd(list)
    for j, c in enumerate(cca):
        cls_a[c].append(j)
    if sorted((k, len(v)) for k, v in cls_a.items()) != \
       sorted((k, len(v)) for k, v in cls_b.items()):
        return False

    total = 1
    for v in cls_a.values():
        total *= _fact(len(v))
        if total > max_total:
            # partition trop grossiere : on tranche par les invariants seuls
            return invariants(a) == invariants(b)

    cible = sorted(tuple(r) for r in b)
    couleurs = sorted(cls_a)
    listes = [list(permutations(cls_b[c])) for c in couleurs]
    for combo in _produit(listes):
        perm = [0] * nc
        for c, ordre in zip(couleurs, combo):
            for src, dst in zip(cls_a[c], ordre):
                perm[src] = dst
        cand = sorted(tuple(r[perm.index(k)] for k in range(nc)) for r in a)
        if cand == cible:
            return True
    return False


def _fact(n):
    r = 1
    for i in range(2, n + 1):
        r *= i
    return r


def _produit(listes):
    from itertools import product as _p
    return _p(*listes) if listes else [()]


def apparier(braun, cicy):
    """
    Renvoie {num_braun: num_cicylist}. Essaie la matrice telle quelle puis
    sa transposee (les deux sources peuvent ne pas orienter pareil).
    """
    index = defaultdict(list)
    for e in cicy:
        index[invariants(e['conf'])].append(e)

    corresp = {}
    ambigus = []
    absents = []
    for b in braun:
        trouve = None
        for conf in (b['conf'], _transposee(b['conf'])):
            cands = index.get(invariants(conf), [])
            exacts = [c for c in cands if meme_matrice(conf, c['conf'])]
            if len(exacts) == 1:
                trouve = exacts[0]['num']; break
            if len(exacts) > 1:
                ambigus.append((b['num'], [c['num'] for c in exacts]))
                trouve = exacts[0]['num']; break
        if trouve is None:
            absents.append(b['num'])
        else:
            corresp[b['num']] = trouve
    return corresp, ambigus, absents


# ======================================================================
# Programme principal
# ======================================================================

_ORDRES_BASE = {'Z': None}


def _ordre_groupe(nom):
    """
    Ordre du groupe a partir de son nom ('Z3', 'Z3 x Z3', 'Z2', 'Q8', ...).

    Les noms cycliques et leurs produits directs se lisent directement. Pour
    les autres, on retombe sur une table courte ; un nom inconnu renvoie
    None et la CICY est alors exportee sans ordre exploitable.
    """
    table = {'Q8': 8, 'D4': 8, 'A4': 12, 'S3': 6, 'S4': 24, 'Dic3': 12,
             'Dic5': 20, 'H8': 8, 'F20': 20}
    n = nom.strip()
    if n in table:
        return table[n]
    facteurs = [f.strip() for f in re.split(r'[x\*]', n)]
    total = 1
    for f in facteurs:
        m = re.fullmatch(r'Z(\d+)', f)
        if not m:
            return None
        total *= int(m.group(1))
    return total if total > 1 else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('braun_m', help="fichier CicyQuotients .m")
    ap.add_argument('cicylist', help="cicylist.txt")
    ap.add_argument('--results', default=None,
                    help="results.jsonl (ou results_clean.jsonl) a croiser")
    ap.add_argument('--export', default='wilson_cicys.json',
                    help="fichier de sortie lu par le pipeline (defaut: "
                         "wilson_cicys.json). Contient, pour chaque CICY de "
                         "TA numerotation, les groupes disponibles et leurs "
                         "ordres.")
    args = ap.parse_args()

    braun = parse_braun(args.braun_m)
    cicy = parse_cicylist(args.cicylist)
    print(f"\n{'='*70}")
    print(f"  Entrees lues : {len(braun)} chez Braun, {len(cicy)} dans cicylist.txt")
    print(f"{'='*70}")

    corresp, ambigus, absents = apparier(braun, cicy)
    print(f"  Apparieees sur la matrice de configuration : {len(corresp)}/{len(braun)}")
    if ambigus:
        print(f"  Ambigues (plusieurs correspondances) : {len(ambigus)}")
        for b, l in ambigus[:5]:
            print(f"    Braun #{b} -> {l}")
    if absents:
        print(f"  Sans correspondance : {len(absents)} -> {absents[:12]}")

    # --- ordres des groupes ---------------------------------------------
    # n_gen(X/Gamma) = n_gen(X) / |Gamma| : pour obtenir 3 generations sur le
    # quotient, il faut 3*|Gamma| generations en amont. On exporte donc les
    # ordres disponibles pour que le prefiltre vise la bonne valeur.
    par_cicy = {}
    for b in braun:
        n = corresp.get(b['num'])
        if n is None:
            continue
        ordres = sorted({_ordre_groupe(s) for s in b['symetries']} - {None})
        e = par_cicy.setdefault(n, {'braun': b['num'], 'groupes': [], 'ordres': []})
        for g in b['symetries']:
            if g not in e['groupes']:
                e['groupes'].append(g)
        for o in ordres:
            if o not in e['ordres']:
                e['ordres'].append(o)
    for e in par_cicy.values():
        e['ordres'] = sorted(e['ordres'])

    with open(args.export, 'w', encoding='utf-8') as f:
        json.dump({str(k): v for k, v in sorted(par_cicy.items())}, f,
                  ensure_ascii=False, indent=1)
    print(f"  Export des groupes -> {args.export}")

    libres = sorted(set(corresp.values()))
    print(f"\n  CICYs de TA numerotation admettant une symetrie librement")
    print(f"  agissante : {len(libres)}")
    print(f"    {libres}")

    if args.results and os.path.exists(args.results):
        rs = [json.loads(l) for l in open(args.results, encoding='utf-8') if l.strip()]
        porteuses = sorted({r['cicy'] for r in rs})
        inter = sorted(set(porteuses) & set(libres))
        print(f"\n{'='*70}")
        print(f"  CROISEMENT AVEC {args.results}")
        print(f"{'='*70}")
        print(f"  Candidats            : {len(rs)}")
        print(f"  CICYs porteuses      : {len(porteuses)}")
        print(f"  Dont a symetrie libre: {len(inter)}   {inter}")

        util = [r for r in rs if r['cicy'] in set(inter)]
        print(f"  Candidats exploitables pour une ligne de Wilson : {len(util)}")
        if util:
            inv = {v: k for k, v in corresp.items()}
            print(f"\n    {'#CICY':>6} {'jauge':>7} {'rk':>2} {'n_anti':>6} "
                  f"{'cohomologie':>16}  symetries")
            for r in sorted(util, key=lambda x: min(x['cohomology'][1],
                                                    x['cohomology'][2]))[:30]:
                h = r['cohomology']
                bnum = inv.get(r['cicy'])
                sym = next((', '.join(sorted(set(b['symetries'])))
                            for b in braun if b['num'] == bnum), '')
                print(f"    {r['cicy']:>6} {r.get('gauge',''):>7} "
                      f"{r.get('rank_V',''):>2} {min(h[1],h[2]):>6} "
                      f"{str(h):>16}  {sym[:40]}")
    print()


if __name__ == '__main__':
    sys.exit(main())
