#!/usr/bin/env python3
"""
portee_substitution.py -- Quelles lignes de verdict le defaut du 5.34 a-t-il
faussees, exactement ?

Le 5.34 a corrige `matrice_substitution` : ses colonnes sortaient dans
l'ordre (sigma^-1(0), ...) alors que `basis_multi` les indexe dans l'ordre
(0, 1, ...). Les deux coincident SI ET SEULEMENT SI sigma est l'identite.
Tout verdict calcule sur une realisation dont sigma PERMUTE deux facteurs de
l'ambiant est donc a jeter.

CE QUE CE SCRIPT CORRIGE DANS LA MESURE DU 5.34
-----------------------------------------------
Le 5.34 annonce « 27 couples (CICY, groupe), 1 224 lignes, 213 SURVIT ». Ce
compte partitionne par NOM DE GROUPE. Or sigma est une propriete de la
REALISATION, pas du nom : une meme CICY peut porter dix realisations de Z2
dont une seule permute. Les couples MIXTES tombent alors du cote « intact »,
et leurs lignes permutantes avec.

Le script mesure les deux granularites et les rapproche :

    couples dont TOUTES les realisations permutent   -> les « 27 » du 5.34
    couples MIXTES (au moins une de chaque)          -> comptes « intacts »
    couples sans aucune realisation permutante       -> reellement intacts
    pseudo-couples (CICY, '-')                       -> aucun calcul

La somme des trois premiers doit rendre les chiffres du 5.34 au chiffre pres.
Si elle ne les rend pas, le script S'ARRETE : ces trois nombres sont la seule
reference EXTERIEURE dont on dispose pour valider la classification, et une
classification qui ne les retrouve pas est fausse, pas interessante.

CE QU'IL PRODUIT
----------------
La liste des LOTS a recalculer -- un lot est touche si sa tranche contient au
moins une realisation permutante. Le lot, et non la ligne : c'est l'unite que
parle le checkpoint, et c'est exact, alors qu'attribuer une ligne a une
realisation demanderait de realigner a posteriori l'ordre dans lequel
`analyser` a rendu ses lignes.

Consequence assumee : un lot touche emmene avec lui ses realisations non
permutantes. Elles seront recalculees pour rien -- et c'est un CONTROLE : le
correctif ne doit RIEN changer sur elles. Toute divergence sur ces lignes-la
signifie que la classification de sigma est fausse, ou qu'autre chose a bouge.

REGLE DES FILTRES (section 8)
-----------------------------
Ce script est lui-meme un filtre : il decide ce qu'on jette. Il declare donc
les DEUX cotes -- ce qu'il retient, ce qu'il ecarte, et ce qu'il ne sait pas
separer (les lignes des couples mixtes, dont une partie seulement est
fausse). Et il refuse de servir si ses trois ancres ne tombent pas juste.

GARDE SUR LES IDENTIFIANTS DE LOT
---------------------------------
Un identifiant de lot est ('T', k, t) ou k est le RANG DE LA TACHE dans `rs`.
Il ne veut rien dire hors de la configuration exacte du run qui l'a produit :
--cicy, --replier-orbites, --controle-orbites et --taille-lot entrent dans
l'empreinte du checkpoint. Le script recalcule cette empreinte et REFUSE de
continuer si elle ne correspond pas a celle du checkpoint : sans cela, il
designerait des lots parfaitement valides comme etant a jeter.

Usage (lecture seule, n'ecrit que son rapport) :
    python -u portee_substitution.py cicyquotients.m cicylist.txt scan_wilson4 \
           --max-realisations 16 --replier-orbites --sortie portee_5_34.json
"""
import os
import sys
import json
import argparse
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import equivariance_f as EF
from cy_landscape.core.braun_symmetry import (parse_symmetries, ordres_rt,
                                              matrice_mod_p)
from cy_landscape.core.gamma_action import choisir_premier, racine_primitive
from cy_landscape.core.covariant_ring import permutation_facteurs_numerique
from cy_landscape.data.parse_oxford import load_oxford_file
from wilson_match import parse_braun, parse_cicylist, apparier

# Les trois ancres du 5.34, cote « a recalculer » et cote « intact ».
ANCRES = {
    'recalculer': {'couples': 27, 'lignes': 1224, 'survit': 213},
    'intact':     {'couples': 102, 'lignes': 173623, 'survit': 32886},
}

ID, PERM, INEXTRACTIBLE, ILLISIBLE = 'id', 'permute', 'non extractible', 'illisible'


# ======================================================================
# 1. sigma, realisation par realisation
# ======================================================================

def statut_sigma(sym, ambient):
    """
    Statut d'UNE realisation : 'id' si tous ses generateurs fixent chaque
    facteur, 'permute' si l'un d'eux permute, 'non extractible' si
    `permutation_facteurs_numerique` refuse (la matrice melangerait des
    facteurs), 'illisible' si les generateurs ne se reduisent pas mod p.

    Les deux derniers cas ne sont PAS ranges avec 'id'. Un statut inconnu
    qu'on compte comme sain est exactement le defaut que la section 8
    appelle « remplacer par zero ce qui n'est pas calcule ».
    """
    ordres = ordres_rt(sym['coord']) | ordres_rt(sym['poly']) | {2}
    try:
        p, _ = choisir_premier(sorted(ordres), minimum=30011)
        rac = {n: racine_primitive(p, n) for n in ordres}
        Mc = [matrice_mod_p(x, p, rac) for x in sym['coord']]
    except Exception:
        return ILLISIBLE
    identite = list(range(len(ambient)))
    sigmas = [permutation_facteurs_numerique(M, ambient, p) for M in Mc]
    if any(s is None for s in sigmas):
        return INEXTRACTIBLE
    return PERM if any(s != identite for s in sigmas) else ID


def classer_realisations(entries, inv, SYM):
    """{cicy: [[indice, nom, statut], ...]} sur toute la liste de Braun."""
    out = {}
    for cicy, num_b in sorted(inv.items()):
        if num_b not in SYM or cicy not in entries:
            continue
        amb = entries[cicy]['ambient']
        out[cicy] = [[n, sy['nom'], statut_sigma(sy, amb)]
                     for n, sy in enumerate(SYM[num_b]['symetries'])]
    return out


# ======================================================================
# 2. reconstruction des lots -- la meme boucle que `equivariance_f.main`
# ======================================================================

def reconstruire_lots(rs, entries, inv, SYM, taches, controles, args,
                      plafond=None):
    """
    Rend [(identifiant, cicy, [indices de realisation])].

    `plafond` = --max-realisations a simuler ; None reprend celui de `args`.
    0 = aucun plafond, c'est-a-dire le SUR-ENSEMBLE de tous les lots que ce
    balayage peut produire.

    LE SUR-ENSEMBLE EST LA BONNE BASE D'INTERPRETATION. Le plafond n'entre
    pas dans l'empreinte du checkpoint -- volontairement, pour que le lever
    n'invalide rien et ne fasse qu'ajouter des lots. La consequence est
    qu'un checkpoint peut MELANGER des lots venus de runs a plafonds
    differents : celui de `scan_wilson4` contient un ('T', 823, 5) que la
    reconstruction a plafond 16 ne connait pas, et toute sa zone heritee a
    ete calculee sans plafond du tout. Interpreter les identifiants avec le
    plafond du prochain run declarerait donc « inexistants » des lots bien
    reels, et laisserait leurs lignes en place sans les voir.

    C'est la boucle de `main` sous une autre forme -- une duplication, donc
    un risque : un controle qui partage le defaut de l'objet controle ne
    controle rien (5.34). Deux garde-fous exterieurs la rattachent au reel :
    l'empreinte du checkpoint, et l'exigence que TOUT lot du checkpoint
    existe dans la reconstruction (`verifier_reconstruction`).
    """
    if plafond is None:
        plafond = args.max_realisations
    lots = []
    for k, (i_rep, membres) in enumerate(taches):
        r = rs[i_rep]
        num_b = inv.get(r['cicy'])
        if entries.get(r['cicy']) is None or num_b not in SYM:
            continue
        g = None if args.tous_groupes else set(r.get('groupes_utiles') or [])
        if g is not None and not g:
            lots.append((('T', k, -1), r['cicy'], []))
            continue
        idx = [n for n, sy in enumerate(SYM[num_b]['symetries'])
               if (not g) or sy['nom'] in g]
        if plafond:
            # plafond ARRONDI AU MULTIPLE SUPERIEUR de la taille de lot,
            # comme dans `main` : sans cet arrondi les tranches ne seraient
            # plus les memes et les identifiants designeraient autre chose.
            pas = max(1, args.taille_lot)
            idx = idx[:-(-plafond // pas) * pas]
        tranches = [idx[t:t + args.taille_lot]
                    for t in range(0, len(idx), args.taille_lot)] or [[]]
        for t, tr in enumerate(tranches):
            lots.append((('T', k, t), r['cicy'], tr))
        for j in controles.get(k, ()):
            for t, tr in enumerate(tranches):
                lots.append((('C', k, j, t), rs[j]['cicy'], tr))
    return lots


def verifier_reconstruction(lots, faits):
    """
    Tout lot du checkpoint doit exister dans la reconstruction. Un seul
    manquant, et les identifiants ne designent pas les memes calculs : on
    s'arrete.

    Les entrees a deux elements ('T', k) sont exclues du test -- elles
    viennent d'un checkpoint sequentiel migre et ne portent pas de rang de
    tranche.
    """
    ids = {l[0] for l in lots}
    absents = [f for f in faits if len(f) > 2 and f not in ids]
    return absents


# ======================================================================
# 3. bilan par couple (CICY, groupe) -- les trois ancres
# ======================================================================

def bilan_couples(chemin_jsonl, statuts):
    """
    Compte lignes et SURVIT par couple, puis range chaque couple dans une
    des quatre categories. Rend (categories, lignes_par_couple, survit_par_couple).
    """
    lignes = Counter()
    survit = Counter()
    par_lot = Counter()
    sans_lot = 0
    with open(chemin_jsonl, encoding='utf-8') as f:
        for L in f:
            if not L.strip():
                continue
            d = json.loads(L)
            cle = (d.get('cicy'), d.get('groupe'))
            lignes[cle] += 1
            survit[cle] += bool(d.get('survit'))
            lot = d.get('_lot')
            if lot is None:
                sans_lot += 1
            else:
                par_lot[tuple(lot)] += 1

    # statut par (cicy, nom de groupe) : le multi-ensemble des statuts de
    # ses realisations.
    par_couple = {}
    for cicy, rows in statuts.items():
        for _, nom, st in rows:
            par_couple.setdefault((cicy, nom), Counter())[st] += 1

    cats = {'recalculer': [], 'mixte': [], 'intact': [], 'pseudo': []}
    for cle in lignes:
        cnt = par_couple.get(cle)
        if cnt is None:
            cats['pseudo'].append(cle)          # (cicy, '-') : aucun calcul
        elif set(cnt) == {PERM}:
            cats['recalculer'].append(cle)
        elif PERM in cnt:
            cats['mixte'].append(cle)
        else:
            cats['intact'].append(cle)
    return cats, lignes, survit, par_lot, sans_lot


def somme(cles, lignes, survit):
    return {'couples': len(cles),
            'lignes': sum(lignes[c] for c in cles),
            'survit': sum(survit[c] for c in cles)}


# ======================================================================

def principal():
    ap = argparse.ArgumentParser()
    ap.add_argument('braun_m')
    ap.add_argument('cicylist')
    ap.add_argument('output_dir')
    ap.add_argument('--input', default='results_equivariant.jsonl')
    ap.add_argument('--sortie', default='portee_5_34.json')
    # Les options qui DEFINISSENT les identifiants de lot. Elles doivent
    # reproduire celles du run qu'on repare, et l'empreinte le verifie.
    ap.add_argument('--cicy', type=int, default=None)
    ap.add_argument('--tous-groupes', action='store_true')
    ap.add_argument('--replier-orbites', action='store_true')
    ap.add_argument('--controle-orbites', type=int, default=20)
    ap.add_argument('--taille-lot', type=int, default=16)
    ap.add_argument('--max-realisations', type=int, default=0)
    ap.add_argument('--sans-ancres', action='store_true',
                    help="ne pas s'arreter si les chiffres du 5.34 ne sont "
                         "pas reproduits (a n'utiliser que sur un AUTRE "
                         "balayage, ou ils n'ont aucune raison de tomber)")
    args = ap.parse_args()

    src = os.path.join(args.output_dir, args.input)
    dst = os.path.join(args.output_dir, 'results_equivariance_f.jsonl')
    chemin_ck = os.path.join(args.output_dir, 'progress_equivariance_f.json')
    for p in (src, dst, chemin_ck):
        if not os.path.exists(p):
            sys.exit(f"  ABSENT : {p}")

    print("\n  chargement de la liste CICY et des symetries de Braun ...")
    entries = {e['num']: e for e in load_oxford_file(args.cicylist)}
    braun = parse_braun(args.braun_m)
    cl = parse_cicylist(args.cicylist)
    corr, _, _ = apparier(braun, cl)
    inv = {v: k for k, v in corr.items()}
    SYM = parse_symmetries(args.braun_m)

    # ---- garde n.1 : l'empreinte -------------------------------------
    empreinte = EF._empreinte(src, (args.cicy, bool(args.replier_orbites),
                                    int(args.controle_orbites),
                                    int(args.taille_lot)))
    with open(chemin_ck, encoding='utf-8') as f:
        ck = json.load(f)
    print(f"\n  empreinte du checkpoint : {ck.get('empreinte')}")
    print(f"  empreinte des options   : {empreinte}")
    if ck.get('empreinte') != empreinte:
        sys.exit("\n  ARRET. Les options passees ne sont pas celles du run qui a\n"
                 "  produit ce checkpoint. Les identifiants de lot ('T', k, t)\n"
                 "  designeraient d'autres calculs, et le script marquerait a\n"
                 "  jeter des lots valides. Verifier --cicy, --replier-orbites,\n"
                 "  --controle-orbites et --taille-lot.")
    print("  -> les identifiants de lot sont ceux de ce run.")
    faits = {tuple(x[0]): x[1] for x in ck['lots']}
    heritees = {f[1] for f in faits if len(f) == 2}

    rs = [json.loads(l) for l in open(src, encoding='utf-8') if l.strip()]
    if args.cicy:
        rs = [r for r in rs if r['cicy'] == args.cicy]
    taches, rapport = EF._construire_taches(rs, entries, args)
    controles = EF._echantillon_controle(taches, args.controle_orbites)
    print(f"\n  {len(rs)} candidats, {len(taches)} taches"
          + (f", repli x{rapport['facteur']:.2f}" if rapport else ""))

    print("\n  classification de sigma, realisation par realisation ...")
    statuts = classer_realisations(entries, inv, SYM)
    tous = Counter(st for rows in statuts.values() for _, _, st in rows)
    print(f"    {sum(tous.values())} realisations dans la liste de Braun : "
          + ", ".join(f"{k} = {v}" for k, v in sorted(tous.items())))

    # SUR-ENSEMBLE (aucun plafond) = base d'interpretation des identifiants.
    # CIBLE (plafond du prochain run) = ce que la relance recalculera.
    lots = reconstruire_lots(rs, entries, inv, SYM, taches, controles, args,
                             plafond=0)
    ids_cible = {l[0] for l in reconstruire_lots(
        rs, entries, inv, SYM, taches, controles, args,
        plafond=args.max_realisations)}
    absents = verifier_reconstruction(lots, faits)
    print(f"\n  {len(lots)} lots au total (aucun plafond), "
          f"{len(ids_cible)} dans la cible "
          f"(--max-realisations {args.max_realisations or 'aucun'}), "
          f"{len(faits)} dans le checkpoint")
    if absents:
        sys.exit(f"\n  ARRET. {len(absents)} lots du checkpoint n'existent pas dans\n"
                 f"  la reconstruction, par exemple {absents[:3]}. La boucle de\n"
                 f"  reconstruction ne reproduit pas celle du run.")
    print("  -> tout lot du checkpoint existe dans la reconstruction.")
    hors_cible = [f for f in faits if len(f) > 2 and f not in ids_cible]
    if hors_cible:
        print(f"  {len(hors_cible)} lots deja faits sont HORS de la cible : "
              f"la relance ne les rejouera pas.")

    # ---- garde n.2 : les trois ancres du 5.34 ------------------------
    cats, lignes, survit, par_lot, sans_lot = bilan_couples(dst, statuts)
    b = {k: somme(v, lignes, survit) for k, v in cats.items()}
    intact_5_34 = {c: b['intact'][c] + b['mixte'][c] + b['pseudo'][c]
                   for c in ('couples', 'lignes', 'survit')}
    intact_5_34['couples'] = b['intact']['couples'] + b['mixte']['couples']

    print(f"\n{'=' * 78}")
    print("  BILAN PAR COUPLE (CICY, groupe)")
    print(f"{'=' * 78}")
    print(f"  {'categorie':<44} {'couples':>8} {'lignes':>8} {'SURVIT':>8}")
    for cle, titre in (
            ('intact', "aucune realisation ne permute -- intacts"),
            ('mixte', "MIXTE : au moins une de chaque -- A RECALCULER EN PARTIE"),
            ('pseudo', "pseudo-couples (CICY, '-') : aucun calcul"),
            ('recalculer', "toutes permutent -- les « 27 » du 5.34")):
        v = b[cle]
        print(f"  {titre:<44} {v['couples']:>8} {v['lignes']:>8} {v['survit']:>8}")
    print(f"  {'-' * 76}")
    print(f"  {'TOTAL':<44} {sum(v['couples'] for v in b.values()):>8} "
          f"{sum(v['lignes'] for v in b.values()):>8} "
          f"{sum(v['survit'] for v in b.values()):>8}")

    ok_r = all(b['recalculer'][c] == ANCRES['recalculer'][c]
               for c in ('couples', 'lignes', 'survit'))
    ok_i = all(intact_5_34[c] == ANCRES['intact'][c]
               for c in ('couples', 'lignes', 'survit'))
    print(f"\n  ancre 5.34 « a recalculer » {ANCRES['recalculer']} : "
          f"{'RETROUVEE' if ok_r else 'NON RETROUVEE -> ' + str(b['recalculer'])}")
    print(f"  ancre 5.34 « intacts »      {ANCRES['intact']} : "
          f"{'RETROUVEE' if ok_i else 'NON RETROUVEE -> ' + str(intact_5_34)}")
    if not (ok_r and ok_i) and not args.sans_ancres:
        sys.exit("\n  ARRET. Les chiffres du 5.34 ne sont pas reproduits : la\n"
                 "  classification de sigma est fausse, ou le fichier a change.\n"
                 "  --sans-ancres pour forcer (sur un autre balayage seulement).")

    # ---- ce que le 5.34 n'avait pas compte --------------------------
    print(f"\n  Les {b['mixte']['couples']} couples MIXTES etaient ranges du cote "
          f"« intact » :")
    for cle in sorted(cats['mixte']):
        cnt = Counter(st for _, nom, st in statuts[cle[0]] if nom == cle[1])
        print(f"    #{cle[0]:<5} {cle[1]:<11} {dict(cnt)}  "
              f"{lignes[cle]:>6} lignes, {survit[cle]:>5} SURVIT")

    # ---- les lots touches -------------------------------------------
    st_par_cicy = {c: {n: st for n, _, st in rows} for c, rows in statuts.items()}
    touches, vues, permutantes = [], 0, 0
    for ident, cicy, tr in lots:
        stp = st_par_cicy.get(cicy, {})
        npm = sum(1 for n in tr if stp.get(n) == PERM)
        if npm:
            touches.append((ident, cicy, tr, npm))
            vues += len(tr)
            permutantes += npm
    ids_t = [t for t in touches if t[0][0] == 'T']
    ids_c = [t for t in touches if t[0][0] == 'C']
    expl = [t for t in touches if t[0] in faits]
    couv = [t for t in touches if t[0] not in faits and ('T', t[0][1]) in faits]
    pend = [t for t in touches if t[0] not in faits and ('T', t[0][1]) not in faits]
    deja = expl + couv

    print(f"\n{'=' * 78}")
    print("  LOTS A RECALCULER (tranche contenant au moins une realisation "
          "permutante)")
    print(f"{'=' * 78}")
    print(f"  {len(touches)} lots : {len(ids_t)} representants, {len(ids_c)} controles")
    print(f"    deja faits, identifiant explicite      : {len(expl)}")
    print(f"    deja faits via une entree heritee      : {len(couv)}")
    print(f"    pas encore faits (rien a retirer)      : {len(pend)}")
    print(f"  {vues} realisations dans ces tranches, dont {permutantes} "
          f"permutent reellement")
    print(f"    les {vues - permutantes} autres seront recalculees pour rien "
          f"-- et doivent revenir IDENTIQUES (controle)")

    # REGLE DES FILTRES : ce que la cible ne rattrape pas doit etre dit.
    # Un lot touche, deja fait, mais hors de la cible du prochain run ne
    # sera JAMAIS recalcule : ses lignes fausses resteraient en place si on
    # ne les retirait pas, et la couverture correspondante disparait.
    dedans = [t for t in deja if t[0] in ids_cible]
    dehors = [t for t in deja if t[0] not in ids_cible]
    print(f"\n  parmi les {len(deja)} lots touches DEJA FAITS :")
    print(f"    dans la cible -- retires puis RECALCULES : {len(dedans)} lots, "
          f"{sum(len(t[2]) for t in dedans)} realisations, "
          f"{sum(t[3] for t in dedans)} permutantes")
    print(f"    hors cible    -- retires et NON recalcules : {len(dehors)} lots, "
          f"{sum(len(t[2]) for t in dehors)} realisations, "
          f"{sum(t[3] for t in dehors)} permutantes")
    if dehors:
        print(f"    ces realisations-la retombent dans les non testees. "
              f"Le champ `realisations`")
        print(f"    des lignes qui survivent sur ces taches sur-declare alors "
              f"la couverture.")

    # lignes a jeter, declarees des deux cotes
    lig_expl = sum(faits.get(t[0], 0) or 0 for t in expl)
    taches_heritees = sorted({t[0][1] for t in couv})
    cles_heritees = set()
    conflits = Counter()
    for k in taches_heritees:
        for j in taches[k][1]:
            cle = (rs[j]['cicy'], json.dumps(rs[j]['b_charges']),
                   json.dumps(rs[j]['c_charges']))
            cles_heritees.add(cle)
            conflits[cle] += 1
    # une meme identite de candidat ne doit appartenir qu'a une tache
    doublons = {c: n for c, n in conflits.items() if n > 1}
    lig_her = 0
    surv_her = 0
    if cles_heritees:
        with open(dst, encoding='utf-8') as f:
            for L in f:
                if not L.strip():
                    continue
                d = json.loads(L)
                if d.get('_lot') is not None:
                    continue
                cle = (d.get('cicy'), json.dumps(d.get('b_charges')),
                       json.dumps(d.get('c_charges')))
                if cle in cles_heritees:
                    lig_her += 1
                    surv_her += bool(d.get('survit'))
    surv_expl = 0
    ens_expl = {t[0] for t in expl}
    with open(dst, encoding='utf-8') as f:
        for L in f:
            if not L.strip():
                continue
            d = json.loads(L)
            lot = d.get('_lot')
            if lot is not None and tuple(lot) in ens_expl:
                surv_expl += bool(d.get('survit'))

    print(f"\n  LIGNES A JETER")
    print(f"    dans les lots a identifiant explicite : {lig_expl:>6}, "
          f"dont {surv_expl} SURVIT")
    print(f"    lignes heritees sans `_lot` ({len(taches_heritees)} taches) : "
          f"{lig_her:>6}, dont {surv_her} SURVIT")
    print(f"    TOTAL : {lig_expl + lig_her} lignes, "
          f"{surv_expl + surv_her} SURVIT")
    print(f"    (le 5.34 en annoncait 1 224 et 213)")
    print(f"    lignes sans `_lot` dans le fichier : {sans_lot} "
          f"-- il en restera {sans_lot - lig_her}")
    if doublons:
        print(f"\n  ATTENTION : {len(doublons)} identites (cicy, b, c) apparaissent "
              f"dans plusieurs taches heritees.")
        print(f"    Supprimer par identite en retirerait plus que prevu. "
              f"A traiter avant d'appliquer quoi que ce soit.")

    rapport_json = {
        'options': {k: getattr(args, k) for k in
                    ('cicy', 'tous_groupes', 'replier_orbites',
                     'controle_orbites', 'taille_lot', 'max_realisations')},
        'empreinte': empreinte,
        'ancres': {'recalculer_retrouvee': ok_r, 'intact_retrouvee': ok_i,
                   'mesure_recalculer': b['recalculer'],
                   'mesure_intact_5_34': intact_5_34},
        'bilan_couples': b,
        'couples': {k: [list(c) for c in sorted(v)] for k, v in cats.items()},
        'statuts_realisations': {str(c): rows for c, rows in statuts.items()},
        'lots_a_recalculer': [[list(i), c, tr, npm] for i, c, tr, npm in touches],
        'lots_touches_hors_cible': [list(t[0]) for t in touches
                                    if t[0] not in ids_cible
                                    and (t[0] in faits
                                         or ('T', t[0][1]) in faits)],
        'lots_deja_faits_hors_cible': [list(f) for f in hors_cible],
        'taches_heritees_a_recalculer': taches_heritees,
        'cles_heritees': [list(c) for c in sorted(cles_heritees)],
        'doublons_identite': [list(c) for c in sorted(doublons)],
        'lignes': {'explicites': lig_expl, 'survit_explicites': surv_expl,
                   'heritees': lig_her, 'survit_heritees': surv_her,
                   'sans_lot_total': sans_lot},
        'realisations': {'dans_les_tranches': vues, 'permutantes': permutantes},
    }
    chemin_sortie = os.path.join(args.output_dir, args.sortie)
    with open(chemin_sortie, 'w', encoding='utf-8') as f:
        json.dump(rapport_json, f, indent=1)
    print(f"\n  rapport ecrit : {chemin_sortie}")
    print("  (ce script n'a rien modifie d'autre)\n")


if __name__ == '__main__':
    EF._sortie_tolerante()
    principal()
