#!/usr/bin/env python3
"""
ancres_port.py -- Le balayage porte-t-il les verdicts des §5.36 a §5.38 ?

A QUOI CA REPOND
----------------
Les verdicts des trois dernieres sections vivent hors du balayage, dans des
fichiers que git protege :

    tous_indetermines.jsonl   les 2 440 candidats et leurs verdicts par lambda
    lieu_de_base_rv3.jsonl    944 lignes : f a un lieu de base SUR Y (§5.37)
    lieu_de_base_rc2.jsonl     34 lignes : idem a rank_C = 2 (§5.38)
    verdict_z4.json           #7745 et #6947, rank_C = 2 (§5.36)

Les porter dans `scan_wilson6` demande ~30 h. Ce script dit AVANT de les
depenser -- sur quelques candidats verifies en `--cicy` -- si le branchement
reproduit ces verdicts, et APRES, sur le balayage entier, s'il les a tous
reproduits. C'est la meme question posee des deux cotes du calcul.

CE QU'IL COMPARE, ET A QUELLE MAILLE
------------------------------------
`tous_indetermines.jsonl` indexe les lambda par leur RANG dans
`out['solutions']` ; le JSONL du balayage ecrit la VALEUR de lambda. Les deux
ne se recollent pas ligne a ligne sans rejouer le calcul -- ce qui reviendrait
a refaire ce qu'on verifie.

On compare donc a la maille de l'IDENTITE (cicy, b, c, groupe), par
COMPTAGES : combien de lambda survivent, combien tombent par Hoppe, combien
sortent non fibres. Un comptage juste sur chaque identite ne laisse pas de
place a une permutation de verdicts entre lambda d'une meme identite -- et
c'est la seule ambiguite que la maille introduit.

LE FACTEUR DE REALISATION, ET POURQUOI IL N'EST PAS UNE COMMODITE
------------------------------------------------------------------
La reference n'a evalue qu'UNE realisation de Gamma par couple (candidat,
groupe) ; le balayage les evalue toutes. A cette maille, le balayage porte
donc r fois plus de lignes lambda que la reference -- r = 4 sur #4078, dont
Braun donne quatre realisations de ZZ2. Le premier essai de ce script
comparait 2 a 8 et criait a l'ecart sur 35 identites sur 35, alors que les
verdicts concordaient parfaitement, realisation par realisation.

r est deduit du nombre TOTAL de lignes, une fois par identite, et impose
ensuite a chaque champ. Un r ajuste champ par champ serait un parametre
libre : il rendrait l'accord automatique et le controle vide. Sous cette
forme, la comparaison teste deux choses a la fois -- que le balayage
reproduit la reference, ET que les r realisations s'accordent entre elles.
Le rapport AFFICHE r : un facteur qu'on divise en silence transforme un
accord sous hypothese en accord.

L'ATTENDU DES DEUX STRATES
--------------------------
Les §5.37 et §5.38 ne concluent pas « la plupart » : ils concluent 944 sur
944 et 34 sur 34. Donc, sur les strates (rank_C, rang_V) = (1, 3) et (2, 3),
TOUTE ligne lambda qui passe Hoppe sans certificat de surjectivite doit
sortir `fibre = False`. C'est un attendu exact, donc falsifiable : s'il en
sort 943, l'ancre tombe.

Ce script VERIFIE D'ABORD SES PROPRES REFERENCES : il recompte 944 et 34 dans
les deux fichiers de lieu de base, et exige que chaque ligne y porte son
temoin et son F.Y. Une ancre tiree d'un fichier qu'on n'a pas relu ne vaut
pas mieux que le fichier.

Usage :
    python -u ancres_port.py                       # les references seules
    python -u ancres_port.py scan_w6_c4078         # un dossier `--cicy`
    python -u ancres_port.py scan_wilson6          # le balayage entier
"""
import os
import sys
import json
import argparse
from collections import Counter, defaultdict

STRATES_LIEU = {(1, 3), (2, 3)}


# ======================================================================
# 1. Les references, relues et recomptees
# ======================================================================

def lire_references(ici):
    """Rend (attendu, bilan). `attendu` : {identite: comptages attendus}."""
    src = os.path.join(ici, 'tous_indetermines.jsonl')
    if not os.path.exists(src):
        sys.exit(f"  ABSENT : {src}")

    attendu = {}
    bilan = Counter()
    for d in map(json.loads, open(src, encoding='utf-8')):
        if d.get('etat') != 'ok':
            bilan[f"candidat sans verdict : {d.get('etat')}"] += 1
            continue
        rc, rv = d.get('rank_C'), d.get('rang_V')
        cle = (d['cicy'], d['b'], d['c'], d['groupe'])
        att = Counter()
        for mz in d.get('mesures', []):
            if mz.get('hoppe') is True and mz.get('surjectif'):
                att['survit'] += 1
            elif mz.get('hoppe') is False:
                # ELIMINEE. Peu importe par quoi : `hoppe = False` couvre
                # ici le p = 1, qui EST h^0(V), et que le balayage attrape
                # une etape plus tot. Voir `lire_scan`.
                att['elimine'] += 1
            elif mz.get('hoppe') is True and not mz.get('surjectif'):
                # C'est ici que les §5.37 / §5.38 ont tranche -- et ils ont
                # tranche INTEGRALEMENT sur les deux strates. Ailleurs, la
                # ligne reste en attente, et c'est le reliquat du §5.38.
                att['elimine' if (rc, rv) in STRATES_LIEU
                    else 'indetermine'] += 1
            else:
                att['indetermine'] += 1
        # Compte a part : il n'entre PAS dans le total qui sert a deduire le
        # facteur de realisation, et il n'est pas exige ligne a ligne -- la
        # route peut differer d'une realisation a l'autre. Il sert au bilan
        # global, ou il est exigible.
        n_lieu = sum(1 for mz in d.get('mesures', [])
                     if mz.get('hoppe') is True and not mz.get('surjectif')
                     and (rc, rv) in STRATES_LIEU)
        attendu[cle] = {'strate': (rc, rv), 'comptes': att,
                        'n_lieu': n_lieu, 'n_lambda': d.get('n_lambda')}
        for k, v in att.items():
            bilan[f"attendu : {k}"] += v
        bilan['attendu : dont non fibrees (lieu de base)'] += n_lieu
    return attendu, bilan


def verifier_lieux(ici):
    """Recompte les deux fichiers de lieu de base, et exige leurs gardes.

    Sans cette passe, les ancres 944 et 34 seraient recopiees du document,
    pas lues -- exactement le defaut trouve au §5.38 en verifiant l'etat de
    la suite : un chiffre ecrit une fois, puis reconduit.
    """
    out = {}
    for nom, champ_temoin in (('lieu_de_base_rv3.jsonl', 'temoin_verifie'),
                              ('lieu_de_base_rc2.jsonl', 'temoin_rang_1')):
        chemin = os.path.join(ici, nom)
        if not os.path.exists(chemin):
            out[nom] = {'absent': True}
            continue
        n = ok = sans_temoin = sans_fy = 0
        completion = 0
        for d in map(json.loads, open(chemin, encoding='utf-8')):
            n += 1
            if d.get('etat') != 'ok':
                continue
            if not d.get(champ_temoin):
                sans_temoin += 1
                continue
            fy = d.get('F_Y')
            if fy is None or fy <= 0:
                sans_fy += 1
                continue
            if isinstance(d.get('dim_F'), int) and isinstance(d.get('K'), int) \
                    and d['dim_F'] > d['K']:
                completion += 1
            ok += 1
        out[nom] = {'lignes': n, 'eliminations': ok,
                    'sans temoin': sans_temoin, 'sans F.Y > 0': sans_fy,
                    'par completion H^(dimF-K)': completion}
    return out


# ======================================================================
# 2. Le balayage, lu a la meme maille
# ======================================================================

def lire_scan(chemin):
    """{identite: comptages observes}, plus le reliquat par strate.

    DEUX NIVEAUX, ET POURQUOI
    -------------------------
    Le NIVEAU 1 est l'ISSUE : survit / elimine / indetermine. C'est ce qui
    doit coincider, et c'est le seul niveau bloquant. Les deux sens du
    §5.34 y sont : un SURVIT usurpe et un elimine a tort y tombent tous
    les deux.

    Le NIVEAU 2 est la ROUTE : par Hoppe, par h^0 equivariant, par le lieu
    de base. Il est DECLARE et non exige, parce que deux routes peuvent
    mener a la meme elimination sans qu'aucune des deux soit fausse :

      - `echantillon_rank_c2.py` appelle Hoppe SANS la garde h^0 qui le
        precede dans le balayage. Un h^0(V) non nul y ressort en
        `hoppe = False, valeurs {1: 1}` -- car le p = 1 de Hoppe EST
        h^0(V) -- la ou le balayage tue la ligne plus tot et laisse
        `hoppe = None`. Meme fait, meme verdict, deux champs ;
      - deux realisations de Gamma ne portent pas le meme sigma (§5.35).
        Sur #4078 a rank_C = 2, deux realisations donnent un espace
        equivariant de dimension 21 et tombent par le lieu de base, deux
        autres de dimension 24 et tombent par Hoppe. Toutes quatre sont
        eliminees ; exiger la meme route reviendrait a exiger que sigma
        n'existe pas.

    Exiger le niveau 2 ferait donc crier l'ancre sur deux accords. Ne PAS
    exiger le niveau 1 la rendrait aveugle. On separe, et on affiche les
    deux -- un desaccord de route est une mesure, pas un silence.
    """
    obs = defaultdict(Counter)
    routes = defaultdict(Counter)
    lots = defaultdict(set)
    par_lot = defaultdict(lambda: defaultdict(Counter))
    reliquat = Counter()
    n = n_controle = 0
    with open(chemin, encoding='utf-8') as f:
        for ligne in f:
            if not ligne.strip():
                continue
            d = json.loads(ligne)
            if d.get('etat') != 'ok':
                continue
            lot = d.get('_lot')
            # --- LES LIGNES DE CONTROLE D'ORBITE SORTENT D'ICI ----------
            # Un lot ('C', k, j, t) reevalue POUR DE VRAI un membre non
            # representant, afin de comparer son verdict a celui qui lui a
            # ete recopie (§5.25). Ses lignes sont donc, PAR CONSTRUCTION,
            # des doublons deliberes : la meme identite y recoit une
            # seconde fois ses verdicts.
            #
            # La reference n'a pas d'equivalent -- `echantillon_rank_c2.py`
            # ne replie rien et ne controle rien. Les compter reviendrait a
            # comparer un fichier a une reference qui n'a jamais eu ces
            # lignes : l'identite paraitrait porter plus de lambda qu'elle
            # n'en a, et le facteur r se decalerait avec.
            #
            # Ce n'est PAS une exemption de controle : ces lignes sont
            # verifiees par le balayage lui-meme, qui compare chaque membre
            # reevalue a son representant et declare le run invalide sur la
            # moindre discordance. Elles sont deja contrelees, ailleurs, et
            # mieux qu'ici. On les COMPTE et on l'affiche.
            if lot and lot[0] == 'C':
                n_controle += 1
                continue
            n += 1
            cle = (d.get('cicy'), json.dumps(d.get('b_charges')),
                   json.dumps(d.get('c_charges')), d.get('groupe'))
            c = obs[cle]
            # Les LOTS qui ont ecrit cette identite. C'est la mesure de
            # completude la plus directe qu'on ait, et elle ne demande
            # aucune reference exterieure : le nombre de tranches ne depend
            # que de (cicy, groupe), donc deux identites d'un meme couple
            # doivent en porter autant l'une que l'autre.
            if lot:
                lots[cle].add(tuple(lot))
            rc = len(d.get('c_charges') or [])
            strate = (rc, d.get('rang_V'))
            cl_lot = tuple(lot) if lot else None
            if d.get('survit'):
                c['survit'] += 1
                par_lot[cle][cl_lot]['survit'] += 1
                routes[cle]['survit'] += 1
            elif d.get('indetermine'):
                c['indetermine'] += 1
                par_lot[cle][cl_lot]['indetermine'] += 1
                routes[cle]['indetermine'] += 1
                # LE RELIQUAT. Sur les deux strates traitees il doit etre
                # NUL : c'est l'enonce meme des §5.37 et §5.38, et il ne
                # depend d'aucune route.
                reliquat[strate] += 1
            else:
                c['elimine'] += 1
                par_lot[cle][cl_lot]['elimine'] += 1
                if d.get('fibre') is False:
                    routes[cle]['par le lieu de base'] += 1
                elif d.get('hoppe_complet') is False:
                    routes[cle]['par Hoppe'] += 1
                elif d.get('h0_equivariant'):
                    routes[cle]['par h0 equivariant'] += 1
                else:
                    routes[cle]['autre'] += 1
    return obs, routes, lots, par_lot, reliquat, n, n_controle


def principal():
    ap = argparse.ArgumentParser()
    ap.add_argument('cible', nargs='?', default=None,
                    help='dossier de scan, ou chemin du JSONL. Sans lui, on '
                         'ne verifie que les references.')
    ap.add_argument('--strate', default=None,
                    help='ne comparer qu une strate, ex. "1,3"')
    args = ap.parse_args()

    ici = os.path.dirname(os.path.abspath(__file__))
    attendu, bilan = lire_references(ici)
    lieux = verifier_lieux(ici)

    print(f"\n{'=' * 78}\n  LES REFERENCES\n{'=' * 78}")
    print(f"    {len(attendu)} identites (cicy, b, c, groupe) dans "
          f"tous_indetermines.jsonl")
    for k, v in sorted(bilan.items()):
        print(f"    {v:>7}  {k}")
    print()
    for nom, d in lieux.items():
        if d.get('absent'):
            print(f"    {nom} : ABSENT -- l'ancre correspondante ne sera pas "
                  f"verifiee")
            continue
        print(f"    {nom} : " + ", ".join(f"{k} = {v}" for k, v in d.items()))
    # Les deux gardes des §5.37 / §5.38, exigees et non supposees.
    for nom, cible in (('lieu_de_base_rv3.jsonl', 944),
                       ('lieu_de_base_rc2.jsonl', 34)):
        d = lieux.get(nom, {})
        if d.get('absent'):
            continue
        if d['eliminations'] != cible or d['lignes'] != cible:
            print(f"\n  /!\\ {nom} : {d['eliminations']} eliminations sur "
                  f"{d['lignes']} lignes, attendu {cible} sur {cible}.")
            print(f"      L'ancre elle-meme a bouge. NE PAS s'en servir "
                  f"avant d'avoir compris pourquoi.")
            return 2

    attendu_non_fibre = sum(v['n_lieu'] for v in attendu.values())
    attendu_survit = sum(v['comptes']['survit'] for v in attendu.values())
    print(f"\n    Attendu du balayage sur ces identites : "
          f"{attendu_survit} SURVIT, {attendu_non_fibre} non fibrees.")
    if attendu_non_fibre != 944 + 34:
        print(f"    /!\\ 944 + 34 = 978 attendu, {attendu_non_fibre} deduit "
              f"de tous_indetermines.jsonl. Les deux sources divergent.")
        return 2

    if args.cible is None:
        print(f"\n  Rien a comparer (aucun dossier passe). Les references "
              f"sont coherentes.\n")
        return 0

    chemin = args.cible
    if os.path.isdir(chemin):
        chemin = os.path.join(chemin, 'results_equivariance_f.jsonl')
    if not os.path.exists(chemin):
        sys.exit(f"  ABSENT : {chemin}")

    obs, routes, lots, par_lot, reliquat, n_lignes, n_controle = lire_scan(chemin)
    filtre = None
    if args.strate:
        filtre = tuple(int(x) for x in args.strate.split(','))

    print(f"\n{'=' * 78}\n  COMPARAISON : {chemin}\n{'=' * 78}")
    print(f"    {n_lignes} lignes `ok`, {len(obs)} identites")
    print(f"    {n_controle} lignes de CONTROLE D'ORBITE ecartees de la "
          f"comparaison (lots 'C') --")
    print(f"      doublons deliberes du §5.25, verifies par le balayage "
          f"lui-meme, absents de la reference\n")

    vus = ecarts = absents = partielles = 0
    detail = Counter()
    detail_partiel = Counter()
    facteurs = Counter()
    lignes_ecart = []

    # ------------------------------------------------------------------
    # LE NIVEAU 1, REECRIT -- son hypothese a ete DEMENTIE par la mesure
    # ------------------------------------------------------------------
    # Il exigeait `observe == attendu x r` : que les r realisations de Gamma
    # rendent toutes le verdict de celle qu'a vue la reference. Sur #480 /
    # ZZ2 x ZZ2, `diag_ecart.py` a mesure le contraire, et sans ambiguite :
    #
    #     23 tranches, TROIS profils distincts
    #        12 tranches -> 16 SURVIT, 48 indetermine
    #         6 tranches -> 28 SURVIT, 36 indetermine
    #         5 tranches -> 64 SURVIT
    #     et DEUX dimensions d'espace equivariant pour la meme identite :
    #        17/67 sur 1104 lignes, 16/67 sur 368
    #
    # Deux realisations, deux quotients X/Gamma, deux verdicts legitimes --
    # c'est le §5.35, cette fois sur des SURVIVANTS. Exiger l'egalite,
    # c'etait exiger que sigma n'existe pas.
    #
    # CE QU'ON PEUT ENCORE EXIGER, ET QUI NE SUPPOSE RIEN
    # ---------------------------------------------------
    # La realisation qu'a vue la reference est l'une des r du balayage. Elle
    # y contribue donc, a elle seule, ses propres comptes. D'ou :
    #
    #     observe[champ] >= attendu[champ],  pour chaque champ
    #
    # Necessaire, non suffisant -- mais falsifiable et sans hypothese : si
    # le branchement cessait d'agir, `elimine` chuterait sous l'attendu et
    # l'ancre tomberait. C'est le sens qui compte, celui du §5.34 : un
    # verdict qui DISPARAIT.
    #
    # L'ecart a `attendu x r` reste MESURE et affiche : c'est la divergence
    # entre realisations, un fait a consigner, pas un defaut.
    #
    # LA COMPLETUDE N'EST PLUS INFEREE ICI
    # ------------------------------------
    # Elle l'etait a partir du nombre de lots par identite, en supposant
    # qu'il ne depend que de (cicy, groupe). Faux : `idx` est filtre par
    # `groupes_utiles`, qui est une propriete du CANDIDAT. Deux candidats
    # d'une meme CICY et d'un meme groupe peuvent donc porter 22 et 23
    # tranches sans que rien ne cloche -- mesure sur #480, ou 44 identites
    # en portent 22 et 2 en portent 23.
    #
    # La completude est etablie ailleurs, et mieux : par le balayage
    # (« Lots : N termines sur N ») et par `retirer_lots.py --verifier`.
    # Ce script ne la mesure plus et le DIT, plutot que de la mesurer mal.
    divergences = Counter()
    n_divergentes = 0

    for cle, att in attendu.items():
        if filtre and att['strate'] != filtre:
            continue
        o = obs.get(cle)
        if o is None:
            absents += 1
            continue
        vus += 1
        a = att['comptes']

        # Le facteur de realisation : la reference n'a evalue qu'UNE
        # realisation de Gamma par couple, le balayage les evalue toutes.
        # r est deduit du nombre TOTAL de lignes, une fois par identite, et
        # il sert a MESURER la divergence -- plus a l'exiger.
        a_total = sum(a.values())
        o_total = sum(o.values())
        if a_total == 0:
            continue
        r = o_total / a_total if a_total else 0
        if float(r).is_integer():
            facteurs[int(r)] += 1

        # NIVEAU 1, BLOQUANT : la realisation de la reference est l'une de
        # celles du balayage ; elle y contribue ses propres comptes. Donc
        # aucun champ ne peut y etre INFERIEUR. Un verdict qui disparait est
        # le sens qui compte (§5.34).
        for champ in ('survit', 'elimine', 'indetermine'):
            if o.get(champ, 0) < a.get(champ, 0):
                ecarts += 1
                detail[f"{champ} : la reference en porte {a.get(champ, 0)}, "
                       f"le balayage {o.get(champ, 0)} -- un verdict a "
                       f"DISPARU"] += 1
                lignes_ecart.append((cle, att['strate'], champ,
                                     a.get(champ, 0), o.get(champ, 0), None))
                break
        else:
            # DECLARE : les realisations divergent-elles, et sur combien de
            # profils distincts ?
            if any(a.get(ch, 0) * r != o.get(ch, 0)
                   for ch in ('survit', 'elimine', 'indetermine')):
                n_divergentes += 1
                profils = {tuple(sorted(p.items()))
                           for p in par_lot.get(cle, {}).values()}
                divergences[f"#{cle[0]} {cle[3]} strate {att['strate']} : "
                            f"{len(profils)} profils distincts sur "
                            f"{o_total} lignes"] += 1

    # REGLE DES FILTRES : les deux cotes. Un rapport qui ne compterait que
    # les accords serait un plaidoyer ; un rapport qui ne compterait que les
    # ecarts ne dirait pas sur quoi ils portent. Et le facteur r se DECLARE :
    # une division silencieuse ferait passer pour un accord ce qui est un
    # accord SOUS HYPOTHESE -- que les r realisations rendent le meme
    # verdict, ce qui est justement l'hypothese que ce controle met a
    # l'epreuve champ par champ.
    print(f"    {vus} identites comparees, {vus - ecarts} d'accord, "
          f"{ecarts} en ECART")
    print(f"    {absents} identites de la reference ABSENTES du balayage")
    print(f"      -- la COMPLETUDE n'est pas mesuree ici : voir « Lots : N "
          f"termines sur N » du balayage")
    print(f"         et `retirer_lots.py --verifier`. Elle l'etait, mal, a "
          f"partir du nombre de lots par identite.")
    if facteurs:
        print(f"    realisations par identite (r) : "
              + ", ".join(f"r={k} sur {v} identites"
                          for k, v in sorted(facteurs.items())))
        if set(facteurs) == {1}:
            print(f"      (r = 1 partout : la reference et le balayage "
                  f"portent le meme nombre de realisations)")
    # Les partielles se NOMMENT. « 66 partielles » sans dire lesquelles ni
    # de combien, c'est un compte qu'on ne peut pas verifier -- et une fois
    # le balayage termine, il doit tomber a zero : c'est la seule facon de
    # voir qu'une identite ne se remplit plus.
    for k, v in detail.most_common(15):
        print(f"      {v:>5}  {k}")
    for cle, strate, champ, a, o, r in lignes_ecart[:15]:
        print(f"      #{cle[0]} {cle[3]} strate {strate} : {champ} "
              f"attendu {a}" + (f" x {r}" if r else "") + f", observe {o}")

    # ---- NIVEAU 2 : la route, declaree et non exigee -------------------
    par_route = Counter()
    for cle in obs:
        for k, v in routes[cle].items():
            par_route[k] += v
    print(f"\n    Realisations DIVERGENTES (mesure du §5.35, non un defaut) : "
          f"{n_divergentes} identites sur {vus}")
    for k, v in divergences.most_common(10):
        print(f"      {k}")
    if len(divergences) > 10:
        print(f"      ... et {len(divergences) - 10} autres")
    print(f"\n    Par quelle route (niveau 2, DECLARE, non exige) :")
    for k, v in par_route.most_common():
        print(f"      {v:>7}  {k}")
    print(f"      Deux routes peuvent mener a la meme elimination sans "
          f"qu'aucune soit fausse :")
    print(f"      la garde h0 precede Hoppe dans le balayage et pas dans la "
          f"reference, et")
    print(f"      deux realisations ne portent pas le meme sigma (§5.35). "
          f"C'est mesure, pas exige.")

    # ---- L'ANCRE QUI NE DEPEND D'AUCUNE ROUTE --------------------------
    # C'est l'enonce meme des §5.37 et §5.38 : sur les deux strates
    # traitees, plus AUCUNE ligne en attente. Il ne se contourne pas en
    # changeant de route, il ne depend pas du facteur r, et il tombe des
    # que le branchement cesse d'agir -- les lignes redeviennent
    # `indetermine`. C'est le controle le plus dur de ce script.
    print(f"\n    Reliquat indetermine par strate (BLOQUANT sur les deux "
          f"strates traitees) :")
    if not reliquat:
        print(f"      aucune ligne indeterminee dans ce lot")
    for strate, v in sorted(reliquat.items()):
        marque = '  <- DOIT ETRE NUL' if strate in STRATES_LIEU else \
                 '  (strate sans module -- le reliquat du §5.38)'
        print(f"      {v:>7}  rank_C={strate[0]}, rang_V={strate[1]}{marque}")
    fuite = sum(v for s, v in reliquat.items() if s in STRATES_LIEU)

    if fuite:
        print(f"\n  ARRET : {fuite} lignes restent indeterminees sur les "
              f"strates que les §5.37 / §5.38 declarent tranchees.")
        print(f"  Le branchement n'agit pas, ou n'agit pas partout.\n")
        return 1
    if ecarts:
        print(f"\n  ARRET : {ecarts} identites PERDENT un verdict que la "
              f"reference porte (niveau 1).\n")
        return 1
    if not vus:
        print(f"\n  /!\\ AUCUNE identite comparee. Ce n'est pas un accord : "
              f"c'est un controle vide.\n")
        return 1
    print(f"\n  Les {vus} identites comparees reproduisent l'issue des "
          f"§5.36 a §5.38, et le reliquat est nul sur les deux strates "
          f"traitees.\n")
    return 0


if __name__ == '__main__':
    sys.exit(principal())
