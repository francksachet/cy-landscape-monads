#!/usr/bin/env python3
"""
retirer_lots.py -- Retire du checkpoint et du JSONL les lots faussees par le
defaut du 5.34, pour que la relance les recalcule.

Prend en entree le rapport de `portee_substitution.py` et n'invente rien : la
liste des lots a retirer vient de lui, avec ses gardes (empreinte du
checkpoint, ancres du 5.34).

CE QU'IL FAIT, ET DANS QUEL ORDRE
---------------------------------
1. Retire du checkpoint les identifiants de lot touches.
2. Retire les entrees HERITEES ('T', k) des taches touchees. Ces entrees-la
   viennent d'un checkpoint sequentiel migre et couvrent la tache ENTIERE.
3. Reecrit le JSONL en ne gardant que
     - les lignes dont le `_lot` est encore dans le checkpoint,
     - les lignes SANS `_lot` dont l'identite (cicy, b, c) n'appartient a
       aucune tache heritee retiree.
4. Recalcule les trois compteurs sur les lignes conservees.

POURQUOI L'ETAPE 3 NE SE DEDUIT PAS DE LA MACHINERIE DE REPRISE
---------------------------------------------------------------
`Progression.charger` filtre deja le JSONL sur les lots valides. Mais sa
passe 2 dit :

    if lot is None or tuple(lot) in prog.faits:   -> on garde

Une ligne SANS `_lot` est donc gardee INCONDITIONNELLEMENT -- l'ancien format
sequentiel ecrivait dans l'ordre et etait tronque a son offset, ce qui
suffisait alors. Retirer une entree heritee ('T', k) sans supprimer ses
lignes soi-meme relancerait donc le calcul de la tache tout en gardant
l'ancien resultat : des DOUBLONS, dont l'un est faux, et rien pour les
signaler. Il y a 55 170 lignes sans `_lot` dans `scan_wilson4`.

POURQUOI L'ETAPE 4 NON PLUS
---------------------------
`charger` relit `compteurs` tel quel et ne les recalcule jamais apres avoir
filtre le fichier. Des lots retires sans retoucher les compteurs laisseraient
des totaux qui comptent des lignes qui n'existent plus.

La regle exacte est celle de `_ecrire`, et elle se verifie : sur le fichier
intact, elle doit rendre les compteurs du checkpoint AU CHIFFRE PRES. Si elle
ne les rend pas, le script s'arrete -- c'est la seule reference exterieure
disponible sur ce point.

    survivants   = lignes avec survit
    indetermines = lignes avec etat == 'ok' ET indetermine
    ecartes      = lignes avec etat != 'ok'

Usage :
    python -u retirer_lots.py scan_wilson4                  # a blanc
    python -u retirer_lots.py scan_wilson4 --appliquer      # ecrit
"""
import os
import sys
import json
import argparse
from collections import Counter


HORS_DOMAINE = 'hors domaine (modele S/I non valide)'


def contradictions(chemin):
    """
    Identites (cicy, b, c) portant A LA FOIS « hors domaine » et « ok ».

    « Hors domaine » dit que le modele S/I ne represente pas H^0 pour ces
    charges ; « ok » dit qu'on a calcule des verdicts dessus. Les deux sur la
    meme identite prouvent, SANS AUCUN RECALCUL et sans reference exterieure,
    que deux versions du code ont ecrit dans ce fichier.

    C'est le seul controle de ce dossier qui ne coute rien : une lecture.
    """
    from collections import defaultdict
    vus = defaultdict(set)
    with open(chemin, encoding='utf-8') as f:
        for ligne in f:
            if not ligne.strip():
                continue
            try:
                d = json.loads(ligne)
            except json.JSONDecodeError:
                continue
            e = d.get('etat')
            if e not in (HORS_DOMAINE, 'ok'):
                continue
            cle = (d.get('cicy'), json.dumps(d.get('b_charges')),
                   json.dumps(d.get('c_charges')))
            vus[cle].add('hors' if e == HORS_DOMAINE else 'ok')
    return sorted(k for k, v in vus.items() if len(v) == 2)


def compter(chemin, garder=None):
    """
    Compte les lignes selon la regle de `_ecrire`. `garder(d)` filtre.
    Rend (compteurs, total, survit_total).
    """
    c = {'survivants': 0, 'indetermines': 0, 'ecartes': 0}
    n = 0
    with open(chemin, encoding='utf-8') as f:
        for ligne in f:
            if not ligne.strip():
                continue
            d = json.loads(ligne)
            if garder is not None and not garder(d):
                continue
            n += 1
            if d.get('etat') != 'ok':
                c['ecartes'] += 1
            else:
                c['indetermines'] += bool(d.get('indetermine'))
            c['survivants'] += bool(d.get('survit'))
    return c, n


def principal():
    ap = argparse.ArgumentParser()
    ap.add_argument('output_dir')
    ap.add_argument('--rapport', default='portee_5_34.json')
    ap.add_argument('--sauvegarde', default=None,
                    help="dossier ou les deux fichiers d'origine ont ete "
                         "copies ; defaut <output_dir>_avant_5.35")
    ap.add_argument('--appliquer', action='store_true',
                    help="sans cette option, rien n'est ecrit")
    ap.add_argument('--verifier', action='store_true',
                    help="ne verifie QUE l'accord entre les compteurs du "
                         "checkpoint et le fichier, et s'arrete. A passer "
                         "apres chaque session de calcul.")
    ap.add_argument('--refaire-controles', action='store_true',
                    help="retirer les lots de controle d'orbite pour qu'ils "
                         "se rejouent contre un JSONL complet. A passer une "
                         "fois tous les lots 'T' termines, si le balayage a "
                         "ete fractionne en plusieurs sessions.")
    ap.add_argument('--refaire', action='store_true',
                    help="reappliquer alors que le marqueur dit que c'est "
                         "deja fait (a n'utiliser qu'apres restauration de "
                         "la sauvegarde)")
    args = ap.parse_args()

    dst = os.path.join(args.output_dir, 'results_equivariance_f.jsonl')
    ck_p = os.path.join(args.output_dir, 'progress_equivariance_f.json')
    rap_p = os.path.join(args.output_dir, args.rapport)
    # `--verifier` ne lit que le fichier et le checkpoint : il doit rester
    # utilisable APRES la reparation, quand le rapport n'a plus lieu d'etre.
    for p in (dst, ck_p):
        if not os.path.exists(p):
            sys.exit(f"  ABSENT : {p}")
    # Seul le retrait des lots fausses par le 5.34 a besoin du rapport de
    # `portee_substitution`. `--verifier` et `--refaire-controles` ne lisent
    # que le fichier et le checkpoint, et doivent rester utilisables sur un
    # dossier qui n'a jamais eu de reparation a faire -- `scan_wilson5`,
    # par exemple, calcule d'un bloc dans un seul etat du code.
    if not (args.verifier or args.refaire_controles) \
            and not os.path.exists(rap_p):
        sys.exit(f"  ABSENT : {rap_p}")

    with open(ck_p, encoding='utf-8') as f:
        ck = json.load(f)

    # ---- mode verification seule --------------------------------------
    # A passer apres chaque session de calcul. `equivariance_f` ecrit ses
    # lignes (flush + fsync) AVANT d'enregistrer le lot, donc une coupure
    # ne peut que laisser des lignes sans lot -- que la reprise jette. Mais
    # si un lot restait enregistre avec des lignes manquantes, la reprise
    # le retirerait de `faits` SANS decrementer les compteurs, qui
    # compteraient alors des lignes qui n'existent plus. Ce mode le voit.
    if args.verifier:
        c, n = compter(dst)
        ref = {k: int(ck['compteurs'].get(k, -1)) for k in c}
        print(f"\n  {n} lignes dans le JSONL, {len(ck['lots'])} lots au "
              f"checkpoint")
        print(f"    checkpoint : {ref}")
        print(f"    recompte   : {c}")
        # Deux controles de plus, qui ne coutent qu'une lecture.
        try:
            from empreinte_code import empreinte_code, repartition
            emp, nf = empreinte_code(os.path.dirname(os.path.abspath(__file__)))
            rep = repartition(dst)
            print(f"\n  Versions du code presentes ({nf} fichiers surveilles, "
                  f"code actuel {emp}) :")
            for k, v in sorted(rep.items(), key=lambda x: -x[1]):
                marque = '  <- code actuel' if k == emp else ''
                print(f"    {v:>8} lignes   "
                      f"{k if k else '(avant le marquage)'}{marque}")
            if len(rep) > 1:
                print("    -> le fichier MELANGE plusieurs versions du code.")
        except ImportError:
            print("\n  (empreinte_code absent : versions non verifiees)")
        contra = contradictions(dst)
        print(f"\n  Identites portant a la fois « hors domaine » et « ok » : "
              f"{len(contra)}")
        if contra:
            print(f"    sur {len({k[0] for k in contra})} CICYs, dont "
                  f"{sorted({k[0] for k in contra})[:6]}")
            print("    -> contradiction interne : deux versions du code ont "
                  "ecrit ici.")
        if c == ref and not contra:
            print("\n    -> ACCORD. Rien a signaler.\n")
            return
        if c == ref:
            print("\n    -> compteurs d'accord, mais le fichier n'est pas "
                  "homogene.\n")
            sys.exit(1)
        print("    -> DESACCORD. Les compteurs du checkpoint ne decrivent "
              "plus le fichier.")
        print("       Les verdicts eux-memes restent bons (chaque ligne "
              "porte le sien) ;")
        print("       ce sont les TOTAUX affiches en fin de run qui sont "
              "faux.\n")
        sys.exit(1)

    # ---- refaire les controles d'orbite --------------------------------
    # Le controle du repli (5.25) compare, en fin de session, les lignes du
    # membre de controle a celles du representant RELUES DANS LE JSONL. Si
    # le lot de controle a tourne dans une session et le representant dans
    # une autre, la comparaison lit une liste vide -- et `if a and a != b_`
    # la laisse passer EN SILENCE. Fractionner un run affaiblit donc le
    # controle sans le dire, ce qui est exactement le defaut du 5.25.
    #
    # Remede : une fois tous les lots 'T' termines, retirer les lots 'C' et
    # relancer une derniere fois. Ils se recalculent alors contre un JSONL
    # complet, et le verdict du controle vaut enfin quelque chose.
    if args.refaire_controles:
        faits = {tuple(x[0]): x[1] for x in ck['lots']}
        ctrl = {f for f in faits if f[0] == 'C'}
        print(f"\n  {len(ctrl)} lots de controle sur {len(faits)} au total")
        if not ctrl:
            print("  Rien a refaire.\n")
            return
        restants = set(faits) - ctrl
        apres, n_apres = compter(dst, lambda d: d.get('_lot') is None
                                 or tuple(d['_lot']) in restants)
        avant, n_avant = compter(dst)
        print(f"  lignes : {n_avant} -> {n_apres} ({n_avant - n_apres} retirees)")
        if not args.appliquer:
            print("  Rien n'a ete ecrit. --appliquer pour executer.\n")
            return
        tmp = dst + '.nouveau'
        with open(dst, encoding='utf-8') as s, open(tmp, 'w', encoding='utf-8') as o:
            for ligne in s:
                if not ligne.strip():
                    continue
                d = json.loads(ligne)
                if d.get('_lot') is None or tuple(d['_lot']) in restants:
                    o.write(ligne)
        ck['lots'] = [[list(k), faits[k]] for k in sorted(restants)]
        ck['compteurs'] = apres
        tmp2 = ck_p + '.tmp'
        with open(tmp2, 'w', encoding='utf-8') as f:
            json.dump(ck, f)
        os.replace(tmp, dst)
        os.replace(tmp2, ck_p)
        print(f"  ECRIT. Relancer le balayage : les {len(ctrl)} controles se "
              f"rejoueront contre un JSONL complet.\n")
        return

    marqueur = os.path.join(args.output_dir, 'retrait_5_35_applique.json')
    if os.path.exists(marqueur) and not args.refaire:
        with open(marqueur, encoding='utf-8') as f:
            m = json.load(f)
        sys.exit(f"\n  DEJA APPLIQUE ({m.get('lignes_jetees')} lignes jetees, "
                 f"{m.get('lots_retires')} lots retires).\n"
                 f"  Rejouer ce script APRES la relance reproposerait de "
                 f"retirer les lots\n"
                 f"  qui viennent d'etre recalcules -- c'est-a-dire de jeter "
                 f"le travail.\n"
                 f"  --verifier pour un controle d'integrite, --refaire pour "
                 f"forcer.")

    with open(rap_p, encoding='utf-8') as f:
        rap = json.load(f)

    # ---- gardes -------------------------------------------------------
    if rap.get('empreinte') != ck.get('empreinte'):
        sys.exit("\n  ARRET. Le rapport et le checkpoint n'ont pas la meme\n"
                 "  empreinte : le rapport a ete produit sur un autre etat.\n"
                 "  Relancer portee_substitution.py.")
    if not (rap['ancres']['recalculer_retrouvee']
            and rap['ancres']['intact_retrouvee']):
        sys.exit("\n  ARRET. Le rapport dit que les ancres du 5.34 ne sont pas\n"
                 "  reproduites. On ne retire rien sur cette base.")
    if rap.get('doublons_identite'):
        sys.exit(f"\n  ARRET. {len(rap['doublons_identite'])} identites "
                 f"(cicy, b, c) apparaissent dans plusieurs taches heritees.\n"
                 f"  Supprimer par identite en retirerait plus que prevu.")

    sauve = args.sauvegarde or (args.output_dir.rstrip('/\\') + '_avant_5.35')
    manquants = [os.path.basename(p) for p in (dst, ck_p)
                 if not os.path.exists(os.path.join(sauve, os.path.basename(p)))]
    if manquants and args.appliquer:
        sys.exit(f"\n  ARRET. Sauvegarde incomplete dans {sauve} : {manquants}.\n"
                 f"  Cette operation n'est pas reversible ; ces fichiers ne\n"
                 f"  sont pas dans Git. Copier les deux avant d'appliquer,\n"
                 f"  ou passer --sauvegarde.")

    # ---- ce qu'on retire ----------------------------------------------
    faits = {tuple(x[0]): x[1] for x in ck['lots']}
    touches = {tuple(t[0]) for t in rap['lots_a_recalculer']}
    a_retirer = {f for f in faits if f in touches}
    her_retirer = {('T', k) for k in rap['taches_heritees_a_recalculer']}
    her_retirer &= set(faits)
    cles = {tuple(c) for c in rap['cles_heritees']}

    print(f"\n{'=' * 74}")
    print(f"  RETRAIT DES LOTS FAUSSES PAR LE 5.34 -- "
          f"{'APPLIQUE' if args.appliquer else 'A BLANC'}")
    print(f"{'=' * 74}")
    print(f"  checkpoint : {len(faits)} lots")
    print(f"    a retirer, identifiant explicite : {len(a_retirer)}")
    print(f"    a retirer, entrees heritees      : {len(her_retirer)}")
    print(f"    conserves                        : "
          f"{len(faits) - len(a_retirer) - len(her_retirer)}")

    # ---- garde n.3 : la regle des compteurs se verifie sur l'intact ----
    print(f"\n  verification de la regle des compteurs sur le fichier intact ...")
    avant, n_avant = compter(dst)
    ref = {k: int(ck['compteurs'].get(k, -1)) for k in avant}
    print(f"    checkpoint : {ref}")
    print(f"    recompte   : {avant}   ({n_avant} lignes)")
    if avant != ref:
        sys.exit("\n  ARRET. La regle de comptage ne reproduit pas les compteurs\n"
                 "  du checkpoint. Recalculer les compteurs sur cette base\n"
                 "  ecrirait des totaux faux.")
    print("    -> regle validee.")

    # ---- reecriture ----------------------------------------------------
    restants = set(faits) - a_retirer - her_retirer

    def garder(d):
        lot = d.get('_lot')
        if lot is None:
            cle = (d.get('cicy'), json.dumps(d.get('b_charges')),
                   json.dumps(d.get('c_charges')))
            return cle not in cles
        return tuple(lot) in restants

    apres, n_apres = compter(dst, garder)
    jetees = n_avant - n_apres
    print(f"\n  LIGNES")
    print(f"    avant   : {n_avant}")
    print(f"    jetees  : {jetees}   (dont "
          f"{avant['survivants'] - apres['survivants']} SURVIT, "
          f"{avant['ecartes'] - apres['ecartes']} ecartees)")
    print(f"    apres   : {n_apres}")
    print(f"\n  COMPTEURS")
    for k in ('survivants', 'indetermines', 'ecartes'):
        print(f"    {k:<13} {avant[k]:>7}  ->  {apres[k]:>7}")

    attendu = rap['lignes']['explicites'] + rap['lignes']['heritees']
    print(f"\n  attendu par le rapport : {attendu} lignes jetees, "
          f"{rap['lignes']['survit_explicites'] + rap['lignes']['survit_heritees']} SURVIT")
    if jetees != attendu:
        print(f"  ECART de {jetees - attendu} lignes avec le rapport.")
        if args.appliquer:
            sys.exit("  ARRET : on n'applique pas une chirurgie dont le compte "
                     "ne tombe pas juste.")

    if not args.appliquer:
        print(f"\n  Rien n'a ete ecrit. --appliquer pour executer.\n")
        return

    # LES DEUX FICHIERS SONT PREPARES AVANT LE MOINDRE REMPLACEMENT.
    # Il reste une fenetre de quelques millisecondes entre les deux
    # `os.replace`. Une coupure pile dedans laisse un etat incoherent que
    # la reprise ne rattrape PAS toute seule : les entrees heritees
    # ('T', k) sont conservees sans condition par `charger`, donc un
    # checkpoint qui les garde alors que leurs lignes sont parties ferait
    # disparaitre 198 taches en silence. D'ou la sauvegarde exigee plus
    # haut, et le message ci-dessous.
    tmp_dst = dst + '.nouveau'
    with open(dst, encoding='utf-8') as src_f, \
            open(tmp_dst, 'w', encoding='utf-8') as out_f:
        for ligne in src_f:
            if not ligne.strip():
                continue
            if garder(json.loads(ligne)):
                out_f.write(ligne)
    # meme forme et meme tri que `Progression.sauver`
    ck['lots'] = [[list(k), faits[k]] for k in sorted(restants)]
    ck['compteurs'] = apres
    tmp_ck = ck_p + '.tmp'
    with open(tmp_ck, 'w', encoding='utf-8') as f:
        json.dump(ck, f)

    print(f"\n  remplacement des deux fichiers (ne pas interrompre ici) ...")
    try:
        os.replace(tmp_dst, dst)
        os.replace(tmp_ck, ck_p)
    except BaseException:
        print(f"\n  /!\\ COUPURE PENDANT LE REMPLACEMENT. L'etat des deux "
              f"fichiers est indetermine.")
        print(f"      Restaurer les deux depuis {sauve}, puis relancer avec "
              f"--refaire.")
        raise

    with open(marqueur, 'w', encoding='utf-8') as f:
        json.dump({'empreinte': ck['empreinte'],
                   'lignes_jetees': jetees,
                   'survit_jetes': avant['survivants'] - apres['survivants'],
                   'lots_retires': len(a_retirer) + len(her_retirer),
                   'lignes_apres': n_apres,
                   'compteurs_apres': apres}, f, indent=1)

    print(f"\n  ECRIT. {len(restants)} lots conserves, {n_apres} lignes.")
    print(f"  Marqueur : {marqueur} -- ce script refusera de se rejouer.")
    print(f"  Relancer maintenant LA MEME commande que le balayage, sans "
          f"plafond :")
    print(f"    python -u equivariance_f.py cicyquotients.m cicylist.txt "
          f"{args.output_dir} --replier-orbites -j 7\n")


if __name__ == '__main__':
    principal()
