#!/usr/bin/env python3
"""
compter_strates.py -- Combien de lignes le branchement du lieu de base
va-t-il reellement toucher, et dans quelles strates ?

POURQUOI CETTE MESURE EXISTE
----------------------------
Le branchement des §5.37 / §5.38 dans `equivariance_f.analyser` s'insere la
ou `f_sans_point_base` echoue. Son cout est donc borne par le nombre de
lignes qui atteignent ce point -- pas par les 505 601 lignes du fichier, et
pas non plus par les 978 verdicts deja obtenus hors balayage.

`lieu_de_base_rv3.analyser` balaie t dans [0, p) avec p ~ 30 011 pour trouver
les racines du cubique : ~30 000 evaluations de monomes par ligne. Si les
lignes concernees se comptent par milliers, c'est negligeable devant les
~30 h du balayage. Si elles se comptent par centaines de milliers, il faut
plafonner. La difference se mesure en une passe de lecture ; elle ne
s'estime pas -- c'est la lecon du §5.37, ou une demi-heure de mesure a
remplace trente heures extrapolees d'une moyenne prise sur les cas legers.

CE QUE CE SCRIPT DECLARE, ET DANS LES DEUX SENS
-----------------------------------------------
Il ne compte pas seulement la cible. Il ventile TOUTES les lignes `ok` selon
ce qui les arrete, de sorte que la cible se lise comme une part d'un tout et
non comme un nombre isole :

    survivantes / tuees par h0 / tuees par Hoppe / CIBLE / autres

et croise la CIBLE par (rank_C, rang_V), parce que seules deux strates ont
un module capable de trancher. Ce qui tombera en dehors est le silence a
surveiller : ces lignes-la recevront `forme inattendue`, et il faut savoir
combien AVANT de lire un compteur a zero comme une absence.

Lecture seule. N'ecrit rien, ne modifie rien.

Usage :
    python -u compter_strates.py scan_wilson5
    python -u compter_strates.py scan_wilson5\\results_equivariance_f.jsonl
"""
import os
import sys
import json
import argparse
from collections import Counter

# Les deux seules strates pour lesquelles un module existe.
#   (rank_C, rang_V) -> module
STRATES_TRAITEES = {
    (1, 3): 'lieu_de_base_rv3.analyser',
    (2, 3): 'lieu_de_base_rc2.analyser_rc2',
}


def compter(chemin):
    tot = Counter()
    cible_par_strate = Counter()
    cible_candidats = set()
    autres_par_strate = Counter()
    n_lignes = n_illisibles = 0

    with open(chemin, encoding='utf-8') as f:
        for ligne in f:
            if not ligne.strip():
                continue
            n_lignes += 1
            try:
                d = json.loads(ligne)
            except json.JSONDecodeError:
                n_illisibles += 1
                continue

            if d.get('etat') != 'ok':
                tot['hors etat ok'] += 1
                continue

            # La strate. rank_C est la longueur de c_charges ; rang_V est
            # ecrit par `analyser`, avec repli sur le champ du candidat.
            c = d.get('c_charges')
            if isinstance(c, str):
                try:
                    c = json.loads(c)
                except json.JSONDecodeError:
                    c = None
            rank_c = len(c) if isinstance(c, list) else None
            rang_v = d.get('rang_V', d.get('rank_V'))
            strate = (rank_c, rang_v)

            if d.get('survit'):
                tot['survivantes'] += 1
            elif d.get('h0_equivariant') not in (0, None) \
                    or d.get('h0_generique') not in (0, None):
                tot['tuees par h0'] += 1
            elif d.get('hoppe_complet') is False:
                tot['tuees par Hoppe'] += 1
            elif d.get('hoppe_complet') is True \
                    and d.get('surjectif_certifie') is not True:
                # LA CIBLE : Hoppe passe, la surjectivite n'est pas
                # certifiee. C'est exactement l'entree du branchement.
                tot['CIBLE (Hoppe ok, surjectivite non certifiee)'] += 1
                cible_par_strate[strate] += 1
                cible_candidats.add((d.get('cicy'),
                                     json.dumps(d.get('b_charges')),
                                     json.dumps(d.get('c_charges'))))
            else:
                tot['autres (Hoppe indetermine, etc.)'] += 1
                autres_par_strate[strate] += 1

    return (tot, cible_par_strate, cible_candidats, autres_par_strate,
            n_lignes, n_illisibles)


def principal():
    ap = argparse.ArgumentParser()
    ap.add_argument('cible', help='dossier de scan, ou chemin du JSONL')
    args = ap.parse_args()

    chemin = args.cible
    if os.path.isdir(chemin):
        chemin = os.path.join(chemin, 'results_equivariance_f.jsonl')
    if not os.path.exists(chemin):
        print(f"  introuvable : {chemin}")
        return 1

    taille = os.path.getsize(chemin) / (1 << 20)
    print(f"  {chemin}  ({taille:.0f} Mo)\n", flush=True)

    (tot, cible, cands, autres, n_lignes, n_illisibles) = compter(chemin)

    print(f"  {n_lignes} lignes lues"
          + (f", {n_illisibles} ILLISIBLES" if n_illisibles else ""))
    print(f"\n{'=' * 70}\n  VENTILATION DES LIGNES `ok`\n{'=' * 70}")
    for k, v in tot.most_common():
        print(f"    {v:>9}  {k}")

    n_cible = sum(cible.values())
    print(f"\n{'=' * 70}\n  LA CIBLE, PAR STRATE (rank_C, rang_V)\n{'=' * 70}")
    if not cible:
        print("    aucune -- le branchement ne serait JAMAIS appele.")
    for (rc, rv), v in sorted(cible.items(), key=lambda x: -x[1]):
        mod = STRATES_TRAITEES.get((rc, rv))
        etiq = mod if mod else 'AUCUN MODULE -> `forme inattendue`'
        print(f"    {v:>9}  rank_C={rc}, rang_V={rv}   {etiq}")

    traitees = sum(v for k, v in cible.items() if k in STRATES_TRAITEES)
    hors = n_cible - traitees
    print(f"\n    {traitees} lignes tombent dans une strate traitee, "
          f"{hors} hors des deux.")
    print(f"    {len(cands)} identites (cicy, b, c) distinctes dans la cible.")

    # Le cout : seules les lignes d'une strate traitee paient le balayage
    # des racines. Les autres sortent au test de forme, en microsecondes.
    print(f"\n{'=' * 70}\n  COUT ATTENDU\n{'=' * 70}")
    print(f"    ~30 000 evaluations de monomes par ligne traitee.")
    print(f"    {traitees} lignes traitees. Les {hors} autres sortent au")
    print(f"    test de forme (porteurs, amb[k], partition des colonnes),")
    print(f"    qui est purement combinatoire.")

    print(f"\n{'=' * 70}\n  LE RESTE, PAR STRATE (ni cible, ni tue, ni survivant)\n{'=' * 70}")
    for (rc, rv), v in sorted(autres.items(), key=lambda x: -x[1])[:10]:
        print(f"    {v:>9}  rank_C={rc}, rang_V={rv}")
    print()
    return 0


if __name__ == '__main__':
    sys.exit(principal())
