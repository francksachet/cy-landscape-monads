#!/usr/bin/env python3
"""
empreinte_code.py -- Quelle version du code a produit cette ligne ?

POURQUOI
--------
L'empreinte du checkpoint de `equivariance_f.py` couvre le fichier d'entree
et quatre options -- pas le CODE. Une reprise reconduit donc sans un mot des
lots calcules par une version anterieure du programme.

Ce n'est pas theorique. Dans `scan_wilson4`, la tranche 0 de plusieurs taches
portait « hors domaine (modele S/I non valide) » pendant que les tranches 1 a
3 de la MEME tache, memes charges, produisaient 192 verdicts « ok ». Ces
candidats ont deux c_charges : ils tombaient sous l'ancien `rank_c_max = 1`,
contrainte levee depuis. Les lignes de la tranche 0 etaient des reliques d'un
code d'avant, reconduites a chaque reprise parce que leur lot etait enregistre
et leur compte de lignes juste. Mesure : 41 identites (CICY, b, c) portant a
la fois « hors domaine » et « ok » -- une contradiction interne, dans un
fichier de resultats.

CE QUE CE MODULE FAIT, ET CE QU'IL NE FAIT PAS
----------------------------------------------
Il DECLARE. Chaque ligne ecrite porte `_code`, et la reprise affiche la
repartition des lots par version.

Il ne REFUSE pas, et l'empreinte du code n'entre surtout PAS dans celle du
checkpoint. Un checkpoint invalide fait EFFACER le JSONL (`equivariance_f`,
branche « checkpoint present mais inutilisable ») : y mettre le code
signifierait qu'une correction d'une ligne de commentaire detruit trente
heures de calcul. Le remede serait pire que le mal. On mesure, on affiche,
et on decide -- `retirer_lots.py` sait retirer selectivement.

NORMALISATION DES FINS DE LIGNE
-------------------------------
Les octets sont normalises CRLF -> LF avant hachage. Sans cela, ce depot --
ou cinq fichiers suivis different de HEAD par leurs seules fins de ligne --
verrait son empreinte changer sans qu'une instruction ait bouge.
"""
import os
import hashlib

# Tout ce dont une modification peut changer un verdict.
DOSSIERS = ('cy_landscape/core',)
FICHIERS = ('equivariance_f.py', 'wilson_match.py',
            'cy_landscape/data/parse_oxford.py')


def fichiers_surveilles(racine=None):
    """Chemins relatifs, tries. Les absents sont ignores en silence ici --
    `empreinte_code` les compte et les rend, c'est la que ca se declare."""
    racine = racine or os.path.dirname(os.path.abspath(__file__))
    out = []
    for d in DOSSIERS:
        plein = os.path.join(racine, d)
        if not os.path.isdir(plein):
            continue
        for nom in os.listdir(plein):
            if nom.endswith('.py'):
                out.append(os.path.join(d, nom).replace('\\', '/'))
    for f in FICHIERS:
        if os.path.exists(os.path.join(racine, f)):
            out.append(f.replace('\\', '/'))
    return sorted(out)


def empreinte_code(racine=None, longueur=12):
    """
    Rend (empreinte, n_fichiers). L'empreinte hache le CHEMIN et le CONTENU
    de chaque fichier surveille : renommer un module la change aussi.
    """
    racine = racine or os.path.dirname(os.path.abspath(__file__))
    h = hashlib.sha256()
    fichiers = fichiers_surveilles(racine)
    for rel in fichiers:
        h.update(rel.encode('utf-8'))
        with open(os.path.join(racine, rel), 'rb') as f:
            h.update(f.read().replace(b'\r\n', b'\n'))
    return h.hexdigest()[:longueur], len(fichiers)


def repartition(chemin_jsonl):
    """
    {empreinte de code: nombre de lignes} sur un JSONL de resultats.
    `None` = ligne ecrite avant l'introduction du marquage.

    C'est la mesure qui manquait : sans elle, un fichier melangeant trois
    versions du code se lit comme un fichier homogene.
    """
    import json
    from collections import Counter
    c = Counter()
    with open(chemin_jsonl, encoding='utf-8') as f:
        for ligne in f:
            if not ligne.strip():
                continue
            try:
                c[json.loads(ligne).get('_code')] += 1
            except json.JSONDecodeError:
                c['__illisible__'] += 1
    return c


if __name__ == '__main__':
    import sys
    emp, n = empreinte_code()
    print(f"  empreinte du code : {emp}  ({n} fichiers surveilles)")
    for chemin in sys.argv[1:]:
        print(f"\n  {chemin}")
        for k, v in sorted(repartition(chemin).items(),
                           key=lambda x: -x[1]):
            etiquette = k if k else '(avant le marquage)'
            print(f"    {v:>8} lignes   {etiquette}")
