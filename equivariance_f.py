#!/usr/bin/env python3
"""
equivariance_f.py -- Test d'equivariance de f : B -> C, et stabilite restreinte.

--------------------------------------------------------------------------
Ce que ce script ajoute a `equivariance.py`
--------------------------------------------------------------------------
`equivariance.py` teste une condition NECESSAIRE portant sur les CHARGES :
sigma doit permuter l'ensemble des b_i et des c_j. Elle ne dit rien des
polynomes, et elle est vide sur les CICYs ou Gamma n'agit que par des phases.

Ce script descend au niveau des polynomes, en trois etapes :

  1. POLYNOMES COVARIANTS. Les coefficients des K polynomes definissants
     sont resolus pour que Gamma preserve l'ideal (`covariant_ring`), au lieu
     d'etre tires au hasard. Sans cela l'action ne descend pas au quotient et
     tout ce qui suit serait un calcul sur une autre variete.

  2. EXISTENCE D'UN f EQUIVARIANT. On resout S_g(f_{j,i}) = lambda_g
     f_{rho(j),pi(i)} pour tous les generateurs (`equivariant_monad`), en
     ENUMERANT les lambda_g admissibles a partir de l'ordre projectif de
     l'operateur -- ce qui traite le relevement projectif au lieu de le
     supposer trivial.

  3. STABILITE RESTREINTE. C'est l'etape qui mord. Sur un Gamma qui agit par
     phases, le sous-espace equivariant represente une fraction 1/|Gamma| de
     l'espace des f et n'est donc JAMAIS vide : conclure de l'etape 2 que le
     fibre descend serait une erreur. On recalcule donc h^0(V) avec un f tire
     dans le sous-espace equivariant. Si h^0(V) devient non nul, V n'est plus
     stable et le candidat tombe.

--------------------------------------------------------------------------
Domaine
--------------------------------------------------------------------------
Les etapes 2 et 3 utilisent le modele R_a = S_a / I_a, qui ne represente
H^0(Y, O(a)) que sur le domaine verifie par `sections.domaine_valide`.
Hors de ce domaine le script n'affiche pas de verdict : il marque
`hors domaine`. Un candidat hors domaine n'est ni retenu ni elimine.

Usage:
    python equivariance_f.py cicyquotients.m cicylist.txt scan_wilson2
    python equivariance_f.py cicyquotients.m cicylist.txt scan_wilson2 --cicy 6947
"""
import os
import sys
import json
import argparse
import hashlib

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cy_landscape.core.braun_symmetry import (parse_symmetries, ordres_rt,
                                              matrice_mod_p)
from cy_landscape.core.gamma_action import choisir_premier, racine_primitive
from cy_landscape.core.covariant_ring import (resoudre_covariants,
                                              tirer_covariants,
                                              verifier_covariance,
                                              CovariantRing)
from cy_landscape.core.equivariant_monad import (espace_f_equivariant,
                                                 h0_V_sur_espace,
                                                 h0_wedge2_V_sur_espace,
                                                 hoppe_sur_espace,
                                                 f_sans_point_base)
from cy_landscape.core.sections import (basis_multi, rref_mod, _mult_matrix,
                                        domaine_valide)
from cy_landscape.data.parse_oxford import load_oxford_file
from wilson_match import parse_braun, parse_cicylist, apparier


def h0_V_generique(anneau, b_charges, c_charges, p, rng, n_essais=5):
    """h^0(V) pour un f generique, dans le MEME anneau covariant.

    C'est la valeur de reference : la comparer a `h0_V_sur_espace` isole
    l'effet de la contrainte d'equivariance, sans melanger l'effet du
    passage aux polynomes covariants.
    """
    m = len(b_charges[0])
    dsrc = sum(anneau.dimY(list(b)) for b in b_charges)
    meilleur = None
    for _ in range(n_essais):
        lignes = []
        for cj in c_charges:
            ddst = anneau.dimY(list(cj))
            ligne = []
            for b in b_charges:
                deg = [cj[k] - b[k] for k in range(m)]
                if any(x < 0 for x in deg):
                    ligne.append(np.zeros((ddst, anneau.dimY(list(b))),
                                          dtype=np.int64))
                    continue
                fb = basis_multi(anneau.amb, deg)
                ligne.append(_mult_matrix(anneau, list(b), deg,
                                          (fb, rng.randint(1, p, size=len(fb))),
                                          list(cj)))
            lignes.append(np.hstack(ligne))
        rang, _ = rref_mod(np.vstack(lignes).T.copy(), p)
        h0 = dsrc - rang
        meilleur = h0 if meilleur is None else min(meilleur, h0)
    return meilleur


def ordre_nom(nom):
    """Ordre du groupe lu sur son nom de Braun ("Z2 x Z2" -> 4)."""
    import re as _re
    n = 1
    for x in _re.findall(r'Z(\d+)', str(nom)):
        n *= int(x)
    return n


def n_gen_quotient(cohomology, nom_groupe):
    """
    Nombre de generations sur X/Gamma : |chi(V)| / |Gamma|.

    Gamma agissant librement, n_gen(X/Gamma) = n_gen(X)/|Gamma|. Sans ce
    nombre, un verdict `SURVIT` ne dit RIEN du contenu physique : un fibre
    peut etre stable, equivariant et surjectif tout en donnant 12
    generations. Cas reel : #21, SU(5) de rang 5, |chi| = 24 avec Gamma =
    Z2 -- douze generations, et pourtant `SURVIT`.
    """
    if not cohomology or len(cohomology) < 4:
        return None
    chi = abs(cohomology[1] - cohomology[2])
    o = ordre_nom(nom_groupe)
    if o <= 0 or chi == 0 or chi % o:
        return None
    return chi // o


def analyser(cicy_num, amb, cfg, b, c, symetries, groupes=None, graine=0):
    """Renvoie une liste de lignes de resultat, une par (symetrie, lambda)."""
    lignes = []
    # rank_c_max=None : on autorise rank_C >= 2. Les fonctions generalisees
    # (h^0(V), decomposition de H^1(V)) s'y appliquent ; celles qui supposent
    # rank_C = 1 -- wedge^p V, surjectivite -- se declarent elles-memes non
    # calculables et le verdict tombe en `indetermine`, jamais en succes.
    if not domaine_valide(amb, cfg, b, c, rank_c_max=None):
        return [{'groupe': '-', 'etat': 'hors domaine (modele S/I non valide)'}]

    for sym in symetries:
        if groupes and sym['nom'] not in groupes:
            continue
        ordres = ordres_rt(sym['coord']) | ordres_rt(sym['poly']) | {2}
        try:
            p, _ = choisir_premier(sorted(ordres), minimum=30011)
            rac = {n: racine_primitive(p, n) for n in ordres}
            Mc = [matrice_mod_p(x, p, rac) for x in sym['coord']]
            Np = [matrice_mod_p(x, p, rac) for x in sym['poly']]
        except Exception as exc:
            lignes.append({'groupe': sym['nom'],
                           'etat': f'generateurs illisibles ({exc})'})
            continue

        res = resoudre_covariants(amb, cfg, Mc, Np, p)
        if res is None:
            lignes.append({'groupe': sym['nom'], 'etat': 'sigma non extractible'})
            continue
        bons = [k for k, v in res['par_convention'].items() if v['non_degenere']]
        if 'N' not in bons:
            lignes.append({'groupe': sym['nom'],
                           'etat': 'aucun ideal covariant non degenere'})
            continue
        v = res['par_convention']['N']
        co = tirer_covariants(v['base'], res['offsets'], res['dims'], p,
                              np.random.RandomState(graine))
        ok, ec = verifier_covariance(amb, cfg, Mc, Np, 'N', co, p)
        if not ok:
            lignes.append({'groupe': sym['nom'],
                           'etat': f'covariance non revérifiée (ecart {ec})'})
            continue

        anneau = CovariantRing(amb, cfg, co, p)
        out = espace_f_equivariant(anneau, amb, b, c, Mc, p)
        if out['etat'] != 'ok':
            lignes.append({'groupe': sym['nom'], 'etat': out['etat']})
            continue
        if not out['solutions']:
            lignes.append({'groupe': sym['nom'], 'etat': 'aucun f equivariant',
                           'dim_totale': out['dim_totale'], 'elimine': True})
            continue

        rang_V = len(b) - len(c)
        # Hoppe, c1(V) = 0 : rk 3 -> h0(V), h3(V) ; rk >= 4 -> + h0(w2V).
        # w2V n'est calculable par ce chemin que si rank_C = 1, comme dans
        # hoppe_fast. Au-dela on ne conclut pas.
        besoin_w2 = rang_V >= 4
        w2_possible = besoin_w2 and len(c) == 1

        h_gen = h0_V_generique(anneau, b, c, p, np.random.RandomState(graine + 5))
        base_tot = np.eye(out['dim_totale'], dtype=np.int64)
        w2_gen = None
        if w2_possible:
            w2_gen, _ = h0_wedge2_V_sur_espace(
                anneau, b, c, base_tot, out['offsets'], out['dims'],
                out['degres'], p, np.random.RandomState(graine + 5))

        for s in out['solutions']:
            h_eq, _ = h0_V_sur_espace(anneau, amb, b, c, s['base'],
                                      out['cases'], out['offsets'], out['dims'],
                                      out['degres'], p,
                                      np.random.RandomState(graine + 5))
            w2_eq = None
            if w2_possible:
                w2_eq, _ = h0_wedge2_V_sur_espace(
                    anneau, b, c, s['base'], out['offsets'], out['dims'],
                    out['degres'], p, np.random.RandomState(graine + 5))
            lam = tuple(int(x) if x < p // 2 else int(x) - p for x in s['lambda'])

            # « survit » exige TOUS les tests disponibles a ce rang. Un test
            # non calculable (w2 hors taille, ou rank_C >= 2) rend le verdict
            # indetermine -- jamais favorable.
            if h_gen != 0 or h_eq != 0:
                survit, indetermine = False, False
            elif not besoin_w2:
                survit, indetermine = True, False
            elif not w2_possible or w2_eq is None or w2_gen is None:
                survit, indetermine = False, True
            else:
                survit = (w2_gen == 0 and w2_eq == 0)
                indetermine = False

            # Surjectivite : V = ker(f) n'est un FIBRE que si f est surjective
            # en tout point. Rien ne le garantit sur le sous-espace
            # equivariant. Le critere est SUFFISANT (voir f_sans_point_base) :
            # un echec ne prouve pas la non-surjectivite, il empeche seulement
            # de conclure -- d'ou `indetermine` et non une elimination.
            # Teste uniquement quand tout le reste passe, pour le cout.
            # Critere de Hoppe COMPLET sur le sous-espace : h0(wedge^p V) = 0
            # pour p = 1..rk-1. Inclut h^3(V) = h0(wedge^{rk-1} V), qui n'etait
            # jusqu'ici traite nulle part sous contrainte. Un p non calculable
            # rend le verdict indetermine, jamais favorable.
            hoppe = None
            if survit and len(c) == 1:
                hoppe = hoppe_sur_espace(
                    anneau, b, c, s['base'], out['offsets'], out['dims'],
                    out['degres'], p, np.random.RandomState(graine + 5))
                if hoppe['stable'] is not True:
                    survit = False
                    indetermine = (hoppe['stable'] is None)

            surj = None
            if survit:
                surj = f_sans_point_base(
                    anneau, b, c, s['base'], out['offsets'], out['dims'],
                    out['degres'], p, np.random.RandomState(graine + 7),
                    n_essais=2, n_degres=8)
                if not surj['certifie']:
                    survit, indetermine = False, True

            lignes.append({
                'groupe': sym['nom'], 'etat': 'ok', 'lambda': lam,
                'dim_equivariant': s['dim'], 'dim_totale': out['dim_totale'],
                'rang_V': rang_V,
                'h0_generique': h_gen, 'h0_equivariant': h_eq,
                'h0w2_generique': w2_gen, 'h0w2_equivariant': w2_eq,
                'hoppe_complet': None if hoppe is None else hoppe['stable'],
                'hoppe_valeurs': None if hoppe is None else
                                 {str(k): v for k, v in hoppe['valeurs'].items()},
                'surjectif_certifie': None if surj is None else bool(surj['certifie']),
                'surjectif_degre': None if surj is None else surj['degre'],
                'survit': bool(survit), 'indetermine': bool(indetermine),
            })
    return lignes


# ======================================================================
# Reprise sur checkpoint
# ======================================================================
#
# POURQUOI
# --------
# La version precedente accumulait toutes les lignes dans une liste et
# n'ecrivait le JSONL qu'apres le DERNIER candidat. Tant que le lot faisait
# 108 candidats (une heure), c'etait sans consequence. Le generateur enumere
# du §5.23 en produit 14 945 : le meme calcul demande une cinquantaine
# d'heures, et une interruption -- Ctrl+C, redemarrage, coupure -- perdait
# la totalite du travail, sans qu'aucun fichier n'en garde trace.
#
# Mesure reelle : deux heures vingt de calcul, zero octet recuperable.
#
# COMMENT
# -------
# Meme principe que `main_optimized` : un JSONL en ecriture « append-only »,
# et un fichier de progression leger ecrit APRES chaque candidat.
#
# Le compteur ne suffit pas. Si le processus meurt AU MILIEU d'un candidat,
# quelques lignes de ce candidat sont deja dans le JSONL : reprendre au
# candidat suivant les laisserait en double, reprendre au meme les
# dupliquerait. On enregistre donc aussi l'OFFSET du fichier apres le
# dernier candidat complet, et la reprise tronque a cet offset. Ce qui est
# relu est alors exactement ce qui a ete valide.
#
# L'empreinte du fichier d'entree (sha256) interdit de reprendre un
# checkpoint sur un autre lot : un decalage d'indices produirait des
# resultats attribues aux mauvais candidats, ce qui est pire que pas de
# reprise du tout.

def _empreinte(chemin, filtre_cicy):
    h = hashlib.sha256()
    with open(chemin, 'rb') as f:
        for bloc in iter(lambda: f.read(1 << 20), b''):
            h.update(bloc)
    h.update(str(filtre_cicy).encode())
    return h.hexdigest()


class Progression:
    """Progression legere : indice du prochain candidat + offset du JSONL."""

    def __init__(self, dossier, empreinte):
        self.chemin = os.path.join(dossier, 'progress_equivariance_f.json')
        self.empreinte = empreinte
        self.fait = 0
        self.offset = 0
        self.compteurs = {'survivants': 0, 'indetermines': 0, 'ecartes': 0}

    def charger(self, taille_jsonl):
        if not os.path.exists(self.chemin):
            return False, "aucun checkpoint"
        try:
            with open(self.chemin, encoding='utf-8') as f:
                d = json.load(f)
        except (json.JSONDecodeError, OSError) as e:
            return False, f"checkpoint illisible ({type(e).__name__})"
        if d.get('empreinte') != self.empreinte:
            return False, ("le fichier d'entree a change depuis le checkpoint "
                           "(ou --cicy differe) -- reprise refusee")
        if d.get('offset', 0) > taille_jsonl:
            return False, ("le JSONL est plus court que le checkpoint "
                           "-- reprise refusee")
        self.fait = int(d.get('fait', 0))
        self.offset = int(d.get('offset', 0))
        self.compteurs.update(d.get('compteurs') or {})
        return True, f"reprise au candidat {self.fait}"

    def sauver(self, fait, offset, compteurs):
        self.fait, self.offset = fait, offset
        self.compteurs = dict(compteurs)
        tmp = self.chemin + '.tmp'
        with open(tmp, 'w', encoding='utf-8') as f:
            json.dump({'empreinte': self.empreinte, 'fait': fait,
                       'offset': offset, 'compteurs': self.compteurs}, f)
        os.replace(tmp, self.chemin)   # remplacement atomique


# ======================================================================
# Repli par orbite sous Aut(matrice de configuration)
# ======================================================================
#
# Le generateur enumere (§5.23) produit la famille COMPLETE des sommes de
# vecteurs unite. Sur une CICY dont la matrice de configuration est
# symetrique, cette famille contient les images les unes des autres : les
# 12 monades survivantes de #6947 sont les 12 arrangements d'un meme motif
# (3, 1, 1, 0) sur ses quatre facteurs P^1, et |Aut| y vaut 24.
#
# Mesure sur les 14 945 candidats de scan_wilson4 : 3 688 orbites, soit un
# facteur 4,05 -- environ 14 h au lieu de 55.
#
# CE QUI EST DEMONTRE, ET CE QUI NE L'EST PAS
# -------------------------------------------
# Qu'une permutation preserve la matrice de configuration n'implique PAS
# qu'elle commute avec l'action de Gamma lue chez Braun, laquelle est
# attachee a des coordonnees precises. L'egalite des verdicts sur une
# orbite est donc une hypothese, verifiee empiriquement (0 discordance sur
# les 919 lignes de #6890, #6947 et #6715) mais non prouvee.
#
# D'ou deux precautions, sans lesquelles ce repli serait une redite exacte
# du §5.23 -- un filtre qui fait disparaitre des candidats en silence :
#
#   1. AUCUNE LIGNE NE DISPARAIT. Le verdict du representant est recopie
#      sur chaque membre, avec `verdict_replique: True`. Le JSONL de sortie
#      compte autant de lignes qu'un balayage complet, et l'aval ne voit
#      aucune difference -- sinon un champ qui dit d'ou vient le verdict.
#   2. LE REPLI SE CONTROLE LUI-MEME. `--controle-orbites N` evalue POUR DE
#      VRAI N membres non representants et compare. Une discordance est
#      affichee, comptee, et inscrite dans le JSONL.
#
# Le repli reste optionnel : par defaut, chaque candidat est evalue.

def _construire_taches(rs, entries, args):
    """
    Rend (taches, rapport). Une tache = (i_representant, [i_membres]).

    Sans `--replier-orbites`, une tache par candidat : la boucle
    principale ne connait que les taches et ne change pas de forme.
    """
    if not getattr(args, 'replier_orbites', False):
        return [(i, [i]) for i in range(len(rs))], None

    from cy_landscape.core.symetrie_config import automorphismes, canonique

    par_cicy = {}
    for i, r in enumerate(rs):
        par_cicy.setdefault(r['cicy'], []).append(i)

    taches = []
    n_aut_trivial = 0
    cache = {}
    for cicy, idx in par_cicy.items():
        e = entries.get(cicy)
        if e is None:
            taches.extend((i, [i]) for i in idx)
            continue
        if cicy not in cache:
            cache[cicy] = automorphismes(e['ambient'], e['config'])
        autos, complet = cache[cicy]
        if len(autos) <= 1:
            n_aut_trivial += 1
            taches.extend((i, [i]) for i in idx)
            continue
        classes = {}
        for i in idx:
            r = rs[i]
            if not r.get('b_charges'):
                classes[('_sans_charges_', i)] = [i]
                continue
            classes.setdefault(
                canonique(r['b_charges'], r['c_charges'], autos), []).append(i)
        for membres in classes.values():
            taches.append((membres[0], membres))

    rapport = {
        'candidats': len(rs),
        'taches': len(taches),
        'cicys_sans_symetrie': n_aut_trivial,
        'facteur': len(rs) / max(1, len(taches)),
    }
    return taches, rapport


def _echantillon_controle(taches, n, graine=0):
    """
    Rend {indice de tache: [membres a reevaluer pour de vrai]}.

    On tire N COUPLES (orbite, membre), pas N orbites. La difference n'est
    pas cosmetique : avec un tirage par orbite, un repli abusif qui range
    TOUS les candidats d'une CICY dans une seule orbite ne recevait qu'UN
    controle -- et passait. Constate en sabotant `canonique` pour qu'elle
    renvoie une constante : 176 candidats replies en 1 tache, 1 controle,
    0 discordance. Le controle validait un repli entierement faux.

    En tirant des couples, une orbite geante recoit une part des controles
    proportionnelle a sa taille, et le sabotage tombe.

    Les orbites a un seul membre sont exclues : les controler ne
    comparerait rien -- le genre de controle qui passe toujours (§8).
    """
    if n <= 0:
        return {}
    couples = [(k, j) for k, (i_rep, membres) in enumerate(taches)
               for j in membres if j != i_rep]
    if not couples:
        return {}
    rng = np.random.RandomState(graine)
    if len(couples) > n:
        pris = [couples[t] for t in
                rng.choice(len(couples), size=n, replace=False)]
    else:
        pris = couples
    sortie = {}
    for k, j in pris:
        sortie.setdefault(k, []).append(j)
    return sortie


def _sortie_tolerante():
    """
    Empeche un plantage d'encodage sur une console Windows.

    Les etiquettes de groupe de jauge contiennent des indices Unicode
    ("E₆"), que la console cp1252 ne sait pas encoder : `print` levait
    UnicodeEncodeError des la premiere ligne de resultat, apres plusieurs
    minutes de calcul et sans rien ecrire. On passe donc stdout en
    errors='replace' -- le caractere devient '?', le calcul continue. Le
    fichier JSONL de sortie, lui, est ecrit en UTF-8 explicite et garde
    l'etiquette exacte.
    """
    for flux in (sys.stdout, sys.stderr):
        try:
            flux.reconfigure(errors='replace')
        except (AttributeError, ValueError):
            pass


def main():
    _sortie_tolerante()
    ap = argparse.ArgumentParser()
    ap.add_argument('braun_m')
    ap.add_argument('cicylist')
    ap.add_argument('output_dir')
    ap.add_argument('--input', default='results_equivariant.jsonl')
    ap.add_argument('--cicy', type=int, default=None)
    ap.add_argument('--tous-groupes', action='store_true',
                    help="Ne pas se limiter aux groupes d'ordre compatible "
                         "avec l'indice (champ groupes_utiles).")
    ap.add_argument('--reset', action='store_true',
                    help="Ignorer le checkpoint et repartir de zero.")
    ap.add_argument('--arret-apres', type=int, default=0,
                    help="S'arreter proprement apres N taches, comme le "
                         "ferait une interruption. Sert aux tests : une "
                         "coupure provoquee par un chronometre tombe presque "
                         "toujours ENTRE deux candidats, et la reprise n'est "
                         "alors pas reellement mise a l'epreuve. 0 = inactif.")
    ap.add_argument('--replier-orbites', action='store_true',
                    help="N'evaluer qu'un representant par orbite sous le "
                         "groupe d'automorphismes de la matrice de "
                         "configuration, et recopier son verdict sur les "
                         "autres membres (champ verdict_replique). Mesure "
                         "sur les 14 945 candidats de scan_wilson4 : 3 688 "
                         "orbites, soit un facteur 4,05. Les lignes de "
                         "sortie restent AU COMPLET.")
    ap.add_argument('--controle-orbites', type=int, default=20,
                    help="Avec --replier-orbites : evaluer pour de vrai N "
                         "membres non representants, et comparer au verdict "
                         "recopie. Sans ce controle le repli serait une "
                         "hypothese invisible -- exactement le mecanisme du "
                         "§5.23. 0 = desactiver.")
    args = ap.parse_args()

    entries = {e['num']: e for e in load_oxford_file(args.cicylist)}
    braun = parse_braun(args.braun_m)
    cl = parse_cicylist(args.cicylist)
    corr, _, _ = apparier(braun, cl)
    inv = {v: k for k, v in corr.items()}
    SYM = parse_symmetries(args.braun_m)

    src = os.path.join(args.output_dir, args.input)
    rs = [json.loads(l) for l in open(src, encoding='utf-8') if l.strip()]
    if args.cicy:
        rs = [r for r in rs if r['cicy'] == args.cicy]

    # ---------------- repli par orbite (optionnel) -------------------
    # Une TACHE = (indice du representant, indices de tous les membres).
    # Sans repli, une tache par candidat : la boucle et le checkpoint ne
    # connaissent que les taches, et rien d'autre ne change.
    taches, rapport_orbites = _construire_taches(rs, entries, args)

    # ---------------- checkpoint ------------------------------------
    dst = os.path.join(args.output_dir, 'results_equivariance_f.jsonl')
    prog = Progression(args.output_dir,
                       _empreinte(src, (args.cicy, bool(args.replier_orbites),
                                        int(args.controle_orbites))))
    depart = 0
    survivants = indetermines = ecartes = 0

    if args.reset:
        for p in (dst, prog.chemin):
            if os.path.exists(p):
                os.remove(p)
        print(f"\n  Checkpoint reinitialise (--reset)")
    else:
        taille = os.path.getsize(dst) if os.path.exists(dst) else 0
        ok, motif = prog.charger(taille)
        if ok:
            depart = prog.fait
            survivants = prog.compteurs['survivants']
            indetermines = prog.compteurs['indetermines']
            ecartes = prog.compteurs['ecartes']
            # Tronque ce qu'un candidat interrompu a pu laisser derriere lui.
            with open(dst, 'r+b') as fh:
                fh.truncate(prog.offset)
            print(f"\n  Reprise : {depart} / {len(taches)} taches deja "
                  f"traitees ({100.0 * depart / max(1, len(taches)):.1f} %)")
        elif os.path.exists(prog.chemin) or taille:
            # Un checkpoint existait mais n'est pas utilisable : le DIRE.
            # Repartir de zero en silence ferait passer un recommencement
            # complet pour une reprise.
            print(f"\n  Checkpoint present mais inutilisable : {motif}")
            print(f"  Reprise depuis le debut ; l'ancien JSONL est ecrase.")
            if os.path.exists(dst):
                os.remove(dst)

    print(f"\n{'=' * 96}")
    print("  EQUIVARIANCE DE f  --  polynomes covariants, puis stabilite restreinte")
    print(f"{'=' * 96}")
    print(f"  {'CICY':>5} {'jauge':>7} {'rk':>2} {'groupe':<11} {'lambda':>14} "
          f"{'dim eq':>6} {'N':>5} {'h0 gen':>6} {'h0 eq':>5} "
          f"{'w2 gen':>6} {'w2 eq':>5}  verdict")

    # Un candidat -> ses lignes de sortie. Extrait de la boucle pour que le
    # point de sauvegarde soit UNIQUE et atteint par tous les chemins, y
    # compris ceux qui n'evaluent rien. Un `continue` qui saute le
    # checkpoint reperdrait le candidat a la reprise suivante.
    def traiter(r):
        lignes_out, d_surv, d_ind, d_ec = [], 0, 0, 0
        e = entries.get(r['cicy'])
        num_b = inv.get(r['cicy'])
        if e is None or num_b not in SYM:
            return lignes_out, d_surv, d_ind, d_ec
        amb, cfg = e['ambient'], np.asarray(e['config'])
        b = [list(x) for x in r['b_charges']]
        c = [list(x) for x in r['c_charges']]
        # PAS DE REPLI SILENCIEUX. La version precedente faisait :
        #     if groupes is not None and not groupes:
        #         groupes = set(r.get('equivariant_possible') or [])
        # c'est-a-dire : faute de groupe d'ordre compatible avec l'indice,
        # on les essayait TOUS, sans le dire. Mesure sur le balayage
        # precedent : 3 892 couples sur 4 076 (95,5 %) avaient un indice
        # incompatible -- ils ne peuvent donner trois generations avec ce
        # Gamma, quel que soit le verdict -- et certains ressortaient
        # etiquetes SURVIT. Un filtre qui devient vide sans le dire est le
        # meme defaut que le « zero exotique » du §4.8.
        #
        # Le candidat est desormais ECARTE, et la raison PERSISTEE : un
        # fichier de resultats doit dire pourquoi un cas n'a pas ete traite.
        groupes = None if args.tous_groupes else set(r.get('groupes_utiles') or [])
        if groupes is not None and not groupes:
            L = {'groupe': '-', 'etat': 'aucun groupe d ordre compatible avec '
                                        'l indice (|chi| != 3.|Gamma|)'}
            print(f"  {r['cicy']:>5} {r.get('gauge', ''):>7} "
                  f"{L['groupe']:<11} {L['etat']}")
            d_ec += 1
            lignes_out.append({**{k: r.get(k) for k in
                                  ('cicy', 'gauge', 'rank_V', 'cohomology',
                                   'b_charges', 'c_charges', 'groupes_utiles',
                                   'equivariant_possible', 'ordres_gamma')},
                               **L, 'survit': False, 'indetermine': True})
            return lignes_out, d_surv, d_ind, d_ec
        lignes = analyser(r['cicy'], amb, cfg, b, c, SYM[num_b]['symetries'],
                          groupes=groupes)
        # Identite du candidat, recopiee sur CHAQUE ligne de sortie.
        ident = {k: r.get(k) for k in
                 ('cicy', 'gauge', 'rank_V', 'cohomology',
                  'b_charges', 'c_charges', 'groupes_utiles',
                  'equivariant_possible', 'ordres_gamma')}

        for L in lignes:
            if L['etat'] != 'ok':
                # ECRITES ELLES AUSSI. La version precedente ne persistait que
                # les lignes 'ok' : les « hors domaine », « charges non
                # permutees » et « espace trop grand » n'existaient que dans la
                # sortie console. Sur le balayage `scan_gros_gamma`, cela a
                # rendu le JSONL trompeur -- il montrait 0 couple sur un groupe
                # d'ordre compatible, alors que 26 candidats en portaient un et
                # avaient simplement ete ecartes en amont, pour une raison que
                # le fichier ne contenait pas. Un fichier de resultats doit
                # dire pourquoi un cas n'a pas ete traite, sinon son silence se
                # lit a tort comme une absence de candidats.
                print(f"  {r['cicy']:>5} {r.get('gauge', ''):>7} "
                      f"{L['groupe']:<11} {L['etat']}")
                d_ec += 1
                lignes_out.append({**ident, **L, 'survit': False,
                                   'indetermine': True})
                continue
            L['n_gen_quotient'] = n_gen_quotient(r.get('cohomology'),
                                                 L['groupe'])
            if L['survit']:
                ng = L['n_gen_quotient']
                verdict = (f"SURVIT -- {ng if ng is not None else '?'} gen "
                           f"sur X/Gamma (Hoppe complet + surjectif en "
                           f"{L['surjectif_degre']})")
            elif L.get('hoppe_complet') is False:
                verdict = f"tue par Hoppe complet : {L.get('hoppe_valeurs')}"
            elif L.get('hoppe_complet') is None and L.get('indetermine') \
                    and L.get('surjectif_certifie') is None:
                verdict = "indetermine : Hoppe complet non calculable"
            elif L.get('surjectif_certifie') is False:
                verdict = "indetermine : surjectivite de f non certifiee"
            elif L.get('indetermine'):
                verdict = "indetermine (w2V non calculable)"
            elif L['h0_generique'] != 0:
                verdict = "deja non stable avec f generique"
            elif L['h0_equivariant'] != 0:
                verdict = "tue par h0(V) equivariant"
            elif L.get('h0w2_generique'):
                verdict = "deja non stable : h0(w2V) generique != 0"
            else:
                verdict = "tue par h0(w2V) equivariant"
            d_surv += bool(L['survit'])
            d_ind += bool(L.get('indetermine'))
            fmt = lambda x: '-' if x is None else str(x)
            print(f"  {r['cicy']:>5} {r.get('gauge', ''):>7} {L['rang_V']:>2} "
                  f"{L['groupe']:<11} {str(L['lambda']):>14} "
                  f"{L['dim_equivariant']:>6} {L['dim_totale']:>5} "
                  f"{L['h0_generique']:>6} {L['h0_equivariant']:>5} "
                  f"{fmt(L['h0w2_generique']):>6} "
                  f"{fmt(L['h0w2_equivariant']):>5}  {verdict}")
            lignes_out.append({**ident, **L})
        return lignes_out, d_surv, d_ind, d_ec

    # ---------------- boucle avec checkpoint -------------------------
    # 'a' et non 'w' : la reprise a deja tronque le fichier a l'offset du
    # dernier candidat complet, donc ouvrir en ecriture l'effacerait.
    interrompu = False
    controles = _echantillon_controle(taches, args.controle_orbites)
    n_controles = n_discordances = 0

    def _cle_verdict(x):
        # Ce qui doit coincider sur une orbite : le verdict, pas les
        # charges ni le degre temoin -- ceux-la sont PERMUTES par
        # l'automorphisme et differeraient legitimement.
        return (str(x.get('groupe')), str(x.get('lambda')),
                bool(x.get('survit')), str(x.get('etat')),
                str(x.get('n_gen_quotient')))

    with open(dst, 'a', encoding='utf-8') as fh:
        for k, (i_rep, membres) in enumerate(taches):
            if k < depart:
                continue
            if args.arret_apres and (k - depart) >= args.arret_apres:
                interrompu = True
                break
            try:
                bloc, ds, di, de = traiter(rs[i_rep])
                # --- controle du repli ---------------------------------
                for j in controles.get(k, ()):
                    bloc_j, _, _, _ = traiter(rs[j])
                    n_controles += 1
                    a = sorted(map(_cle_verdict, bloc))
                    b_ = sorted(map(_cle_verdict, bloc_j))
                    if a != b_:
                        n_discordances += 1
                        print(f"  !! DISCORDANCE D'ORBITE sur "
                              f"#{rs[i_rep]['cicy']} : le membre de controle "
                              f"ne donne pas le meme verdict que le "
                              f"representant.")
                        print(f"     representant : {a[:2]}")
                        print(f"     membre       : {b_[:2]}")
            except KeyboardInterrupt:
                interrompu = True
                break

            # --- ecriture : une ligne par MEMBRE, jamais une de moins ---
            for j in membres:
                replique = (j != i_rep)
                ident_j = {c: rs[j].get(c) for c in
                           ('cicy', 'gauge', 'rank_V', 'cohomology',
                            'b_charges', 'c_charges', 'groupes_utiles',
                            'equivariant_possible', 'ordres_gamma')}
                for x in bloc:
                    y = dict(x)
                    if replique:
                        # Charges du membre, verdict du representant. Le
                        # degre temoin, lui, reste celui du representant :
                        # il est permute par l'automorphisme, et le
                        # recopier tel quel serait faux -- on le marque.
                        y.update(ident_j)
                        y['verdict_replique'] = True
                        y['representant'] = rs[i_rep].get('b_charges')
                    else:
                        y['verdict_replique'] = False
                    fh.write(json.dumps(y, default=int) + '\n')
                if replique:
                    survivants += ds
                    indetermines += di
                    ecartes += de
            survivants += ds
            indetermines += di
            ecartes += de
            # L'ordre compte : vider le tampon, forcer l'ecriture disque,
            # PUIS enregistrer l'offset. L'inverse laisserait un checkpoint
            # pointant au-dela de ce qui est reellement sur le disque.
            fh.flush()
            os.fsync(fh.fileno())
            prog.sauver(k + 1, fh.tell(),
                        {'survivants': survivants,
                         'indetermines': indetermines, 'ecartes': ecartes})

    if interrompu:
        print(f"\n  INTERROMPU apres {prog.fait} / {len(taches)} taches.")
        print(f"  Relancer la meme commande reprend a cet endroit.")
    if rapport_orbites:
        print(f"\n  Repli par orbite : {rapport_orbites['candidats']} candidats "
              f"-> {rapport_orbites['taches']} taches "
              f"(facteur {rapport_orbites['facteur']:.2f})")
        print(f"    {rapport_orbites['cicys_sans_symetrie']} CICYs sans "
              f"symetrie de configuration : rien n'y est replie.")
        print(f"    Controle : {n_controles} membres non representants "
              f"reevalues pour de vrai, {n_discordances} discordance(s).")
        if n_discordances:
            print(f"    /!\\ Le repli est INVALIDE sur ce lot. Relancer sans "
                  f"--replier-orbites.")
        elif not n_controles:
            print(f"    /!\\ AUCUN controle effectue : le repli n'est pas "
                  f"verifie sur ce lot (--controle-orbites 0 ?).")
    print(f"\n  Couples (candidat, lambda) qui survivent   : {survivants}")
    print(f"  Indetermines (un test non calculable)     : {indetermines}"
          f"   <- ni retenus ni elimines")
    print(f"  Ecartes avant evaluation                  : {ecartes}"
          f"   <- hors domaine, charges non permutees, etc.")
    print(f"  Toutes ces lignes sont dans le JSONL, champ 'etat'.")
    print(f"  Ecrit : {dst}")
    print(f"{'=' * 96}\n")


if __name__ == '__main__':
    sys.exit(main())
