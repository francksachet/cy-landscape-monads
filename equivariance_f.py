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
import time
import argparse
import hashlib
import multiprocessing as mp

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
    """
    Progression : l'ENSEMBLE des lots termines, et non plus un compteur.

    Le compteur suffisait tant que les lots etaient traites dans l'ordre.
    Avec plusieurs workers, le lot 12 peut finir avant le lot 7 : « les
    n premiers sont faits » n'a plus de sens, et le premier lot manquant
    ferait rejouer tout ce qui le suit. On enregistre donc les identifiants.

    Il n'y a plus d'offset : le JSONL est filtre a la reprise (les lignes
    portent leur `_lot`), ce qui elimine exactement les lignes des lots non
    valides -- y compris celles d'un lot ecrit a moitie.
    """

    def __init__(self, dossier, empreinte, empreinte_heritee=None):
        self.chemin = os.path.join(dossier, 'progress_equivariance_f.json')
        self.empreinte = empreinte
        # Empreinte du format SEQUENTIEL, calculee sans la taille de lot
        # (qui n'existait pas). Sans elle, passer a la version parallele
        # jetterait un checkpoint parfaitement valide : dans le cas reel,
        # 823 taches et six heures de calcul.
        self.empreinte_heritee = empreinte_heritee
        # {identifiant de lot: nombre de lignes ecrites}. Le compte permet
        # de reperer un lot ecrit a moitie (cf. `charger`).
        self.faits = {}
        self.offset_herite = None
        self.compteurs = {'survivants': 0, 'indetermines': 0, 'ecartes': 0}

    def charger(self, taille_jsonl=0):
        if not os.path.exists(self.chemin):
            return False, "aucun checkpoint"
        try:
            with open(self.chemin, encoding='utf-8') as f:
                d = json.load(f)
        except (json.JSONDecodeError, OSError) as e:
            return False, f"checkpoint illisible ({type(e).__name__})"
        herite = (self.empreinte_heritee is not None
                  and d.get('empreinte') == self.empreinte_heritee
                  and 'lots' not in d)
        if d.get('empreinte') != self.empreinte and not herite:
            return False, ("le fichier d'entree a change depuis le checkpoint "
                           "(ou --cicy / les options de repli / la taille de "
                           "lot different) -- reprise refusee")
        self.compteurs.update(d.get('compteurs') or {})
        if 'lots' in d:
            # [identifiant, nombre de lignes ecrites]. Le compte est ce qui
            # permet de detecter un lot ecrit A MOITIE : sans lui, un JSONL
            # tronque en plein milieu d'un lot laissait ce lot marque
            # « fait » avec la moitie de ses lignes, definitivement.
            self.faits = {tuple(x[0]): x[1] for x in d['lots']}
            return True, f"reprise sur {len(self.faits)} lots"
        # MIGRATION d'un checkpoint sequentiel (`fait` = nombre de taches
        # terminees dans l'ordre). Sans elle, passer a la version parallele
        # jetterait le travail deja fait -- 823 taches, six heures dans le
        # cas reel qui a motive ce changement.
        if 'fait' in d:
            # None = compte inconnu (l'ancien format ne l'enregistrait pas) ;
            # la validation par l'offset le remplace.
            self.faits = {('T', k): None for k in range(int(d['fait']))}
            # L'ancien format ecrivait dans l'ordre et garantissait que tout
            # ce qui precede `offset` appartient a une tache terminee. On
            # tronque donc la, comme le faisait la version sequentielle ;
            # au-dela il ne peut y avoir que les restes d'une tache
            # interrompue.
            self.offset_herite = int(d.get('offset', 0))
            return True, (f"reprise d'un checkpoint sequentiel : "
                          f"{int(d['fait'])} taches")
        return False, "checkpoint sans lots ni compteur"

    def sauver(self, compteurs):
        self.compteurs = dict(compteurs)
        tmp = self.chemin + '.tmp'
        with open(tmp, 'w', encoding='utf-8') as f:
            json.dump({'empreinte': self.empreinte,
                       'lots': [[list(k), v]
                                for k, v in sorted(self.faits.items())],
                       'compteurs': self.compteurs}, f)
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



# ======================================================================
# Parallelisme : un LOT = (tache, tranche de realisations de symetrie)
# ======================================================================
#
# POURQUOI PAS « UNE TACHE = UN LOT »
# -----------------------------------
# Une tache, ce n'est pas un calcul : c'est un candidat confronte a TOUTES
# les realisations que Braun donne de son groupe. #480 en a 368 pour
# Z2 x Z2 -- mesure : 24,1 s chacune, soit 2 h 28 pour la seule tache 823.
# Prendre la tache comme unite laisserait un worker occupe deux heures et
# demie pendant que les autres finissent, et surtout ne permettrait pas de
# sauvegarder au milieu. On decoupe donc en TRANCHES de realisations.
#
# Consequences, toutes voulues :
#   - repartition : les grosses taches se repartissent sur les workers ;
#   - checkpoint : la granularite descend de 2 h 28 a quelques minutes ;
#   - affichage : chaque lot rendu produit ses lignes, donc l'ecran vit.
#
# Mesure sur la machine de reference : 12 % de CPU sur 8 coeurs, autrement
# dit UN seul coeur occupe. `main_optimized` distribue depuis toujours ;
# ce script-ci ne l'avait jamais fait, parce qu'il tournait en une heure
# sur 108 candidats. Sur 3 698 taches, cet oubli vaut un facteur 7.

_CTX = {}


def _init_worker(braun_m, cicylist):
    """Charge une fois par worker ce qui est lourd et immuable.

    Sous Windows (methode `spawn`), rien n'est herite du parent : chaque
    worker relit les 7 890 CICYs et les symetries de Braun. ~4 s une fois,
    contre autant par tache si on le refaisait a chaque lot.
    """
    global _CTX
    entries = {e['num']: e for e in load_oxford_file(cicylist)}
    braun = parse_braun(braun_m)
    cl = parse_cicylist(cicylist)
    corr, _, _ = apparier(braun, cl)
    _CTX = {'entries': entries,
            'inv': {v: k for k, v in corr.items()},
            'SYM': parse_symmetries(braun_m)}


def _travail(item):
    """Evalue une tranche de realisations. Rend (id_lot, lignes).

    Les exceptions sont capturees et renvoyees comme une ligne d'etat : un
    worker qui meurt en silence ferait disparaitre un candidat sans laisser
    de trace -- le defaut du §5.23, sous une autre forme.
    """
    id_lot, cicy, b, c, groupes, idx_sym, premier = item
    try:
        e = _CTX['entries'].get(cicy)
        num_b = _CTX['inv'].get(cicy)
        if e is None or num_b not in _CTX['SYM']:
            return id_lot, []
        amb, cfg = e['ambient'], np.asarray(e['config'])
        if idx_sym is None:
            # Aucun groupe d'ordre compatible avec l'indice : le candidat
            # est ECARTE sans calcul, et la raison persistee. Ce cas est
            # decide en amont, AVANT le test de domaine -- l'inverse
            # etiquetterait « hors domaine » des candidats qui n'ont
            # simplement pas de Gamma utilisable (18 cas sur le lot de
            # #6947 lors du premier essai).
            return id_lot, [{'groupe': '-',
                             'etat': 'aucun groupe d ordre compatible avec '
                                     'l indice (|chi| != 3.|Gamma|)'}]
        if not domaine_valide(amb, cfg, b, c, rank_c_max=None):
            # Une seule fois par tache, pas une fois par tranche.
            return id_lot, ([{'groupe': '-',
                              'etat': 'hors domaine (modele S/I non valide)'}]
                            if premier else [])
        syms = [_CTX['SYM'][num_b]['symetries'][i] for i in idx_sym]
        return id_lot, analyser(cicy, amb, cfg, b, c, syms,
                                groupes=set(groupes) if groupes else None)
    except Exception as exc:
        return id_lot, [{'groupe': '-',
                         'etat': f'erreur worker ({type(exc).__name__}: {exc})'}]


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
    ap.add_argument('-j', '--jobs', type=int, default=None,
                    help="Nombre de workers. Defaut : nombre de coeurs - 1. "
                         "Ce script etait mono-coeur jusqu'ici, alors que "
                         "`main_optimized` distribue depuis toujours : sur "
                         "une machine a 8 coeurs, cet oubli valait un "
                         "facteur 7. 1 = sequentiel (utile pour deboguer).")
    ap.add_argument('--taille-lot', type=int, default=16,
                    help="Nombre de realisations de symetrie par lot. C'est "
                         "la granularite du checkpoint ET de l'affichage. "
                         "#480 a 368 realisations pour un seul candidat, a "
                         "24 s chacune : sans decoupage, un worker resterait "
                         "2 h 28 sans rien rendre ni rien sauvegarder. "
                         "Change les identifiants de lot, donc invalide un "
                         "checkpoint pris avec une autre valeur.")
    ap.add_argument('--arret-apres', type=int, default=0,
                    help="S'arreter proprement apres N lots, comme le "
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
    prog = Progression(
        args.output_dir,
        _empreinte(src, (args.cicy, bool(args.replier_orbites),
                         int(args.controle_orbites), int(args.taille_lot))),
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
            survivants = prog.compteurs['survivants']
            indetermines = prog.compteurs['indetermines']
            ecartes = prog.compteurs['ecartes']
            # FILTRAGE plutot que troncature : avec plusieurs workers les
            # lots ne finissent pas dans l'ordre, donc « tout ce qui est
            # avant l'offset est valide » est faux. Chaque ligne porte son
            # `_lot` ; on ne garde que celles des lots valides, ce qui
            # elimine aussi les lignes d'un lot ecrit a moitie.
            if prog.offset_herite is not None and taille:
                with open(dst, 'r+b') as fh:
                    fh.truncate(prog.offset_herite)
                taille = prog.offset_herite
            # LE FICHIER FAIT FOI, en DEUX passes.
            #
            # Un seul parcours ne suffit pas : il faudrait connaitre le
            # compte de lignes de chaque lot AVANT de decider si on le
            # garde. Filtrer d'abord, restreindre ensuite, laissait dans le
            # fichier les lignes d'un lot a moitie ecrit tout en le
            # recalculant -- donc en double. Constate : 41 lignes la ou le
            # run d'un trait en donne 40.
            #
            # Passe 1 : compter par lot. Passe 2 : ne reecrire que les lots
            # dont le compte correspond exactement au checkpoint.
            n_avant = n_apres = 0
            vus = {}
            if taille:
                with open(dst, encoding='utf-8') as f:
                    for ligne in f:
                        if not ligne.strip():
                            continue
                        n_avant += 1
                        try:
                            lot = json.loads(ligne).get('_lot')
                        except json.JSONDecodeError:
                            continue
                        if lot is not None:
                            k_ = tuple(lot)
                            vus[k_] = vus.get(k_, 0) + 1
            avant_restriction = len(prog.faits)
            prog.faits = {
                f: n for f, n in prog.faits.items()
                if len(f) == 2                       # herite du sequentiel
                or (f[0] == 'T' and vus.get(f, 0) == n)}
            if len(prog.faits) != avant_restriction:
                print(f"  Checkpoint restreint : {avant_restriction} -> "
                      f"{len(prog.faits)} lots (lignes absentes ou "
                      f"incompletes dans le JSONL, ou lots de controle "
                      f"-- tous seront recalcules)")
            if taille:
                tmp = dst + '.filtre'
                with open(dst, encoding='utf-8') as src_f, \
                        open(tmp, 'w', encoding='utf-8') as out_f:
                    for ligne in src_f:
                        if not ligne.strip():
                            continue
                        try:
                            lot = json.loads(ligne).get('_lot')
                        except json.JSONDecodeError:
                            continue          # ligne tronquee par une coupure
                        # Les lignes d'un checkpoint SEQUENTIEL n'ont pas de
                        # `_lot` : l'ancien format ecrivait dans l'ordre et
                        # a deja ete tronque a son offset ci-dessus.
                        if lot is None or tuple(lot) in prog.faits:
                            out_f.write(ligne)
                            n_apres += 1
                os.replace(tmp, dst)
            print(f"\n  Reprise : {motif}")
            if n_avant != n_apres:
                print(f"  JSONL filtre : {n_avant} lignes -> {n_apres} "
                      f"({n_avant - n_apres} appartenaient a des lots non "
                      f"termines)")
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

    # ---------------- lots, pool, ecriture ---------------------------
    interrompu = False
    controles = _echantillon_controle(taches, args.controle_orbites)
    n_controles = n_discordances = 0

    def _cle_verdict(x):
        # Ce qui doit coincider sur une orbite : le verdict, pas les
        # charges ni le degre temoin -- ceux-la sont PERMUTES par
        # l'automorphisme et differeraient legitimement.
        #
        # `_norm` n'est pas cosmetique. Le representant est RELU depuis le
        # JSONL (ou un tuple est devenu une liste) tandis que le membre de
        # controle est encore en memoire : `str((1, 1))` != `str([1, 1])`,
        # et les 18 controles rendaient 12 discordances entierement
        # fictives. Une garde qui crie a tort est aussi inutile qu'une
        # garde muette.
        def _norm(v):
            if isinstance(v, (list, tuple)):
                return tuple(_norm(u) for u in v)
            return v
        return (str(x.get('groupe')), str(_norm(x.get('lambda'))),
                bool(x.get('survit')), str(x.get('etat')),
                str(x.get('n_gen_quotient')))

    # --- construction des lots ---------------------------------------
    # Un lot = ('T', tache, rang de la tranche) pour un representant,
    #          ('C', tache, membre, rang) pour un controle.
    # L'identifiant est STABLE d'un lancement a l'autre : c'est lui qui
    # permet a la reprise de savoir ce qui est deja fait, quel que soit
    # l'ordre dans lequel les workers ont rendu.
    lots = []
    for k, (i_rep, membres) in enumerate(taches):
        r = rs[i_rep]
        num_b = inv.get(r['cicy'])
        if entries.get(r['cicy']) is None or num_b not in SYM:
            continue
        g = None if args.tous_groupes else set(r.get('groupes_utiles') or [])
        if g is not None and not g:
            # « aucun groupe d'ordre compatible » : pas de calcul, mais une
            # ligne quand meme -- un fichier doit dire pourquoi un cas n'a
            # pas ete traite.
            lots.append((('T', k, -1), r['cicy'], r['b_charges'],
                         r['c_charges'], None, None, True))
            continue
        idx = [n for n, sy in enumerate(SYM[num_b]['symetries'])
               if (not g) or sy['nom'] in g]
        tranches = [idx[t:t + args.taille_lot]
                    for t in range(0, len(idx), args.taille_lot)] or [[]]
        for t, tr in enumerate(tranches):
            lots.append((('T', k, t), r['cicy'], r['b_charges'],
                         r['c_charges'], sorted(g) if g else None, tr, t == 0))
        for j in controles.get(k, ()):
            rj = rs[j]
            for t, tr in enumerate(tranches):
                lots.append((('C', k, j, t), rj['cicy'], rj['b_charges'],
                             rj['c_charges'], sorted(g) if g else None,
                             tr, t == 0))

    # `('T', k)` sans rang de tranche = tache entiere d'un checkpoint
    # sequentiel migre : tous ses lots sont consideres faits.
    a_faire = [x for x in lots
               if x[0] not in prog.faits and ('T', x[0][1]) not in prog.faits]
    print(f"\n  {len(lots)} lots (tranches de {args.taille_lot} realisations), "
          f"{len(lots) - len(a_faire)} deja faits, {len(a_faire)} a traiter")
    # `--arret-apres` simule une interruption : le run doit se terminer
    # comme s'il avait ete coupe, INTERROMPU compris, sinon un test de
    # reprise ne saurait pas qu'il vient d'etre coupe.
    tronque_volontairement = False
    if args.arret_apres and len(a_faire) > args.arret_apres:
        a_faire = a_faire[:args.arret_apres]
        tronque_volontairement = True
        print(f"  --arret-apres : on s'arretera apres {len(a_faire)} lots")

    # --- execution ----------------------------------------------------
    n_jobs = args.jobs if args.jobs else max(1, (os.cpu_count() or 2) - 1)
    resultats_ctrl = {}
    t0 = time.time()
    n_faits = 0

    def _ecrire(fh, id_lot, lignes):
        """Ecrit les lignes d'un lot, repliquees sur les membres. Rend leur nombre."""
        nonlocal survivants, indetermines, ecartes
        n_ecrites = 0
        k = id_lot[1]
        i_rep, membres = taches[k]
        for L in lignes:
            if L.get('etat') != 'ok':
                ecartes += len(membres)
            else:
                survivants += bool(L.get('survit')) * len(membres)
                indetermines += bool(L.get('indetermine')) * len(membres)
                L['n_gen_quotient'] = n_gen_quotient(
                    rs[i_rep].get('cohomology'), L.get('groupe'))
        for j in membres:
            replique = (j != i_rep)
            ident_j = {c: rs[j].get(c) for c in
                       ('cicy', 'gauge', 'rank_V', 'cohomology',
                        'b_charges', 'c_charges', 'groupes_utiles',
                        'equivariant_possible', 'ordres_gamma')}
            for x in lignes:
                y = dict(x)
                y.update(ident_j)
                # Charges du membre, verdict du representant. Le degre
                # temoin reste celui du representant : il est permute par
                # l'automorphisme -- d'ou le marquage.
                y['verdict_replique'] = bool(replique)
                if replique:
                    y['representant'] = rs[i_rep].get('b_charges')
                y['_lot'] = list(id_lot)
                if x.get('etat') != 'ok':
                    y.setdefault('survit', False)
                    y.setdefault('indetermine', True)
                fh.write(json.dumps(y, default=int) + '\n')
                n_ecrites += 1
        return n_ecrites

    def _afficher(cicy, gauge, lignes):
        for L in lignes:
            if L.get('etat') != 'ok':
                print(f"  {cicy:>5} {gauge:>7} {L.get('groupe', '-'):<11} "
                      f"{L.get('etat')}")
                continue
            if L['survit']:
                ng = L.get('n_gen_quotient')
                v = (f"SURVIT -- {ng if ng is not None else '?'} gen sur "
                     f"X/Gamma (Hoppe complet + surjectif en "
                     f"{L['surjectif_degre']})")
            elif L.get('hoppe_complet') is False:
                v = f"tue par Hoppe complet : {L.get('hoppe_valeurs')}"
            elif L.get('hoppe_complet') is None and L.get('indetermine') \
                    and L.get('surjectif_certifie') is None:
                v = "indetermine : Hoppe complet non calculable"
            elif L.get('surjectif_certifie') is False:
                v = "indetermine : surjectivite de f non certifiee"
            elif L.get('indetermine'):
                v = "indetermine (w2V non calculable)"
            elif L['h0_generique'] != 0:
                v = "deja non stable avec f generique"
            elif L['h0_equivariant'] != 0:
                v = "tue par h0(V) equivariant"
            elif L.get('h0w2_generique'):
                v = "deja non stable : h0(w2V) generique != 0"
            else:
                v = "tue par h0(w2V) equivariant"
            f_ = lambda x: '-' if x is None else str(x)
            print(f"  {cicy:>5} {gauge:>7} {L['rang_V']:>2} "
                  f"{L['groupe']:<11} {str(L['lambda']):>14} "
                  f"{L['dim_equivariant']:>6} {L['dim_totale']:>5} "
                  f"{L['h0_generique']:>6} {L['h0_equivariant']:>5} "
                  f"{f_(L['h0w2_generique']):>6} "
                  f"{f_(L['h0w2_equivariant']):>5}  {v}")

    print(f"  {n_jobs} worker(s)\n")
    with open(dst, 'a', encoding='utf-8') as fh:
        if n_jobs <= 1 or not a_faire:
            _init_worker(args.braun_m, args.cicylist)
            flux = ((it[0], _travail(it)[1]) for it in a_faire)
            pool = None
        else:
            pool = mp.Pool(n_jobs, initializer=_init_worker,
                           initargs=(args.braun_m, args.cicylist))
            flux = pool.imap_unordered(_travail, a_faire)
        try:
            for id_lot, lignes in flux:
                if id_lot[0] == 'C':
                    # `n_gen_quotient` est calcule dans `_ecrire` pour les
                    # lignes du representant ; il faut le calculer AUSSI
                    # ici, sinon la comparaison oppose une valeur a None et
                    # declare une discordance a chaque controle. Premier
                    # essai : 12 discordances sur 18, toutes fictives.
                    coh = rs[id_lot[2]].get('cohomology')
                    for L in lignes:
                        if L.get('etat') == 'ok':
                            L['n_gen_quotient'] = n_gen_quotient(
                                coh, L.get('groupe'))
                    resultats_ctrl.setdefault((id_lot[1], id_lot[2]),
                                              []).extend(lignes)
                    n_lignes_ecrites = 0
                else:
                    k = id_lot[1]
                    _afficher(rs[taches[k][0]]['cicy'],
                              rs[taches[k][0]].get('gauge', ''), lignes)
                    n_lignes_ecrites = _ecrire(fh, id_lot, lignes)
                    fh.flush()
                    os.fsync(fh.fileno())
                prog.faits[id_lot] = n_lignes_ecrites
                prog.sauver({'survivants': survivants,
                             'indetermines': indetermines,
                             'ecartes': ecartes})
                n_faits += 1
                if n_faits % 25 == 0:
                    ec = time.time() - t0
                    reste = (len(a_faire) - n_faits) * ec / n_faits
                    print(f"    [{n_faits}/{len(a_faire)} lots, {ec/60:.0f} min "
                          f"ecoulees, ~{reste/3600:.1f} h restantes]")
        except KeyboardInterrupt:
            interrompu = True
            if pool is not None:
                pool.terminate()
        finally:
            if pool is not None:
                pool.close()
                pool.join()

    # --- controle du repli, une fois les lots rassembles ---------------
    n_controles = len(resultats_ctrl)
    if resultats_ctrl:
        # Le representant a ete ecrit ligne par ligne ; on relit ses lignes
        # depuis le JSONL plutot que de les garder en memoire, ce qui
        # marcherait mal sur un lot de 47 Mo.
        par_tache = {}
        with open(dst, encoding='utf-8') as f:
            for ligne in f:
                if not ligne.strip():
                    continue
                x = json.loads(ligne)
                lot = x.get('_lot')
                if not lot or lot[0] != 'T' or x.get('verdict_replique'):
                    continue
                par_tache.setdefault(lot[1], []).append(x)
        n_discordances = 0
        for (k, j), bloc_j in resultats_ctrl.items():
            a = sorted(map(_cle_verdict, par_tache.get(k, [])))
            b_ = sorted(map(_cle_verdict, bloc_j))
            if a and a != b_:
                n_discordances += 1
                print(f"  !! DISCORDANCE D'ORBITE sur "
                      f"#{rs[taches[k][0]]['cicy']} : le membre de controle "
                      f"ne donne pas le meme verdict que le representant.")
                print(f"     representant : {a[:2]}")
                print(f"     membre       : {b_[:2]}")

    interrompu = interrompu or tronque_volontairement
    if interrompu:
        print(f"\n  INTERROMPU apres {n_faits} / {len(a_faire)} lots de cette "
              f"session ({len(prog.faits)} lots faits au total).")
        print(f"  Relancer la meme commande reprend a cet endroit.")
    if rapport_orbites:
        print(f"\n  Repli par orbite : {rapport_orbites['candidats']} candidats "
              f"-> {rapport_orbites['taches']} taches "
              f"(facteur {rapport_orbites['facteur']:.2f})")
        print(f"    {rapport_orbites['cicys_sans_symetrie']} CICYs sans "
              f"symetrie de configuration : rien n'y est replie.")
        print(f"    Controle : {n_controles} membres non representants "
              f"reevalues pour de vrai, {n_discordances} discordance(s).")
        print(f"    Lots : {len(prog.faits)} termines sur {len(lots)}.")
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
