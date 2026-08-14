#!/usr/bin/env python3
"""
Pipeline optimise avec reprise sur checkpoint LEGER.

Architecture du checkpoint (2 fichiers separes) :
  - progress.pkl   : leger (numeros de CICYs traitees, temps, parametres).
                     Charge integralement a chaque demarrage -> toujours rapide.
  - results.jsonl  : lourd, en ecriture "append-only" (une ligne JSON par resultat).
                     JAMAIS relu pendant le scan. Uniquement lu a l'export final.

Distribution des taches : chunksize=1 dans imap_unordered, pour que chaque
CICY soit renvoyee individuellement des qu'elle est terminee. Avec un
chunksize plus grand, un worker ne renvoie ses resultats qu'une fois TOUT
son lot termine -- une seule CICY lente dans un lot de 4 bloque
l'affichage de tout le lot, meme si les 3 autres sont deja finies.

Un "heartbeat" avertit toutes les 30s si aucun resultat n'arrive, pour
distinguer un scan simplement lent d'un scan bloque.

Usage:
  python -m cy_landscape.main_optimized cicylist.txt --max-ps 6 -j 8
"""
import os, sys, json, argparse, time, signal, pickle
import multiprocessing as mp
import numpy as np
from multiprocessing import Pool, cpu_count

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def process_cicy(args):
    """Worker : traite une CICY."""
    (c, max_charge, n_random, seed, sampling_threshold, sampling_budget,
     with_extensions, cibles_gen, ordres_gamma, positive_only,
     exhaustif_max, ext_exhaustif_max) = args

    from cy_landscape.core.intersection import (
        compute_intersection_numbers, compute_euler_from_intersection, compute_c2_tangent)
    from cy_landscape.core.bundles import CICYGeometry, GAUGE_GROUP_TABLE
    from cy_landscape.core.monads import (
        MonadBundle, compute_monad_cohomology, compute_monad_cohomology_ex,
        generate_monads, check_map_exists, check_monad_nondegenerate)
    from cy_landscape.core.positive_monads import (
        generate_positive_monads, GENERATOR_VERSION)
    from cy_landscape.core.chi_exact import ChiCalculator
    from cy_landscape.core.exact_cohomology import koszul_cohomology_ex
    from cy_landscape.core.extensions import (
        ExtensionBundle, check_extension_exists, chi_extension,
        hoppe_extension, cohomology_extension_ex, generate_extensions,
        pente_extension, ContextePente)
    from cy_landscape.core.cache import set_geometry
    from cy_landscape.core.hoppe_fast import hoppe_fast, vecteur_D
    from cy_landscape.core.monad_wedge import cohomology_wedge2_V
    from cy_landscape.core.cohomology import (
        extract_spectrum_su5, extract_spectrum_so10, extract_spectrum_e6)

    rng = np.random.RandomState(seed + c['num'])
    m = len(c['ambient'])

    # 'auto' = 2 x m (m = nombre de facteurs projectifs de l'espace ambiant),
    # resolu ici car m depend de la CICY courante.
    if sampling_threshold == 'auto':
        thr = 2 * m
    elif sampling_threshold in (None, 'off'):
        thr = None
    else:
        thr = int(sampling_threshold)

    try:
        d = compute_intersection_numbers(c['ambient'], c['config'])
        chi = compute_euler_from_intersection(c['ambient'], c['config'], d)
        if chi != c['chi']:
            # chi recalcule != chi du fichier Oxford : geometrie non validee,
            # la CICY est ecartee. Tracee explicitement pour ne pas la
            # confondre avec une CICY simplement sans resultat.
            return {'cicy': c['num'], 'results': [], 'skipped': 'chi_mismatch'}
        c2 = compute_c2_tangent(c['ambient'], c['config'], d)
        chical = ChiCalculator(c['ambient'], d, c2)
        geom = CICYGeometry(
            ambient_dims=c['ambient'], config_matrix=c['config'],
            h11=c['h11'], h21=c['h21'],
            intersection_numbers=d, c2_tangent=c2)
    except Exception as e:
        return {'cicy': c['num'], 'results': [],
                'skipped': f'geometry_error:{type(e).__name__}'}

    set_geometry(c['ambient'], c['config'])
    results = []
    # [monades vues, passant le prefiltre chi, ecartees faute de certification]
    n_prefiltre = [0, 0, 0]
    n_degeneres = [0]
    # [extensions vues, passant le prefiltre chi, non eliminees par la
    #  pente, Hoppe-stables]
    n_ext = [0, 0, 0, 0]
    # Classes de Kahler candidates : ne dependent que de la geometrie, donc
    # construites une fois par CICY et non par extension.
    ctx_pente = None
    D_kahler = None
    seen = set()
    stats_gen = {}

    for rank_V in [3, 4, 5]:
        gauge = GAUGE_GROUP_TABLE.get(rank_V, {}).get("group", "?")
        all_monads = []

        # `--positive-only` restreint au generateur positif. Ce n'est pas une
        # optimisation : les monades classiques produisent des charges b_i a
        # composantes NEGATIVES, pour lesquelles le modele R_a = S_a / I_a
        # n'a aucun monome, donc `sections.domaine_valide` les rejette et
        # toute la chaine d'equivariance (equivariance_f.py) les laisse
        # « hors domaine » -- ni retenues ni eliminees. Sur scan_wilson2,
        # 37 candidats sur 108 sont dans ce cas, dont 23 des 24 qui portent
        # un Gamma d'ordre >= 4. Se restreindre au generateur positif produit
        # donc des candidats TESTABLES de bout en bout.
        generateurs = [
            ('pos_monad', generate_positive_monads(
                m, rank_V, max_charge=max_charge, n_systematic=n_random, rng=rng,
                sampling_threshold=thr, sampling_budget=sampling_budget,
                seed=seed, exhaustif_max=exhaustif_max, stats=stats_gen)),
        ]
        if not positive_only:
            generateurs.append(
                ('monad', generate_monads(m, rank_V, max_charge=3,
                                          n_random=n_random, rng=rng)))
        for kind, monads in generateurs:
            for monad in monads:
                if not monad.c1_vanishes:
                    continue
                sig = (kind, tuple(tuple(b) for b in monad.b_charges),
                       tuple(tuple(cc) for cc in monad.c_charges))
                if sig in seen:
                    continue
                seen.add(sig)
                all_monads.append((kind, monad))

        # ------------------------------------------------------------------
        # BRANCHE EXTENSION -- chemin propre (§5.10)
        # ------------------------------------------------------------------
        # L'ancien code construisait ici une pseudo-monade B = F1 (+) F2,
        # C = F2 pour reutiliser le chemin des monades. Or le noyau de
        # F1(+)F2 -> F2 est de rang rank(F1), tandis que le fibre
        # d'extension 0 -> F1 -> V -> F2 -> 0 est de rang
        # rank(F1) + rank(F2) : cohomologie, Hoppe et groupe de jauge
        # portaient sur un AUTRE objet que celui inscrit dans le resultat
        # (1571 entrees sur 1571 en incoherence de rang sur test_v3).
        #
        # Ce chemin-ci n'utilise plus la pseudo-monade. Il ne passe pas non
        # plus par la boucle `all_monads` ci-dessous : un fibre d'extension
        # n'est pas un noyau de monade et rien de ce qui suit ne
        # s'appliquerait. Ordre des filtres, du moins cher au plus cher :
        #
        #   1. chi(V) = chi(F1) + chi(F2), arithmetique pure, EXACT ;
        #   2. pente : les sous-faisceaux lisibles sur la filtration
        #      doivent pouvoir etre de degre negatif pour une classe de
        #      Kahler. On n'ecarte QUE sur certificat -- l'echec de la
        #      recherche d'un temoin ne demontre rien (voir extensions.py) ;
        #   3. Ext^1(F2, F1) != 0 -- sinon seule l'extension scindee
        #      existe, et une somme directe n'est jamais stable ;
        #   4. Hoppe par borne superieure sur les quotients gradues -- la
        #      borne ne donne pas de faux positif, mais le critere lui-meme
        #      n'est qu'une condition NECESSAIRE des que Pic(X) n'est pas
        #      de rang 1 : d'ou l'etape 2, qui voit ce que Hoppe ne voit
        #      pas ;
        #   5. cohomologie par bornes rigoureuses, sans hypothese de rang
        #      maximal -- l'hypothese qui avait fausse §4.3 et §4.4.
        stats_ext = {}
        if with_extensions:
            for ext in generate_extensions(
                    m, rank_V, max_charge=max_charge, n_random=n_random,
                    seed=seed, exhaustif_max=ext_exhaustif_max,
                    stats=stats_ext):
                try:
                    n_ext[0] += 1
                    # 1. prefiltre chi -- exact, quelques multiplications
                    chi_V = chi_extension(ext, chical)
                    if abs(chi_V) not in cibles_gen:
                        continue
                    n_ext[1] += 1

                    # 2. pente. Arithmetique pure sur les d_ijk, donc placee
                    #    avant toute cohomologie. On n'ecarte que sur
                    #    `False`, qui est DEMONTRE ; `None` est conserve et
                    #    trace, parce qu'un echec de recherche de temoin
                    #    suit le budget de recherche et non la geometrie.
                    if ctx_pente is None:
                        ctx_pente = ContextePente(d, m)
                    pen = pente_extension(ext, ctx=ctx_pente)
                    if pen['stable_possible'] is False:
                        continue
                    n_ext[2] += 1

                    # 3. l'extension non scindee doit exister
                    exists, h1_ext = check_extension_exists(
                        ext, c['ambient'], c['config'])
                    if not exists:
                        continue

                    # 4. Hoppe. La BORNE ne donne pas de faux positif ; le
                    #    critere, lui, n'est qu'une condition necessaire des
                    #    que Pic(X) n'est pas de rang 1. `stable` vaut True
                    #    ou None ; on n'inscrit que les True.
                    hop = hoppe_extension(ext, c['ambient'], c['config'])
                    if hop['stable'] is not True:
                        continue
                    n_ext[3] += 1

                    # 5. cohomologie par bornes ; None = non certifie, on
                    #    n'invente rien.
                    coh_ext = cohomology_extension_ex(
                        ext, c['ambient'], c['config'], chical=chical)
                    if coh_ext is None:
                        continue

                    b1 = coh_ext['bounds'][1]
                    b2 = coh_ext['bounds'][2]
                    results.append({
                        'type': 'extension', 'stable': True,
                        'cicy': c['num'], 'h11': c['h11'], 'h21': c['h21'],
                        'chi': c['chi'],
                        'ambient': "x".join(f"P{n}" for n in c['ambient']),
                        'rank_V': ext.rank_V, 'gauge': gauge,
                        # n_gen = |chi(V)| est LEGITIME ici : la stabilite est
                        # prouvee, donc h0 = h3 = 0.
                        'n_gen': int(abs(chi_V)),
                        'n_gen_amont': int(abs(chi_V)),
                        'chi_V': int(chi_V),
                        # h1 et h2 ne sont inscrits que s'ils sont DETERMINES.
                        'cohomology': ([0, coh_ext[1], coh_ext[2], 0]
                                       if coh_ext['determine'][1] else None),
                        'coh_bounds': {str(k): list(v) for k, v in
                                       coh_ext['bounds'].items()},
                        # chi confronte a la somme alternee des h^i quand
                        # les quatre degres sont certifies (deux chemins).
                        'chi_recoupe': coh_ext['chi_recoupe'],
                        'f1_charges': [list(x) for x in ext.f1_charges],
                        'f2_charges': [list(x) for x in ext.f2_charges],
                        # Pas de (B, C) : un fibre d'extension n'est pas un
                        # noyau de monade. Les laisser vides plutot que d'y
                        # remettre la pseudo-monade du defaut 4.7.
                        'b_charges': [], 'c_charges': [],
                        'hoppe': hop['etat'],
                        # `stable` ici = « non elimine par Hoppe », pas
                        # « stable » : voir extensions.py.
                        'pente_etat': pen['etat'],
                        'pente_verdict': pen['stable_possible'],
                        'pente_temoin_J': pen['temoin'],
                        'pente_J_exhaustif': pen['J_exhaustif'],
                        'hoppe_bornes': {str(k): int(v)
                                         for k, v in hop['bornes'].items()},
                        'ext1': int(h1_ext),
                        # Le spectre detaille (Higgs, exotiques) demanderait
                        # H^1(w^2 V) : non calcule sur cette branche.
                        'higgs': 0, 'higgs_certifie': False,
                        'exotics': None, 'singlets': None, 'reps': {},
                        'score': 30 + 25 + 2,
                        'gen': GENERATOR_VERSION,
                        # DEMONTRE sur le domaine, ou simple sondage : sans
                        # ce champ un resultat d'absence serait ininterpretable
                        # (§5.11).
                        'ext_mode': ('exhaustif'
                                     if stats_ext.get('ext_exhaustifs')
                                     else 'echantillonne'),
                        'ext_max_charge': int(max_charge),
                        'ordres_gamma': sorted(ordres_gamma),
                        'cibles_chi': sorted(cibles_gen),
                    })
                except Exception:
                    continue

        for kind, monad in all_monads:
            try:
                # --- PREFILTRE chi (exact, arithmetique pure) -------------
                # Pour V stable de pente nulle, h0(V) = h3(V) = 0, donc
                # n_gen = |h1 - h2| = |chi(V)|. La condition |chi(V)| = 3 est
                # NECESSAIRE et coute quelques multiplications, contre un
                # parcours de 2^K sous-ensembles de Koszul pour la cohomologie.
                # Selectivite mesuree : 0,01 % de monades retenues sur 797 027.
                # C'est placee AVANT check_map_exists, plus couteux.
                n_prefiltre[0] += 1

                # --- NON-DEGENERESCENCE (comparaisons de charges seules) ---
                # Rejette les monades scindees : un b_i egal a un c_j, une
                # colonne de f identiquement nulle, ou un rang structurel de f
                # inferieur a rank_C. Dans ces cas V est une somme directe de
                # fibres en droites, donc jamais stable -- et le critere de
                # Hoppe ne le voit pas quand les charges ont des signes
                # melangees. Place en tete : c'est le test le moins cher.
                _ok_nd, _motif_nd = check_monad_nondegenerate(monad)
                if not _ok_nd:
                    n_degeneres[0] += 1
                    continue

                chi_V = abs(chical.monad(monad.b_charges, monad.c_charges))
                if chi_V not in cibles_gen:
                    continue
                n_prefiltre[1] += 1

                if kind != 'extension':
                    map_ok, _ = check_map_exists(monad, c['ambient'], c['config'])
                    if not map_ok: continue

                # --- CERTIFICATION de la cohomologie ----------------------
                # koszul_cohomology ne donne des h^i exacts que si la suite
                # spectrale de Koszul degenere a E_1. Sans cette verification,
                # les h^i sont faux dans environ 60 % des cas. On exige que
                # tous les fibres en droites de la monade soient certifies.
                # NB : cela ne certifie pas la suite exacte longue dont on
                # deduit H^i(V) -- condition necessaire, pas suffisante.
                # On n'exige la certification que des degres 1 et 2 : le
                # comptage des generations ne depend que d'eux
                # (n_gen = |h1 - h2| pour V stable). La certification
                # globale rejetterait 3 fois plus de cas sans rien
                # apporter au critere de selection.
                _cert_ok = True
                for ch in list(monad.b_charges) + list(monad.c_charges):
                    _cd = koszul_cohomology_ex(
                        c['ambient'], c['config'], ch)['certified_by_degree']
                    if not (_cd[1] and _cd[2]):
                        _cert_ok = False
                        break
                if not _cert_ok:
                    n_prefiltre[2] += 1
                    continue

                # Cohomologie de V par BORNES RIGOUREUSES (voir
                # compute_monad_cohomology_ex). L'ancienne version supposait
                # les applications de rang maximal : sur les cas ou la
                # version rigoureuse determine h1 et h2, elle donnait une
                # autre reponse dans 100 % des cas testes. n_gen restait juste
                # (il vaut |chi(V)|), mais la repartition entre h1 et h2 --
                # donc n_anti = min(h1,h2), donc le classement des candidats
                # « propres » -- etait fausse.
                coh_ex = compute_monad_cohomology_ex(monad, c['ambient'], c['config'])
                if coh_ex is None: continue
                coh_cert = coh_ex['certified_by_degree']
                if not (coh_cert[1] and coh_cert[2]):
                    n_prefiltre[2] += 1
                    continue
                cohom = {i: coh_ex[i] for i in range(4)}
                if abs(cohom[1] - cohom[2]) not in cibles_gen: continue

                # D = D_k(J) a J = (1,...,1). Passe a hoppe_fast pour
                # activer la phase des twists, qui voit ce que H = 0 et
                # H = e_i ne voient pas (§5.15).
                if D_kahler is None:
                    D_kahler = [int(x) for x in vecteur_D(d, [1] * m)]
                hoppe = hoppe_fast(c['ambient'], c['config'], monad, max_H=1,
                                   D=(D_kahler if all(x > 0 for x in D_kahler)
                                      else None))
                if not hoppe['stable']: continue

                # --- wedge^2 V : uniquement si DETERMINE ------------------
                # monad_wedge v2 ne renvoie que les degres prouves. Un degre
                # absent signifie « inconnu », pas « zero ». L'ancienne
                # heuristique cohom[1]*(rV-1)//2 pour rank_C >= 2 est
                # supprimee : elle fabriquait un nombre sans fondement (elle
                # expliquait a elle seule les H=4 des candidats [0,3,0,0]).
                w2 = cohomology_wedge2_V(
                    c['ambient'], c['config'], monad.b_charges, monad.c_charges,
                    cohom_V={i: cohom[i] for i in range(4)})
                w2_cert = w2.get('certified') or {}
                w2V = dict(w2.get('wedge2V') or {})
                higgs_connu = bool(w2_cert.get(1))
                # extract_spectrum attend les quatre degres : on complete par
                # zero pour ne pas casser l'appel, mais `higgs_certifie` dans
                # le resultat dit si le nombre de Higgs veut dire quelque chose.
                w2V = {i: int(w2V.get(i, 0)) for i in range(4)}

                # end_V = None : h^1(End V) n'est PAS calcule. L'ancienne
                # valeur de remplissage, rank_V^2 - 1 codee en dur, produisait
                # un nombre de singlets invente -- et un nombre invente n'est
                # pas un nombre. `extract_spectrum_*` renvoie desormais
                # n_singlets = None, et le score ne lui accorde rien.
                fmt = {"V": {i: cohom[i] for i in range(4)},
                       "V_dual": {i: cohom[3-i] for i in range(4)},
                       "wedge2V": w2V, "end_V": None}

                if gauge == "SU(5)": sp = extract_spectrum_su5(fmt)
                elif gauge == "SO(10)": sp = extract_spectrum_so10(fmt)
                else: sp = extract_spectrum_e6(fmt)

                score = sp.sm_compatibility * 0.5 + 30 + 25
                if kind == 'pos_monad': score += 3
                if kind == 'extension': score += 2

                results.append({
                    'type': kind, 'stable': True,
                    'cicy': c['num'], 'h11': c['h11'], 'h21': c['h21'], 'chi': c['chi'],
                    'ambient': "x".join(f"P{n}" for n in c['ambient']),
                    'rank_V': rank_V, 'gauge': gauge,
                    'n_gen': int(sp.n_generations), 'higgs': int(sp.n_higgs_candidates),
                    # None = non calcule (§4.8). Jamais 0 : un zero se lirait
                    # comme « pas d'exotiques », donc comme une qualite.
                    'exotics': None if sp.n_exotics is None else int(sp.n_exotics),
                    'singlets': None if sp.n_singlets is None else int(sp.n_singlets),
                    'score': round(score, 1),
                    'cohomology': [int(cohom[i]) for i in range(4)],
                    'reps': {k: int(v) for k, v in sp.representations.items() if v > 0},
                    'b_charges': [list(b) for b in monad.b_charges],
                    'c_charges': [list(cc) for cc in monad.c_charges],
                    'hoppe': hoppe['reason'],
                    'gen': GENERATOR_VERSION,
                    'higgs_certifie': higgs_connu,
                    'coh_bounds': {str(k): list(v) for k, v in coh_ex['bounds'].items()},
                    'chi_V': coh_ex['chi'],
                    'n_gen_amont': abs(cohom[1] - cohom[2]),
                    # ordres_gamma = ORDRES des groupes (|Gamma|)
                    # cibles_chi   = valeurs visees de |chi(V)| = n_gen * |Gamma|
                    # Les deux etaient confondus dans la version precedente :
                    # le champ contenait les cibles sous le nom des ordres, et
                    # l'audit remultipliait par n_gen, rejetant a tort 90 % des
                    # resultats.
                    'ordres_gamma': sorted(ordres_gamma),
                    'cibles_chi': sorted(cibles_gen),
                    'chi_wedge2V': w2.get('chi_wedge2V'),
                    'wedge2V_bounds': {str(k): list(v) for k, v in
                                       (w2.get('bounds') or {}).items()},
                })
            except Exception:
                continue

    return {'cicy': c['num'], 'results': results, 'skipped': None,
            'prefilter': tuple(n_prefiltre), 'degeneres': n_degeneres[0],
            'extensions': tuple(n_ext)}


class ProgressFile:
    """
    Fichier de progression LEGER : uniquement les numeros de CICYs traitees,
    le temps ecoule, et les parametres. Jamais les resultats eux-memes.
    Toujours rapide a charger, quel que soit le nombre de resultats trouves.
    """

    def __init__(self, path):
        self.path = path
        self.done_cicys = set()
        self.start_time = time.time()
        self.elapsed_before = 0.0
        self.params = None

    def load(self):
        if not os.path.exists(self.path):
            return False
        try:
            with open(self.path, 'rb') as f:
                data = pickle.load(f)
            self.done_cicys = data['done_cicys']
            self.elapsed_before = data.get('elapsed', 0.0)
            self.params = data.get('params', None)
            return True
        except Exception as e:
            print(f"  Progress file corrompu ({e}), on repart de zero")
            return False

    def save(self):
        elapsed = self.elapsed_before + (time.time() - self.start_time)
        data = {
            'done_cicys': self.done_cicys,
            'elapsed': elapsed,
            'saved_at': time.time(),
            'params': self.params,
        }
        tmp = self.path + '.tmp'
        for attempt in range(3):
            try:
                with open(tmp, 'wb') as f:
                    pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
                os.replace(tmp, self.path)
                return
            except (PermissionError, OSError) as e:
                if attempt < 2:
                    time.sleep(0.5)
                else:
                    print(f"  ATTENTION : sauvegarde progress.pkl impossible ({e})")

    def total_elapsed(self):
        return self.elapsed_before + (time.time() - self.start_time)


class ResultsWriter:
    """
    Fichier de resultats en mode APPEND-ONLY (JSON Lines : un objet JSON par ligne).

    Jamais relu integralement pendant le scan -- chaque nouveau resultat est
    simplement ajoute a la fin du fichier. Le fichier n'est parcouru en entier
    qu'une seule fois, a l'export final (ou par export_checkpoint.py).
    """

    def __init__(self, path):
        self.path = path
        self._fh = None

    def open(self):
        self._fh = open(self.path, 'a', encoding='utf-8')

    def append(self, result_dict):
        line = json.dumps(result_dict, cls=_NpEncoder)
        self._fh.write(line + '\n')

    def flush(self):
        if self._fh:
            self._fh.flush()
            os.fsync(self._fh.fileno())

    def close(self):
        if self._fh:
            self._fh.close()
            self._fh = None

    def read_all(self):
        """Lit tout le fichier -- utilise UNIQUEMENT a l'export final."""
        results = []
        if not os.path.exists(self.path):
            return results
        with open(self.path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    results.append(json.loads(line))
                except json.JSONDecodeError:
                    continue
        return results


class _NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.integer,)): return int(obj)
        if isinstance(obj, (np.floating,)): return float(obj)
        if isinstance(obj, np.ndarray): return obj.tolist()
        return super().default(obj)


def _param_differences(saved, current, defaults):
    """
    Compare les parametres du checkpoint et ceux du lancement courant.

    Une cle absente du checkpoint (parametre introduit apres coup) n'est
    signalee que si sa valeur courante n'est PAS la valeur par defaut :
    sinon le comportement est inchange et la reprise reste legitime.
    """
    diffs = []
    for key in sorted(set(saved) | set(current)):
        new = current.get(key, defaults.get(key))
        if key not in saved:
            if new != defaults.get(key):
                diffs.append((key, '<absent>', new))
            continue
        old = saved[key]
        if old != new:
            diffs.append((key, old, new))
    return diffs


def cle_identite(r):
    """
    Identifiant d'un candidat : (cicy, type, charges).

    Une monade est identifiee par (B, C), un fibre d'extension par
    (F1, F2). Une cle basee sur les seuls `b_charges`/`c_charges`
    replierait TOUTES les extensions d'une meme CICY sur la cle vide et
    n'en garderait qu'une -- mesure sur un scan max_ps <= 3 : 2 647
    extensions reduites a 132. Le repli est SILENCIEUX : rien ne
    distinguerait « une seule extension trouvee » de « 2 646 ecrasees ».
    """
    b = r.get('b_charges') or r.get('f1_charges') or []
    c = r.get('c_charges') or r.get('f2_charges') or []
    return (r.get('cicy'), r.get('type'),
            tuple(tuple(x) for x in b), tuple(tuple(x) for x in c))


def deduplicate_results(results):
    """Supprime les doublons eventuels."""
    seen = set()
    unique = []
    for r in results:
        key = cle_identite(r)
        if key in seen:
            continue
        seen.add(key)
        unique.append(r)
    return unique


def _ordres_pour(num, wilson):
    """Ordres |Gamma| des groupes librement agissants sur cette CICY."""
    if not wilson:
        return []
    return sorted((wilson.get(num) or {}).get('ordres') or [])


def _cibles_pour(num, wilson, n_gen):
    """
    Valeurs de |chi(V)| acceptables pour cette CICY.

    Sans mode Wilson : {n_gen}. Avec : {n_gen * |Gamma|} pour chaque ordre
    de groupe disponible sur la variete -- plusieurs groupes peuvent
    coexister, chacun donnant une cible differente.
    """
    if not wilson:
        return {n_gen}
    info = wilson.get(num) or {}
    ordres = info.get('ordres') or []
    cibles = {n_gen * o for o in ordres if o}
    return cibles or {n_gen}


def main():
    parser = argparse.ArgumentParser(
        description="CY Landscape — Scan avec checkpoint leger (progress + resultats separes)")
    parser.add_argument('file', nargs='?', default=None)
    parser.add_argument('--max-ps', type=int, default=None)
    parser.add_argument('--chi', type=int, default=None)
    parser.add_argument('--n-random', type=int, default=150)
    parser.add_argument('--max-charge', type=int, default=4)
    parser.add_argument('--output', type=str, default='output_optimized')
    parser.add_argument('-j', '--jobs', type=int, default=None)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--checkpoint-every', type=int, default=25,
                       help="Sauvegarder progress.pkl tous les N CICYs (defaut: 25)")
    parser.add_argument('--reset', action='store_true',
                       help="Ignorer le checkpoint existant et repartir de zero")
    parser.add_argument('--sampling-threshold', type=str, default='off',
                       help="Echantillonnage hybride des monades positives : "
                            "'off' (defaut, enumeration exhaustive), 'auto' "
                            "(seuil = 2*m), ou un entier (seuil sur sum(c)). "
                            "Au-dela du seuil, tirage aleatoire plafonne.")
    parser.add_argument('--sampling-budget', type=int, default=100,
                       help="Nombre max de vecteurs c tires par charge totale "
                            "au-dela du seuil (defaut: 100). Sans effet si "
                            "--sampling-threshold off.")
    parser.add_argument('--wilson', default=None,
                       help="wilson_cicys.json produit par wilson_match.py. "
                            "Restreint le scan aux CICYs a symetrie librement "
                            "agissante et vise |chi(V)| = n_gen * |Gamma| au "
                            "lieu de n_gen.")
    parser.add_argument('--n-gen', type=int, default=3,
                       help="Nombre de generations VOULU (sur le quotient si "
                            "--wilson est utilise). Defaut: 3.")
    # RENOMME. `--with-extensions` designait l'ancien chemin par
    # pseudo-monade, qui evaluait un fibre de rang different de celui qu'il
    # annoncait (defaut 4.7, 1571 entrees sur 1571 incoherentes). Ce chemin
    # n'existe plus : le laisser sous le meme nom ferait croire qu'on
    # reactive l'ancien comportement, et un ancien script le rappellerait
    # en silence. Le nouveau nom force a relire ce qu'on active.
    parser.add_argument('--extensions', action='store_true',
                       help="Active la branche `extension` sur le chemin "
                            "propre (§5.10) : chi additif exact, Ext^1 != 0, "
                            "Hoppe par bornes sur les quotients gradues "
                            "(suffisant, jamais de faux positif), cohomologie "
                            "par bornes rigoureuses. Ne construit plus de "
                            "pseudo-monade.")
    parser.add_argument('--with-extensions', action='store_true',
                       help=argparse.SUPPRESS)
    parser.add_argument('--ext-exhaustif-max', type=int, default=200000,
                       help="Plafond du nombre de tuples ORDONNES au-dela "
                            "duquel la branche extension retombe sur "
                            "l'echantillonnage. Par defaut le domaine est "
                            "ENUMERE : le tirage n'est pas monotone en "
                            "max_charge (216 pertes sur 222 en passant de 2 "
                            "a 3), donc aucun enonce d'absence n'en decoule. "
                            "0 = forcer l'echantillonnage.")
    parser.add_argument('--positive-only', action='store_true',
                       help="N'utiliser que le generateur positif. Les monades "
                            "classiques produisent des charges negatives, hors "
                            "du domaine du modele S/I : la chaine d'equivariance "
                            "les laisse indeterminees. Cette option ne produit "
                            "que des candidats testables de bout en bout.")
    parser.add_argument('--exhaustif-max', type=int, default=0,
                       help="Pour tout vecteur c dont le nombre de B "
                            "admissibles ne depasse pas cette valeur, ENUMERER "
                            "tous les B au lieu d'en tirer 50 au hasard. 0 = "
                            "desactive (comportement historique). Mesure : les "
                            "50 tirages ne couvrent que ~0,04 %% d'un espace de "
                            "3.10^5, sans saturation, donc aucun resultat "
                            "d'absence n'en decoule. L'enumeration donne un "
                            "enonce demontre sur le domaine enumere.")
    parser.add_argument('--heartbeat', type=int, default=30,
                       help="Signal de vie toutes les N secondes si rien n'arrive (defaut: 30, 0=desactive)")
    args = parser.parse_args()

    # `--with-extensions` designait l'ancien chemin par pseudo-monade, qui
    # n'existe plus. On refuse plutot que d'activer silencieusement autre
    # chose que ce que l'utilisateur croit demander.
    if args.with_extensions:
        parser.error(
            "--with-extensions n'existe plus : il activait le chemin par "
            "pseudo-monade du defaut 4.7 (rang et chi d'un AUTRE fibre). "
            "Utiliser --extensions, qui active le chemin propre du §5.10.")

    n_jobs = args.jobs or cpu_count()
    out = args.output
    os.makedirs(out, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"  CY LANDSCAPE — SCAN AVEC CHECKPOINT LEGER")
    print(f"  {n_jobs} coeur(s) | cache | Hoppe rapide")
    print(f"{'='*60}")

    if args.file:
        from cy_landscape.data.parse_oxford import load_oxford_file
        entries = load_oxford_file(args.file)
        print(f"\n  {len(entries)} CICYs chargees depuis {args.file}")
    else:
        from cy_landscape.data.oxford_cicys import get_all_oxford
        entries = get_all_oxford()
        print(f"\n  {len(entries)} CICYs embarquees")

    if args.max_ps:
        entries = [e for e in entries if len(e['ambient']) <= args.max_ps]
        print(f"  Filtre max_ps <= {args.max_ps} : {len(entries)}")
    if args.chi is not None:
        entries = [e for e in entries if e['chi'] == args.chi]
        print(f"  Filtre chi = {args.chi} : {len(entries)}")

    progress_path = os.path.join(out, 'progress.pkl')
    results_path = os.path.join(out, 'results.jsonl')
    progress = ProgressFile(progress_path)

    current_params = {
        'max_ps': args.max_ps, 'chi': args.chi,
        'n_random': args.n_random, 'max_charge': args.max_charge,
        'seed': args.seed,
        'sampling_threshold': args.sampling_threshold,
        'sampling_budget': args.sampling_budget,
        'with_extensions': args.extensions,
        'wilson': args.wilson, 'n_gen': args.n_gen,
        'positive_only': args.positive_only,
        'exhaustif_max': args.exhaustif_max,
        'ext_exhaustif_max': args.ext_exhaustif_max,
    }

    # Valeurs par defaut des cles introduites APRES la creation d'un
    # checkpoint : un checkpoint ancien ne les contient pas. On ne considere
    # donc pas leur absence comme une divergence tant que la valeur courante
    # est la valeur par defaut -- sans quoi un scan de plusieurs jours
    # deviendrait impossible a reprendre apres cette mise a jour.
    PARAM_DEFAULTS = {'sampling_threshold': 'off', 'sampling_budget': 100,
                      'with_extensions': False, 'wilson': None, 'n_gen': 3,
                      'ext_exhaustif_max': 200000}

    if args.reset:
        for p in (progress_path, results_path):
            if os.path.exists(p):
                os.remove(p)
        print(f"\n  Checkpoint reinitialise (--reset)")
    elif progress.load():
        print(f"\n  Reprise : {len(progress.done_cicys)} CICYs deja traitees")
        print(f"  Temps deja consomme : {progress.elapsed_before/60:.1f} min")

        if progress.params is not None:
            diffs = _param_differences(progress.params, current_params, PARAM_DEFAULTS)
            if diffs:
                print(f"\n  ATTENTION : les parametres actuels different du checkpoint !")
                for key, old, new in diffs:
                    print(f"    {key} : checkpoint={old!r}  actuel={new!r}")
                print(f"  Deux options :")
                print(f"    1. Relancer avec les MEMES parametres qu'avant")
                print(f"    2. Utiliser --reset pour tout recommencer")
                return
            if args.sampling_threshold != 'off' and 'sampling_threshold' not in progress.params:
                print(f"\n  ATTENTION : ce checkpoint a ete demarre AVANT l'ajout de")
                print(f"  l'echantillonnage hybride. Reprendre avec")
                print(f"  --sampling-threshold {args.sampling_threshold} melangerait des CICYs")
                print(f"  scannees exhaustivement et des CICYs echantillonnees.")
                print(f"  Utiliser --reset, ou relancer sans --sampling-threshold.")
                return

    progress.params = current_params

    # ------------------------------------------------------------------
    # Mode « lignes de Wilson »
    # ------------------------------------------------------------------
    # Une ligne de Wilson demande un quotient X/Gamma par une symetrie
    # librement agissante. Seules 195 des 7890 CICYs en possedent une
    # (Braun, JHEP 1104 (2011) 005). Sur ces varietes, l'indice se divise :
    #     n_gen(X/Gamma) = n_gen(X) / |Gamma|
    # Chercher |chi(V)| = 3 en amont donnerait donc 3/|Gamma| generations en
    # aval -- pas un entier, donc aucun fibre exploitable. Il faut viser
    # |chi(V)| = 3 * |Gamma|, valeur qui depend de la CICY.
    wilson = None
    if args.wilson:
        with open(args.wilson, encoding='utf-8') as f:
            wilson = {int(k): v for k, v in json.load(f).items()}
        avant = len(entries)
        entries = [c for c in entries if c['num'] in wilson]
        print(f"\n  Mode Wilson : {len(entries)} / {avant} CICYs retenues "
              f"(symetrie librement agissante)")
        ordres = sorted({o for v in wilson.values() for o in v.get('ordres', [])})
        print(f"  Ordres de groupe rencontres : {ordres}")
        print(f"  Cible du prefiltre : |chi(V)| = {args.n_gen} x |Gamma|")

    to_process = [c for c in entries if c['num'] not in progress.done_cicys]
    print(f"\n  A traiter : {len(to_process)} / {len(entries)} CICYs")

    writer = ResultsWriter(results_path)

    if not to_process:
        print(f"\n  Deja tout traite. Passage direct a l'export.")
    else:
        print(f"  n_random = {args.n_random}, progress sauvegarde tous les {args.checkpoint_every} CICYs")
        print(f"  Lancement...\n")

        writer.open()
        worker_args = [(c, args.max_charge, args.n_random, args.seed,
                        args.sampling_threshold, args.sampling_budget,
                        args.extensions,
                        _cibles_pour(c['num'], wilson, args.n_gen),
                        _ordres_pour(c['num'], wilson),
                        args.positive_only,
                        args.exhaustif_max,
                        args.ext_exhaustif_max)
                       for c in to_process]

        interrupted = {'flag': False}
        def handle_sigint(signum, frame):
            if interrupted['flag']:
                print("\n  Interruption forcee.")
                sys.exit(1)
            interrupted['flag'] = True
            print("\n\n  Interruption demandee, sauvegarde en cours...")
        signal.signal(signal.SIGINT, handle_sigint)

        n_done_session = 0
        skipped = {}
        prefilter_tot = [0, 0, 0]
        degeneres_tot = [0]
        ext_tot = [0, 0, 0, 0]
        pool = None
        try:
            if n_jobs == 1:
                iterator = (process_cicy(wa) for wa in worker_args)
                use_timeout = False
            else:
                pool = Pool(n_jobs)
                # chunksize=1 : chaque CICY est renvoyee individuellement des
                # qu'elle est terminee. Avec un chunksize plus grand, un
                # worker ne renvoie ses resultats qu'une fois TOUT son lot
                # fini -- une seule CICY lente bloquerait l'affichage de
                # tout le lot, meme si les autres sont deja terminees.
                iterator = pool.imap_unordered(process_cicy, worker_args, chunksize=1)
                use_timeout = True

            print(f"  En attente des premiers resultats "
                  f"(peut prendre du temps sur les CICYs les plus grandes)...")

            last_result_time = time.time()

            while True:
                try:
                    if use_timeout and args.heartbeat > 0:
                        result = iterator.next(timeout=args.heartbeat)
                    else:
                        result = next(iterator)
                except StopIteration:
                    break
                except mp.TimeoutError:
                    waited = time.time() - last_result_time
                    print(f"    ... toujours en cours ({waited:.0f}s sans nouveau "
                          f"resultat, {n_jobs} worker(s) actifs, pas bloque)")
                    if interrupted['flag']:
                        break
                    continue

                last_result_time = time.time()
                n_done_session += 1
                cicy_num = result['cicy']

                pf = result.get('prefilter')
                if pf:
                    for _i in range(3):
                        prefilter_tot[_i] += pf[_i]
                degeneres_tot[0] += result.get('degeneres', 0)
                ex = result.get('extensions')
                if ex:
                    for _i in range(len(ext_tot)):
                        ext_tot[_i] += ex[_i]

                reason = result.get('skipped')
                if reason:
                    skipped[reason] = skipped.get(reason, 0) + 1

                for r in result['results']:
                    writer.append(r)
                if result['results']:
                    print(f"  [{len(progress.done_cicys)+1:>4}/{len(entries)}] "
                          f"CICY #{cicy_num:>4}: {len(result['results'])} stables")

                progress.done_cicys.add(cicy_num)

                if n_done_session % args.checkpoint_every == 0:
                    writer.flush()
                    progress.save()
                    remain = len(to_process) - n_done_session
                    rate = n_done_session / (time.time() - progress.start_time)
                    eta_min = remain / rate / 60 if rate > 0 else 0
                    print(f"    -- progress sauvegarde "
                          f"({n_done_session}/{len(to_process)}, "
                          f"ETA {eta_min:.0f} min) --")

                if interrupted['flag']:
                    break

        except KeyboardInterrupt:
            interrupted['flag'] = True
        finally:
            if pool is not None:
                try:
                    pool.terminate()
                    pool.join()
                except Exception:
                    pass
            writer.flush()
            writer.close()

        progress.save()

        if prefilter_tot[0]:
            vus, passes, non_cert = prefilter_tot
            _cible = (f"{args.n_gen} x |Gamma|" if args.wilson else str(args.n_gen))
            print(f"\n  Prefiltre : {vus} monades vues, {passes} avec "
                  f"|chi(V)| = {_cible} "
                  f"({100.0*passes/max(1,vus):.3f} %),")
            print(f"              {non_cert} ecartees faute de cohomologie certifiee, "
                  f"{passes - non_cert} evaluees.")
            print(f"              {degeneres_tot[0]} monades degenerees rejetees "
                  f"(scindees : b_i = c_j, colonne nulle, ou rang de f < rank_C)")

        if ext_tot[0]:
            vus, passes, pente_ok, stables = ext_tot
            print(f"\n  Extensions : {vus} vues, {passes} avec |chi(V)| a la "
                  f"cible ({100.0*passes/max(1,vus):.3f} %),")
            print(f"              {passes - pente_ok} ecartees par un "
                  f"certificat d'instabilite de pente, {stables} retenues "
                  f"apres Hoppe.")
            print(f"              Rappel : `stable` = non elimine. Hoppe est "
                  f"NECESSAIRE, pas suffisant, des que Pic(X) n'est pas de "
                  f"rang 1.")

        if skipped:
            print(f"\n  CICYs ecartees pendant cette session "
                  f"({sum(skipped.values())} au total) :")
            for reason, n in sorted(skipped.items(), key=lambda kv: -kv[1]):
                print(f"    {reason:<28} {n}")

        if interrupted['flag']:
            print(f"\n  Sauvegarde effectuee : {len(progress.done_cicys)} CICYs traitees.")
            print(f"  Relance la meme commande pour reprendre.")
            print(f"  Progress : {progress_path}")
            print(f"  Resultats : {results_path}")
            return

    print(f"\n  Lecture des resultats accumules...")
    all_results = writer.read_all()
    all_results = deduplicate_results(all_results)
    all_results.sort(key=lambda r: r['score'], reverse=True)
    n_stable = len(all_results)
    dt = progress.total_elapsed()

    print(f"\n{'='*60}")
    print(f"  RESULTATS")
    print(f"{'='*60}")
    print(f"  CICYs traitees        : {len(progress.done_cicys)}")
    print(f"  Fibres Hoppe-stables  : {n_stable}")
    print(f"  Temps total (cumule)  : {dt/60:.1f} min")

    if all_results:
        print(f"\n  Top 10 :")
        for i, r in enumerate(all_results[:10]):
            print(f"  {i+1:>2} {r['type']:<12} #{r['cicy']:>4} "
                  f"({r['h11']:>2},{r['h21']:>2}) rk{r['rank_V']} {r['gauge']:>7} "
                  f"H={r['higgs']:>2} "
                  f"E={('?' if r.get('exotics') is None else r['exotics']):>1} "
                  f"score={r['score']}")

    export = {
        'parameters': vars(args),
        'n_cicys': len(progress.done_cicys),
        'n_stable': n_stable,
        'time_seconds': round(dt, 1),
        'results': all_results[:500],
    }
    json_path = os.path.join(out, 'results.json')
    with open(json_path, 'w') as f:
        json.dump(export, f, indent=2, cls=_NpEncoder)
    print(f"\n  Export : {json_path}")

    try:
        from cy_landscape.core.visualize import generate_all_plots
        print(f"\n  Generation des graphiques...")
        generate_all_plots(json_path, out)
    except Exception as e:
        print(f"  Graphiques ignores : {e}")

    print(f"{'='*60}")


if __name__ == "__main__":
    main()
