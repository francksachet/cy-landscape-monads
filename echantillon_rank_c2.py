#!/usr/bin/env python3
"""
echantillon_rank_c2.py -- Ce que couterait de decider les candidats sans
verdict, mesure sur un echantillon plutot qu'estime.

POURQUOI CE SCRIPT EXISTE
-------------------------
`scan_wilson5` porte, en candidats distincts (representants d'orbite) :

    rank_C = 1, rang_V = 3 :   472 candidats,   472 indetermines,   0 survivant
    rank_C = 1, rang_V = 4 :   968 candidats,     0 indetermine,  691 survivants
    rank_C = 2, rang_V = 3 : 1 131 candidats, 1 000 indetermines,   0 survivant

Les 691 survivants du catalogue sortent d'une seule strate sur trois. Les
1 000 de rank_C = 2 etaient bloques par les deux gardes `len(c) == 1` levees
au 5.36 ; les 472 de rank_C = 1 / rang_V = 3 le sont pour une AUTRE raison,
qui n'est pas identifiee. Ce script mesure ce que couterait de les decider,
et ce que ca donnerait.

CE QU'IL DECLARE
----------------
L'echantillon est STRATIFIE et tire avec une graine fixe ; le script affiche
la population de chaque strate et combien il en a pris. Il mesure surtout le
RAPPORT entre strates : c'est lui qui transfere d'une machine a l'autre, pas
les secondes.

DEUX LECONS DE LA VERSION PRECEDENTE, CORRIGEES ICI
---------------------------------------------------
1. Le plafond de temps etait teste ENTRE deux lambda, donc il ne protegeait
   pas de ce qui coute vraiment : la construction de l'anneau covariant, qui
   sur les grosses CICYs prend des minutes. Le cout est desormais MESURE et
   rapporte separement (`montage_s`), et l'ordre de traitement ne cache plus
   les cas lourds derriere les cas legers.
2. Rien n'etait ecrit avant la fin. Chaque candidat est maintenant ecrit
   IMMEDIATEMENT dans un JSONL, flush + fsync : une interruption ne coute que
   le candidat en cours. C'est la lecon du §5.24, appliquee a un script de
   mesure.

Usage (PowerShell, 7 coeurs) :
    python -u echantillon_rank_c2.py cicyquotients.m cicylist.txt -j 7 `
           --par-strate 40 | Tee-Object -FilePath echantillon.log

    Ctrl-C quand vous en avez assez : le JSONL est complet a la ligne pres,
    et --resume repart de ce qui manque.
"""
import os
import sys
import json
import time
import argparse
import multiprocessing as mp
from collections import defaultdict, Counter

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

GRAINE = 0
_CTX = None


def _init_worker(braun_m, cicylist):
    """Charge une fois par worker ce qui est lourd et immuable (§5.26)."""
    global _CTX
    from cy_landscape.core.braun_symmetry import parse_symmetries
    from cy_landscape.data.parse_oxford import load_oxford_file
    from wilson_match import parse_braun, parse_cicylist, apparier
    entries = {e['num']: e for e in load_oxford_file(cicylist)}
    corr, _, _ = apparier(parse_braun(braun_m), parse_cicylist(cicylist))
    _CTX = {'entries': entries,
            'inv': {v: k for k, v in corr.items()},
            'SYM': parse_symmetries(braun_m)}


def _travail(item):
    """Un candidat -> une mesure. Les exceptions reviennent comme un etat."""
    from cy_landscape.core.braun_symmetry import ordres_rt, matrice_mod_p
    from cy_landscape.core.gamma_action import (choisir_premier,
                                                racine_primitive)
    from cy_landscape.core.covariant_ring import (resoudre_covariants,
                                                  tirer_covariants,
                                                  CovariantRing,
                                                  verifier_covariance)
    from cy_landscape.core.equivariant_monad import (espace_f_equivariant,
                                                     hoppe_sur_espace,
                                                     f_sans_point_base)
    rc, rv, cicy, bs, cs, groupe = item
    fiche = {'cicy': cicy, 'rank_C': rc, 'rang_V': rv, 'groupe': groupe,
             'b': bs, 'c': cs}
    try:
        e = _CTX['entries'].get(cicy)
        nb = _CTX['inv'].get(cicy)
        if e is None or nb not in _CTX['SYM']:
            return dict(fiche, etat='cicy absente')
        b, c = json.loads(bs), json.loads(cs)
        amb, cfg = e['ambient'], np.asarray(e['config'])
        for sym in _CTX['SYM'][nb]['symetries']:
            if sym['nom'] != groupe:
                continue
            t0 = time.time()
            ordres = ordres_rt(sym['coord']) | ordres_rt(sym['poly']) | {2}
            p, _ = choisir_premier(sorted(ordres), minimum=30011)
            rac = {n: racine_primitive(p, n) for n in ordres}
            Mc = [matrice_mod_p(x, p, rac) for x in sym['coord']]
            Np = [matrice_mod_p(x, p, rac) for x in sym['poly']]
            res = resoudre_covariants(amb, cfg, Mc, Np, p)
            if res is None:
                continue
            v = res['par_convention'].get('N')
            if v is None or not v['non_degenere']:
                continue
            co = tirer_covariants(v['base'], res['offsets'], res['dims'], p,
                                  np.random.RandomState(GRAINE))
            ok, _ = verifier_covariance(amb, cfg, Mc, Np, 'N', co, p)
            if not ok:
                continue
            anneau = CovariantRing(amb, cfg, co, p)
            out = espace_f_equivariant(anneau, amb, b, c, Mc, p)
            # LE COUT EST ICI, pas dans le chemin wedge : sur les grosses
            # CICYs `resoudre_covariants` + `CovariantRing` prennent des
            # minutes. On le mesure separement au lieu de le noyer dans un
            # total qui laisserait croire que c'est le critere qui coute.
            t_montage = time.time() - t0
            if out['etat'] != 'ok' or not out['solutions']:
                return dict(fiche, etat=out['etat'], montage_s=t_montage,
                            n_lambda=0)
            mes = []
            for s in out['solutions']:
                t1 = time.time()
                hop = hoppe_sur_espace(anneau, b, c, s['base'],
                                       out['offsets'], out['dims'],
                                       out['degres'], p,
                                       np.random.RandomState(GRAINE + 5))
                t_h = time.time() - t1
                t2 = time.time()
                surj = f_sans_point_base(anneau, b, c, s['base'],
                                         out['offsets'], out['dims'],
                                         out['degres'], p,
                                         np.random.RandomState(GRAINE + 7),
                                         n_essais=2, n_degres=8)
                t_s = time.time() - t2
                mes.append({'hoppe_s': t_h, 'surj_s': t_s,
                            'hoppe': hop['stable'],
                            'hoppe_motif': hop['motif'],
                            'hoppe_valeurs': {str(k): v
                                              for k, v in hop['valeurs'].items()},
                            'surjectif': bool(surj['certifie'])})
            return dict(fiche, etat='ok', montage_s=t_montage,
                        n_lambda=len(mes), mesures=mes)
        return dict(fiche, etat='groupe introuvable')
    except Exception as exc:
        return dict(fiche, etat=f'erreur ({type(exc).__name__}: {exc})')


def collecter(chemin):
    vus = {}
    with open(chemin, encoding='utf-8') as f:
        for L in f:
            if not L.strip():
                continue
            d = json.loads(L)
            if d.get('etat') != 'ok' or d.get('verdict_replique'):
                continue
            cle = (d['cicy'], json.dumps(d['b_charges']),
                   json.dumps(d['c_charges']), d.get('groupe'))
            if cle not in vus:
                vus[cle] = (len(d.get('c_charges') or []), d.get('rang_V'),
                            bool(d.get('indetermine')), bool(d.get('survit')))
    par = defaultdict(list)
    for cle, (rc, rv, ind, surv) in vus.items():
        par[(rc, rv)].append((cle, ind, surv))
    return par


def resume(fiches):
    print(f"\n{'=' * 74}")
    print("  COUT MESURE")
    print(f"{'=' * 74}")
    ok = [f for f in fiches if f.get('etat') == 'ok' and f.get('n_lambda')]
    par_rc = defaultdict(list)
    for f in ok:
        par_rc[(f['rank_C'], f['rang_V'])].append(f)
    cout = {}
    for cle in sorted(par_rc):
        lot = par_rc[cle]
        nl = sum(f['n_lambda'] for f in lot)
        th = sum(m['hoppe_s'] for f in lot for m in f['mesures'])
        ts = sum(m['surj_s'] for f in lot for m in f['mesures'])
        tm = sum(f['montage_s'] for f in lot)
        cout[cle] = (th + ts) / max(1, nl)
        mont = sorted(f['montage_s'] for f in lot)
        print(f"  rank_C={cle[0]} rang_V={cle[1]} : {len(lot)} candidats, "
              f"{nl} lambda")
        print(f"      critere   : {th + ts:7.0f} s au total, "
              f"{cout[cle]:6.2f} s par lambda")
        print(f"      montage   : {tm:7.0f} s au total, "
              f"median {mont[len(mont) // 2]:6.1f} s, "
              f"max {mont[-1]:6.1f} s   <- c'est ici que ca coute")
    if (2, 3) in cout and (1, 4) in cout:
        print(f"\n  RAPPORT du critere, rank_C=2 / rank_C=1 : "
              f"{cout[(2, 3)] / max(1e-9, cout[(1, 4)]):.2f}")
    for cle in sorted(par_rc):
        v = Counter(str(m['hoppe']) for f in par_rc[cle] for m in f['mesures'])
        s = Counter(bool(m['surjectif'])
                    for f in par_rc[cle] for m in f['mesures'])
        print(f"\n  VERDICTS rank_C={cle[0]} rang_V={cle[1]}")
        print(f"    Hoppe : {dict(v)}   "
              f"(True = stable, False = elimine, None = non calculable)")
        print(f"    surjectivite certifiee : {dict(s)}")
        if str(None) in v:
            motifs = Counter(m['hoppe_motif'] for f in par_rc[cle]
                             for m in f['mesures'] if m['hoppe'] is None)
            for mo, n in motifs.most_common(3):
                print(f"      motif : {mo}  x{n}")
    rates = [f for f in fiches if f.get('etat') != 'ok']
    if rates:
        print(f"\n  {len(rates)} candidats non evalues : "
              f"{dict(Counter(f['etat'] for f in rates).most_common(5))}")


def principal():
    ap = argparse.ArgumentParser()
    ap.add_argument('braun_m')
    ap.add_argument('cicylist')
    ap.add_argument('--dossier', default='scan_wilson5')
    ap.add_argument('--par-strate', type=int, default=40)
    ap.add_argument('-j', '--jobs', type=int, default=None)
    ap.add_argument('--sortie', default='echantillon_rank_c2.jsonl')
    ap.add_argument('--resume', action='store_true',
                    help="repartir de ce qui manque dans --sortie")
    args = ap.parse_args()

    src = os.path.join(args.dossier, 'results_equivariance_f.jsonl')
    print(f"  lecture de {src} ...", flush=True)
    strates = collecter(src)

    print(f"\n{'=' * 74}")
    print("  POPULATION (candidats distincts, representants d'orbite seuls)")
    print(f"{'=' * 74}")
    print(f"  {'rank_C':>7} {'rang_V':>7} {'total':>8} {'indetermines':>13} "
          f"{'survivants':>11}")
    for cle in sorted(strates, key=lambda x: (x[0], x[1] or 0)):
        lst = strates[cle]
        print(f"  {cle[0]:>7} {str(cle[1]):>7} {len(lst):>8} "
              f"{sum(1 for _, i, _ in lst if i):>13} "
              f"{sum(1 for _, _, s in lst if s):>11}")

    # ECHANTILLON. Les strates SANS verdict sont l'objet de la mesure ; une
    # strate deja decidee sert d'etalon, et c'est le rapport qui transfere.
    rng = np.random.RandomState(11)
    choisis = []
    for cle in sorted(strates, key=lambda x: (x[0], x[1] or 0)):
        lst = strates[cle]
        ind = [x for x in lst if x[1]]
        pool = ind if ind else [x for x in lst if not x[1]]
        if not pool:
            continue
        for i in rng.permutation(len(pool))[:args.par_strate]:
            c = pool[int(i)][0]
            choisis.append((cle[0], cle[1], c[0], c[1], c[2], c[3]))
    print(f"\n  echantillon : {len(choisis)} candidats "
          f"({args.par_strate} par strate au plus, graine 11)")
    print("    strates sans verdict : on tire dans les INDETERMINES ;")
    print("    strate deja decidee  : on tire dans les DECIDES, comme etalon.")

    deja = {}
    if args.resume and os.path.exists(args.sortie):
        with open(args.sortie, encoding='utf-8') as f:
            for L in f:
                if L.strip():
                    d = json.loads(L)
                    deja[(d['cicy'], d['b'], d['c'], d['groupe'])] = d
        print(f"  --resume : {len(deja)} deja mesures")
    a_faire = [x for x in choisis
               if (x[2], x[3], x[4], x[5]) not in deja]
    print(f"  {len(a_faire)} a mesurer\n", flush=True)

    n_jobs = args.jobs or max(1, (os.cpu_count() or 2) - 1)
    fiches = list(deja.values())
    t0 = time.time()
    interrompu = False
    # ECRITURE IMMEDIATE. Une interruption ne coute que le candidat en cours.
    with open(args.sortie, 'a', encoding='utf-8') as fh:
        if n_jobs <= 1 or not a_faire:
            _init_worker(args.braun_m, args.cicylist)
            flux = (_travail(x) for x in a_faire)
            pool = None
        else:
            pool = mp.Pool(n_jobs, initializer=_init_worker,
                           initargs=(args.braun_m, args.cicylist))
            flux = pool.imap_unordered(_travail, a_faire)
        try:
            for n, d in enumerate(flux, 1):
                fh.write(json.dumps(d, default=str) + '\n')
                fh.flush()
                os.fsync(fh.fileno())
                fiches.append(d)
                if d.get('etat') == 'ok' and d.get('n_lambda'):
                    v = Counter(str(m['hoppe']) for m in d['mesures'])
                    print(f"  [{n}/{len(a_faire)}] #{d['cicy']:<5} "
                          f"rC={d['rank_C']} rV={d['rang_V']} "
                          f"{str(d['groupe']):<9} {d['n_lambda']} lambda   "
                          f"montage {d['montage_s']:6.1f} s   "
                          f"critere {sum(m['hoppe_s'] + m['surj_s'] for m in d['mesures']):6.1f} s"
                          f"   -> {dict(v)}", flush=True)
                else:
                    print(f"  [{n}/{len(a_faire)}] #{d['cicy']:<5} "
                          f"rC={d['rank_C']} : {d.get('etat')}", flush=True)
        except KeyboardInterrupt:
            interrompu = True
            if pool is not None:
                pool.terminate()
        finally:
            if pool is not None:
                pool.close()
                pool.join()

    if interrompu:
        print(f"\n  INTERROMPU. {len(fiches)} mesures dans {args.sortie} ; "
              f"--resume reprend le reste.")
    print(f"  {time.time() - t0:.0f} s a {n_jobs} coeurs")
    resume(fiches)
    print(f"\n  mesures : {args.sortie}")


if __name__ == '__main__':
    mp.freeze_support()
    principal()
