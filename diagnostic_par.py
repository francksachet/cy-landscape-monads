#!/usr/bin/env python3
"""
diagnostic_par.py -- Pourquoi aucun lot n'est rendu ?

Repond a trois questions, dans cet ordre, en quelques minutes :

  1. Combien coute UNE realisation sur CETTE machine ?
     (reference : 6,7 s sur un conteneur au repos)
  2. Le Pool multiprocessing fonctionne-t-il ici ?
     (test minimal, deux lots, avec le VRAI code des workers)
  3. Les workers demarrent-ils ?
     (le chargement du contexte est chronometre separement)

N'ecrit rien, ne touche a aucun fichier de resultat : il peut tourner
pendant le scan. Il prend un coeur, donc autant arreter le scan avant.

Usage:
    python diagnostic_par.py cicyquotients.m cicylist.txt scan_wilson4
"""
import os
import sys
import json
import time
import multiprocessing as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def principal():
    braun_m, cicylist, dossier = sys.argv[1], sys.argv[2], sys.argv[3]
    import equivariance_f as EF
    from cy_landscape.data.parse_oxford import load_oxford_file

    print(f"\n{'=' * 70}")
    print("  DIAGNOSTIC DU PARALLELISME")
    print(f"{'=' * 70}")
    print(f"  coeurs vus par Python : {os.cpu_count()}")
    print(f"  methode de demarrage  : {mp.get_start_method()}")
    for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        print(f"  {v:<22}: {os.environ.get(v)}")

    # ---- contexte ---------------------------------------------------
    t = time.time()
    EF._init_worker(braun_m, cicylist)
    t_ctx = time.time() - t
    print(f"\n  [1] Chargement du contexte : {t_ctx:.1f} s")
    print(f"      (c'est ce que chaque worker paie une fois au demarrage)")
    if t_ctx > 60:
        print(f"      /!\\ TRES LONG. Sept workers le paient en parallele,")
        print(f"          mais s'il depasse plusieurs minutes, c'est ici que")
        print(f"          le temps passe.")

    # ---- reconstitution d'un lot reel -------------------------------
    import types
    E = {e['num']: e for e in load_oxford_file(cicylist)}
    src = os.path.join(dossier, 'results_equivariant.jsonl')
    rs = [json.loads(l) for l in open(src, encoding='utf-8') if l.strip()]
    args = types.SimpleNamespace(replier_orbites=True, cicy=None,
                                 controle_orbites=20)
    taches, _ = EF._construire_taches(rs, E, args)

    # la tache 823 est celle qui bloque ; a defaut, la premiere venue
    k = 823 if len(taches) > 823 else 0
    i_rep, membres = taches[k]
    r = rs[i_rep]
    num_b = EF._CTX['inv'].get(r['cicy'])
    g = set(r.get('groupes_utiles') or [])
    idx = [n for n, sy in enumerate(EF._CTX['SYM'][num_b]['symetries'])
           if (not g) or sy['nom'] in g]
    print(f"\n  Tache {k} : CICY #{r['cicy']}, {len(idx)} realisations de "
          f"{sorted(g)}")

    # ---- [2] cout d'une realisation ---------------------------------
    print(f"\n  [2] Cout d'une realisation (3 mesures)")
    for n in range(3):
        item = (('D', k, n), r['cicy'], r['b_charges'], r['c_charges'],
                sorted(g) if g else None, idx[n:n + 1], n == 0)
        t = time.time()
        _, lignes = EF._travail(item)
        print(f"      realisation {n + 1} : {time.time() - t:>7.1f} s  "
              f"({len(lignes)} lignes)")
    print(f"      reference conteneur au repos : 6,7 s")

    # ---- [3] le Pool fonctionne-t-il ? ------------------------------
    print(f"\n  [3] Pool a 2 workers, 2 lots de 2 realisations")
    lots = [(('D', k, 100 + j), r['cicy'], r['b_charges'], r['c_charges'],
             sorted(g) if g else None, idx[2 * j:2 * j + 2], False)
            for j in range(2)]
    t = time.time()
    try:
        with mp.Pool(2, initializer=EF._init_worker,
                     initargs=(braun_m, cicylist)) as pool:
            n_rendus = 0
            for id_lot, lignes in pool.imap_unordered(EF._travail, lots):
                n_rendus += 1
                print(f"      lot {id_lot[2]} rendu apres "
                      f"{time.time() - t:>6.1f} s ({len(lignes)} lignes)")
        print(f"      -> {n_rendus}/2 lots rendus en {time.time() - t:.1f} s")
        if n_rendus == 2:
            print(f"      Le Pool fonctionne. Le probleme n'est pas la.")
    except Exception as exc:
        print(f"      /!\\ LE POOL A ECHOUE : {type(exc).__name__}: {exc}")
        print(f"      C'est la cause. Relancer le scan avec -j 1.")

    print(f"\n{'=' * 70}")
    print("  Envoyer cette sortie en entier.")
    print(f"{'=' * 70}\n")


if __name__ == '__main__':
    mp.freeze_support()
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)
    principal()
