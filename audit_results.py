#!/usr/bin/env python3
"""
audit_results.py -- Triage a posteriori d'un scan deja termine.

Ne recalcule RIEN : lit uniquement `results.jsonl` et applique des tests de
coherence interne sur les donnees deja stockees (cohomologie, charges,
rangs). Permet de separer les candidats exploitables des artefacts sans
relancer le moindre calcul.

Usage:
    python audit_results.py output_optimized
    python audit_results.py output_optimized --top 20

Produit dans le meme dossier :
    results_clean.jsonl    candidats passant tous les tests
    results_flagged.jsonl  candidats rejetes, avec le champ "flags"

--------------------------------------------------------------------------
Tests appliques
--------------------------------------------------------------------------

[rank_mismatch]  len(b_charges) - len(c_charges) != rank_V
    Le rang du fibre reellement evalue ne correspond pas au rang annonce,
    donc pas non plus au groupe de jauge annonce. Touche les entrees de
    type "extension" : le pipeline construit une pseudo-monade
    B = F1 (+) F2, C = F2, dont le noyau est de rang rank(F1), et non le
    fibre d'extension de rang rank(F1)+rank(F2). Cohomologie, test de
    Hoppe et groupe de jauge portent donc sur un autre objet que celui
    annonce.

[h0_nonzero]  h0(V) != 0
    Pour V stable de pente nulle (c1 = 0), une section globale non nulle
    fournirait un sous-faisceau O de meme pente, ce qui contredit la
    stabilite. h0(V) > 0 implique donc que V n'est PAS stable.

[h3_nonzero]  h3(V) != 0
    Par dualite de Serre sur un CY3 (K trivial), h3(V) = h0(V*). V* est
    stable de pente nulle si V l'est, donc h3(V) = 0. Meme conclusion.

[index_mismatch]  |h1 - h2| != 3
    Verification de coherence du filtre de selection lui-meme.

[n_gen_not_3]  n_gen != 3

[wedge2_heuristique]  rank_C >= 2 (hors extensions)
    La cohomologie de wedge^2 V n'est calculee exactement que pour
    rank_C = 1 (module monad_wedge). Au-dela, le pipeline utilise
    l'estimation h1 * (rk-1) / 2. Le nombre de generations reste valide
    (il ne depend que de h1 et h2 de V), mais le nombre de Higgs ne l'est
    pas. Signale en AVERTISSEMENT, pas en rejet.

[doublon]  meme (cicy, b_charges, c_charges) qu'une entree precedente
    La deduplication du pipeline inclut le champ "type" dans sa cle : une
    monade produite a la fois par le generateur positif et par le
    generateur classique est comptee deux fois. Sans effet sur la validite
    de l'entree, mais gonfle les totaux.

--------------------------------------------------------------------------
Ce que ce script ne peut PAS verifier
--------------------------------------------------------------------------
Le critere de Hoppe tel qu'implemente teste h0(wedge^p V(-H)) = 0 pour
H = e_i uniquement. Le cas H = 0 -- c'est-a-dire h0(wedge^p V) = 0, qui
EST l'enonce du critere de Hoppe pour c1(V) = 0 -- n'est pas teste. Comme
H0(wedge^p V(-H)) est un sous-espace de H0(wedge^p V) pour H ample, le
test implemente est strictement plus FAIBLE que le critere : il laisse
passer des fibres non stables. Les drapeaux h0_nonzero / h3_nonzero
ci-dessus en detectent une partie (le cas p = 1), mais pas les echecs aux
rangs p >= 2, qui demandent un recalcul.
"""
import os, sys, json, argparse
from collections import Counter, defaultdict


def _indices_attendus(r, n_gen=3):
    """
    Valeurs acceptables de |h1 - h2| pour CETTE entree.

    Sur un scan ordinaire : {n_gen}. Sur un scan en mode Wilson, l'indice
    se divise au passage au quotient -- n_gen(X/Gamma) = n_gen(X)/|Gamma| --
    donc l'indice AMONT vaut n_gen * |Gamma|, valeur qui depend de la CICY
    et des groupes qu'elle porte. Les champs `ordres_gamma` et
    `n_gen_amont`, ecrits par le pipeline, portent cette information.

    Sans ces champs (anciens fichiers), on retombe sur {n_gen}.
    """
    cibles = r.get('cibles_chi')
    if cibles:
        return {int(x) for x in cibles}
    ordres = r.get('ordres_gamma')
    if ordres:
        return {n_gen * int(o) for o in ordres}
    amont = r.get('n_gen_amont')
    if amont:
        return {int(amont)}
    return {n_gen}


def audit_record(r, seen_keys, n_gen=3):
    flags = []
    warns = []

    b = r.get('b_charges') or []
    c = r.get('c_charges') or []
    f1 = r.get('f1_charges') or []
    f2 = r.get('f2_charges') or []
    # `cohomology` vaut None sur les extensions dont h1 n'est pas DETERMINE
    # par les bornes. On ne le remplace pas par [0,0,0,0] : ce serait
    # signaler h0/h3 non nuls et un index faux sur une valeur inexistante.
    coh = r.get('cohomology')

    if len(b) and len(c) and (len(b) - len(c)) != r.get('rank_V'):
        flags.append('rank_mismatch')
    # Meme controle pour les extensions : rk(V) = rk(F1) + rk(F2). C'est
    # lui qui avait revele le defaut 4.7 (1571 entrees sur 1571).
    if len(f1) and len(f2) and (len(f1) + len(f2)) != r.get('rank_V'):
        flags.append('rank_mismatch')
    if coh and len(coh) >= 4:
        if coh[0] != 0:
            flags.append('h0_nonzero')
        if coh[3] != 0:
            flags.append('h3_nonzero')
        if abs(coh[1] - coh[2]) not in _indices_attendus(r, n_gen):
            flags.append('index_mismatch')
    if r.get('n_gen') not in _indices_attendus(r, n_gen):
        flags.append('n_gen_not_3')

    if r.get('type') != 'extension' and len(c) >= 2:
        warns.append('wedge2_heuristique')

    # Cle d'identite : (B, C) pour une monade, (F1, F2) pour une extension.
    # Sans le repli sur F1/F2, toutes les extensions d'une meme CICY
    # partageraient la cle vide et seraient signalees doublons.
    key = (r.get('cicy'), r.get('type'),
           tuple(tuple(x) for x in (b or f1)),
           tuple(tuple(x) for x in (c or f2)))
    if key in seen_keys:
        warns.append('doublon')
    seen_keys.add(key)

    return flags, warns


def _sortie_tolerante():
    """
    Empeche un plantage d'encodage sur une console Windows.

    Les etiquettes de groupe de jauge contiennent des indices Unicode
    ("E\u2086"), que la console cp1252 ne sait pas encoder. `print` levait
    alors UnicodeEncodeError, le script mourait sans rien ecrire, et la
    commande suivante de l'enchainement tournait sur un fichier absent.
    Cas reel : scan_exh2, ou `equivariance.py` est mort ainsi et a fait
    echouer `equivariance_f.py` derriere lui.
    """
    import sys as _sys
    for flux in (_sys.stdout, _sys.stderr):
        try:
            flux.reconfigure(errors='replace')
        except (AttributeError, ValueError):
            pass


def main():
    _sortie_tolerante()
    ap = argparse.ArgumentParser(description="Triage d'un results.jsonl deja produit")
    ap.add_argument('output_dir')
    ap.add_argument('--n-gen', type=int, default=3,
                    help="Nombre de generations VOULU (sur le quotient en "
                         "mode Wilson). Defaut: 3.")
    ap.add_argument('--top', type=int, default=10,
                    help="Nombre de candidats propres a afficher (defaut: 10)")
    args = ap.parse_args()

    src = os.path.join(args.output_dir, 'results.jsonl')
    if not os.path.exists(src):
        print(f"Introuvable : {src}")
        return 1

    records = []
    n_lines = n_bad = 0
    with open(src, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            n_lines += 1
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError:
                n_bad += 1

    print(f"\n{'='*66}")
    print(f"  AUDIT DE {src}")
    print(f"{'='*66}")
    print(f"  Lignes lues            : {n_lines}"
          + (f"  ({n_bad} illisibles)" if n_bad else ""))

    by_type = Counter(r.get('type') for r in records)
    print(f"  Par type               : "
          + ", ".join(f"{k}={v}" for k, v in by_type.most_common()))
    ordres = Counter(tuple(r.get('ordres_gamma') or []) for r in records)
    if any(k for k in ordres):
        print(f"  Ordres |Gamma|         : "
              + ", ".join(f"{list(k) or 'aucun'}={v}" for k, v in ordres.most_common(6)))
    cib = Counter(tuple(r.get('cibles_chi') or []) for r in records)
    if any(k for k in cib):
        print(f"  Cibles |chi(V)|        : "
              + ", ".join(f"{list(k)}={v}" for k, v in cib.most_common(6)))

    gens = Counter(r.get('gen', 1) for r in records)
    print(f"  Version de generateur  : "
          + ", ".join(f"v{k}={v}" for k, v in sorted(gens.items())))

    seen = set()
    clean, flagged = [], []
    flag_counts = Counter()
    warn_counts = Counter()
    flag_by_type = defaultdict(Counter)

    for r in records:
        flags, warns = audit_record(r, seen, args.n_gen)
        for fl in flags:
            flag_counts[fl] += 1
            flag_by_type[r.get('type')][fl] += 1
        for w in warns:
            warn_counts[w] += 1
        if flags or 'doublon' in warns:
            r2 = dict(r); r2['flags'] = flags + warns
            flagged.append(r2)
        else:
            r2 = dict(r)
            if warns:
                r2['warnings'] = warns
            clean.append(r2)

    print(f"\n  --- Rejets ---")
    if flag_counts:
        for k, v in flag_counts.most_common():
            print(f"    {k:<20} {v:>7}  ({100*v/max(1,len(records)):.1f} %)")
    else:
        print("    aucun")

    print(f"\n  --- Avertissements (entree conservee) ---")
    if warn_counts:
        for k, v in warn_counts.most_common():
            print(f"    {k:<20} {v:>7}")
    else:
        print("    aucun")

    if flag_by_type:
        print(f"\n  --- Repartition des rejets par type de fibre ---")
        for t in by_type:
            fc = flag_by_type.get(t)
            if fc:
                print(f"    {t:<12} " + ", ".join(f"{k}={v}" for k, v in fc.most_common()))

    print(f"\n{'='*66}")
    print(f"  Candidats retenus      : {len(clean)} / {len(records)}"
          f"  ({100*len(clean)/max(1,len(records)):.1f} %)")
    print(f"  Candidats ecartes      : {len(flagged)}")

    if clean:
        cicys = sorted({r['cicy'] for r in clean})
        print(f"  CICYs concernees       : {len(cicys)}")
        print(f"  Groupes de jauge       : "
              + ", ".join(f"{k}={v}" for k, v in
                          Counter(r['gauge'] for r in clean).most_common()))

        exact = [r for r in clean if 'warnings' not in r]
        print(f"  Dont Higgs fiable      : {len(exact)}"
              f"  (rank_C = 1, wedge2 exact)")

        pool = exact or clean
        # `exotics` peut valoir None (non calcule). Le tri le renvoie en
        # queue plutot que de le confondre avec 0, qui vaudrait « propre ».
        pool.sort(key=lambda r: (r.get('higgs', 99) if r.get('higgs', 0) > 0 else 99,
                                 -r.get('score', 0)))
        print(f"\n  Top {args.top} (Higgs faible d'abord, Higgs > 0 requis "
              f"pour briser SO(10)/SU(5)) :")
        print(f"    {'#CICY':>6} {'type':<10} {'jauge':>7} {'rk':>2} "
              f"{'H':>4} {'exo':>3} {'cohomologie':>16} {'score':>6}")
        for r in pool[:args.top]:
            print(f"    {r['cicy']:>6} {r['type']:<10} {r['gauge']:>7} "
                  f"{r['rank_V']:>2} {r.get('higgs',0):>4} "
                  f"{('?' if r.get('exotics') is None else r.get('exotics',0)):>3} "
                  f"{str(r.get('cohomology')):>16} {r.get('score',0):>6}")

    clean_path = os.path.join(args.output_dir, 'results_clean.jsonl')
    flag_path = os.path.join(args.output_dir, 'results_flagged.jsonl')
    with open(clean_path, 'w', encoding='utf-8') as f:
        for r in clean:
            f.write(json.dumps(r) + '\n')
    with open(flag_path, 'w', encoding='utf-8') as f:
        for r in flagged:
            f.write(json.dumps(r) + '\n')
    print(f"\n  Ecrit : {clean_path}")
    print(f"  Ecrit : {flag_path}")
    print(f"{'='*66}\n")
    return 0


if __name__ == '__main__':
    sys.exit(main())
