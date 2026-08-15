#!/usr/bin/env python3
"""
resume_cible.py -- resume les sorties de `equivariance_f.py --cicy N`.

Un run cible produit des milliers de lignes identiques : autant de monades
distinctes qui rendent le meme verdict. Ce script les replie et repond a la
seule question posee -- le candidat survit-il, avec quel Gamma, quel lambda,
et combien de generations sur le quotient.

Usage:
    python resume_cible.py scan_w4_c6890 scan_w4_c6947 scan_w4_c6715
"""
import os, sys, json, collections


def _sortie_tolerante():
    for flux in (sys.stdout, sys.stderr):
        try:
            flux.reconfigure(errors='replace')
        except (AttributeError, ValueError):
            pass


def resume(dossier):
    src = os.path.join(dossier, 'results_equivariance_f.jsonl')
    if not os.path.exists(src):
        print(f"\n  {dossier} : {src} absent "
              f"(le run n'est pas alle jusqu'au bout -- rien n'est ecrit "
              f"avant la derniere ligne)")
        return
    rs = [json.loads(l) for l in open(src, encoding='utf-8') if l.strip()]
    print(f"\n{'=' * 78}")
    print(f"  {dossier}   ({len(rs)} lignes)")
    print(f"{'=' * 78}")

    # Repli : ce qui distingue deux lignes, c'est le VERDICT, pas les charges.
    # Le nombre de monades derriere chaque verdict est reporte -- il dit si le
    # resultat tient a une configuration ou a des centaines.
    cle = lambda r: (r.get('cicy'), r.get('gauge'), r.get('rank_V'),
                     r.get('groupe'), str(r.get('lambda')),
                     bool(r.get('survit')), bool(r.get('indetermine')),
                     r.get('n_gen_quotient'), r.get('etat'),
                     str(r.get('surjectif_degre')))
    c = collections.Counter(cle(r) for r in rs)

    surv = [(k, v) for k, v in c.items() if k[5]]
    if surv:
        print(f"\n  SURVIVANTS")
        for k, v in sorted(surv, key=lambda kv: -kv[1]):
            print(f"    #{k[0]}  {k[1]}  rang {k[2]}  Gamma={k[3]}  "
                  f"lambda={k[4]}  ->  {k[7]} gen sur X/Gamma")
            print(f"       {v} monades, degre temoin de surjectivite {k[9]}")
    else:
        print(f"\n  AUCUN SURVIVANT")

    print(f"\n  Repartition des verdicts :")
    for k, v in c.most_common(12):
        # `etat == 'ok'` ne veut pas dire « retenu » : il veut dire
        # « evalue jusqu'au bout ». Le distinguer evite de lire
        # « elimine : ok », qui ne dit rien.
        motif = '' if k[8] == 'ok' else f" : {k[8]}"
        if k[5]:
            etat = f"SURVIT ({k[7]} gen sur X/Gamma)"
        elif k[6]:
            etat = f"indetermine{motif or ' (un test non calculable)'}"
        else:
            etat = f"elimine{motif or ' (stabilite equivariante en defaut)'}"
        print(f"    {v:>6} x  Gamma={str(k[3]):<10} lambda={k[4]:<8} {etat}")

    n_surv = sum(v for k, v in c.items() if k[5])
    n_ind = sum(v for k, v in c.items() if k[6] and not k[5])
    print(f"\n  Lignes : {n_surv} survivantes, {n_ind} indeterminees, "
          f"{len(rs) - n_surv - n_ind} eliminees")
    print(f"  Verdicts DISTINCTS : {len(c)}  <- si ce nombre est tres petit")
    print(f"  devant le nombre de lignes, les monades enumerees sont")
    print(f"  redondantes entre elles (meme fibre a symetrie pres).")


def main():
    _sortie_tolerante()
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    for d in sys.argv[1:]:
        resume(d)
    print()
    return 0


if __name__ == '__main__':
    sys.exit(main())
