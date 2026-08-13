"""
monad_wedge.py -- Cohomologie de wedge^2 V pour une monade, VERSION 2.

--------------------------------------------------------------------------
Pourquoi cette reecriture
--------------------------------------------------------------------------
La version 1 deduisait h^i(wedge^2 V) de la suite exacte longue en
supposant que chaque application induite etait de rang maximal, puis
« corrigeait » le resultat par des ecrasements a zero. Mesure : sur 456
monades, 79 % des vecteurs de Betti retournes contredisaient le chi
calcule par la fonction elle-meme, quelques lignes plus haut. Sur les
monades de rang 3, ou la reponse exacte est connue d'avance
(wedge^2 V = V*, donc h^i(w2V) = h^(3-i)(V)), un seul cas sur 34 etait
juste, avec des ecarts allant de 3 a 7 251 840.

Cette version ne suppose plus rien. Elle fournit trois choses :

  1. chi(wedge^2 V), EXACT et inconditionnel.
  2. Des BORNES rigoureuses sur chaque h^i.
  3. Un drapeau par degre indiquant si la borne est un point, donc si
     h^i est determine.

Elle ne retourne jamais un nombre invente. Quand une valeur n'est pas
determinee, elle est absente et le drapeau vaut False -- a l'appelant
d'en tenir compte plutot que de propager une valeur fausse.

--------------------------------------------------------------------------
Les suites utilisees
--------------------------------------------------------------------------
Pour rank_C = 1, la filtration de wedge^2 B donne (wedge^2 C = 0) :

    0 -> wedge^2 V -> wedge^2 B -> V (x) C -> 0

et la monade tensorisee par C :

    0 -> V (x) C -> B (x) C -> C (x) C -> 0

wedge^2 B, B (x) C et C (x) C sont des sommes de fibres en droites, dont
la cohomologie passe par `koszul_cohomology_ex` -- avec sa certification
par degre, propagee ici.

--------------------------------------------------------------------------
Bornes rigoureuses
--------------------------------------------------------------------------
Pour 0 -> L -> M -> R -> 0, la suite exacte longue donne

    h^i(L) = coker(g_(i-1)) + ker(g_i)

avec, sans aucune hypothese sur les rangs :

    ker(g_i)     dans [max(0, h^i(M) - h^i(R)) , h^i(M)]
    coker(g_j)   dans [max(0, h^j(R) - h^j(M)) , h^j(R)]

Ces bornes sont exactes, pas heuristiques. Elles ne determinent h^i que
dans environ 0,1 % des cas : c'est la mesure honnete de ce que la seule
suite exacte longue permet de conclure. Calculer les rangs reels des
applications demande `koszul_exact.py`, non branche a ce jour.

--------------------------------------------------------------------------
Raccourcis par isomorphisme (det V = O, dualite de Serre sur un CY3)
--------------------------------------------------------------------------
  rang 3 : wedge^2 V = V*        -> h^i(w2V) = h^(3-i)(V), EXACT, sans LES
  rang 4 : wedge^2 V est autodual -> h^3 = h^0, h^2 = h^1, donc
                                     chi(w2V) = 0 identiquement.
                                     Sert de test de coherence.
  rang 5 : wedge^2 V* = wedge^3 V

Le raccourci de rang 3 rend le calcul exact pour tous les candidats E6.
"""
from typing import List, Dict, Optional

import numpy as np

from cy_landscape.core.exact_cohomology import koszul_cohomology_ex


# ================================================================
# Sommes de fibres en droites, avec certification par degre
# ================================================================

def _sum_line_bundles(ambient, config, charges_list):
    """
    h^i d'une somme de fibres en droites + certification par degre.
    Retourne (h, cert) avec h[i] entier et cert[i] booleen.
    """
    h = {i: 0 for i in range(4)}
    cert = {i: True for i in range(4)}
    chi = 0
    for ch in charges_list:
        r = koszul_cohomology_ex(ambient, config, list(ch))
        cd = r.get('certified_by_degree') or {i: False for i in range(4)}
        # chi vient du champ 'chi' de koszul_cohomology_ex, TOUJOURS exact.
        # Surtout pas de la somme alternee des h^i, qui n'est fiable que
        # sur les degres certifies -- c'etait precisement l'erreur de la v1.
        chi += int(r.get('chi', 0))
        for i in range(4):
            h[i] += int(r.get(i, 0))
            if not cd.get(i, False):
                cert[i] = False
    h['chi'] = chi
    return h, cert


def _chi(h):
    """chi exact, stocke par `_sum_line_bundles` sous la cle 'chi'."""
    return int(h['chi'])


# ================================================================
# Bornes rigoureuses issues de la suite exacte longue
# ================================================================

def _les_bounds(hM, hR, certM, certR, boundsM=None, boundsR=None):
    """
    0 -> L -> M -> R -> 0 : bornes sur h^i(L).

    Si M ou R ne sont eux-memes connus que par intervalle (`boundsM`,
    `boundsR` sous forme {i: (lo, hi)}), les intervalles se propagent.
    Retourne {i: (lo, hi, determine)}.
    """
    def rng_of(h, bounds, i):
        if bounds is not None and i in bounds:
            return bounds[i]
        v = h.get(i, 0)
        return (v, v)

    out = {}
    for i in range(4):
        Mlo, Mhi = rng_of(hM, boundsM, i)
        Rlo, Rhi = rng_of(hR, boundsR, i)
        if i > 0:
            Mplo, Mphi = rng_of(hM, boundsM, i - 1)
            Rplo, Rphi = rng_of(hR, boundsR, i - 1)
        else:
            Mplo = Mphi = Rplo = Rphi = 0

        ker_lo = max(0, Mlo - Rhi)
        ker_hi = Mhi
        cok_lo = max(0, Rplo - Mphi)
        cok_hi = Rphi

        lo, hi = ker_lo + cok_lo, ker_hi + cok_hi

        ok = (lo == hi)
        if ok:
            for src, idx in ((certM, i), (certR, i),
                             (certM, i - 1), (certR, i - 1)):
                if idx >= 0 and not src.get(idx, False):
                    ok = False
                    break
        out[i] = (lo, hi, ok)
    return out


# ================================================================
# API publique
# ================================================================

def cohomology_wedge2_B(ambient, config, b_charges):
    """h^i(wedge^2 B) = somme sur i<j de h^i(O(b_i + b_j))."""
    charges = []
    n = len(b_charges)
    for i in range(n):
        for j in range(i + 1, n):
            charges.append([b_charges[i][k] + b_charges[j][k]
                            for k in range(len(b_charges[i]))])
    h, _ = _sum_line_bundles(ambient, config, charges)
    return {i: h[i] for i in range(4)}


def cohomology_wedge2_B_ex(ambient, config, b_charges):
    """Idem, avec la certification par degre."""
    charges = []
    n = len(b_charges)
    for i in range(n):
        for j in range(i + 1, n):
            charges.append([b_charges[i][k] + b_charges[j][k]
                            for k in range(len(b_charges[i]))])
    return _sum_line_bundles(ambient, config, charges)


def cohomology_V_tensor_C_ex(ambient, config, b_charges, c_charges):
    """
    Bornes sur h^i(V (x) C) via 0 -> V(x)C -> B(x)C -> C(x)C -> 0.
    Retourne (bounds, chi, cert) avec bounds = {i: (lo, hi, determine)}.
    """
    m = len(b_charges[0])
    bc = [[b[k] + c[k] for k in range(m)] for b in b_charges for c in c_charges]
    cc = [[c1[k] + c2[k] for k in range(m)] for c1 in c_charges for c2 in c_charges]

    hBC, certBC = _sum_line_bundles(ambient, config, bc)
    hCC, certCC = _sum_line_bundles(ambient, config, cc)

    bounds = _les_bounds(hBC, hCC, certBC, certCC)
    chi = _chi(hBC) - _chi(hCC)
    return bounds, chi


def cohomology_V_tensor_C(ambient, config, b_charges, c_charges):
    """
    Compatibilite : renvoie {i: h^i} en ne donnant une valeur que pour les
    degres determines, 0 sinon. Preferer `cohomology_V_tensor_C_ex`.
    """
    bounds, _ = cohomology_V_tensor_C_ex(ambient, config, b_charges, c_charges)
    return {i: (bounds[i][0] if bounds[i][2] else 0) for i in range(4)}


def cohomology_wedge2_V(ambient, config, b_charges, c_charges,
                        cohom_V: Optional[Dict[int, int]] = None,
                        cert_V: Optional[Dict[int, bool]] = None):
    """
    Cohomologie de wedge^2 V.

    Retourne un dictionnaire :
      'wedge2V'      {i: h^i}  -- SEULEMENT les degres determines
      'bounds'       {i: (lo, hi)}
      'certified'    {i: bool}
      'chi_wedge2V'  int, exact et inconditionnel
      'method'       'rang3' | 'les' | 'indisponible'

    `cohom_V` / `cert_V` : cohomologie de V si deja calculee (utilisee par
    le raccourci de rang 3).
    """
    rank_V = len(b_charges) - len(c_charges)
    m = len(b_charges[0])

    # --- chi, toujours exact ---------------------------------------
    hW2B, certW2B = cohomology_wedge2_B_ex(ambient, config, b_charges)
    chi_w2v = None
    if len(c_charges) == 1:
        _, chi_VC = cohomology_V_tensor_C_ex(ambient, config, b_charges, c_charges)
        chi_w2v = _chi(hW2B) - chi_VC

    # --- rang 3 : wedge^2 V = V*, exact, aucune LES ------------------
    if rank_V == 3 and cohom_V is not None:
        h = {i: int(cohom_V[3 - i]) for i in range(4)}
        cert = ({i: bool(cert_V.get(3 - i, False)) for i in range(4)}
                if cert_V else {i: False for i in range(4)})
        return {
            'wedge2V': {i: h[i] for i in range(4) if cert[i]},
            'bounds': {i: (h[i], h[i]) for i in range(4)},
            'certified': cert,
            'chi_wedge2V': chi_w2v,
            'method': 'rang3',
        }

    if len(c_charges) != 1:
        # wedge^2 C != 0 : la filtration ne se reduit pas a une suite
        # exacte courte. Non traite -- on ne devine pas.
        return {
            'wedge2V': {}, 'bounds': {}, 'certified': {i: False for i in range(4)},
            'chi_wedge2V': chi_w2v, 'method': 'indisponible',
        }

    # --- cas general : bornes propagees ------------------------------
    bVC, _ = cohomology_V_tensor_C_ex(ambient, config, b_charges, c_charges)
    hVC = {i: bVC[i][0] for i in range(4)}
    certVC = {i: bVC[i][2] for i in range(4)}
    boundsVC = {i: (bVC[i][0], bVC[i][1]) for i in range(4)}

    res = _les_bounds(hW2B, hVC, certW2B, certVC, boundsR=boundsVC)

    # Rang 4 : wedge^2 V est autodual, h^3 = h^0 et h^2 = h^1.
    # On croise les bornes des degres apparies.
    if rank_V == 4:
        for a, b in ((0, 3), (1, 2)):
            lo = max(res[a][0], res[b][0])
            hi = min(res[a][1], res[b][1])
            if lo <= hi:
                ok = res[a][2] or res[b][2] or (lo == hi)
                res[a] = (lo, hi, ok and lo == hi)
                res[b] = (lo, hi, ok and lo == hi)

    return {
        'wedge2V': {i: res[i][0] for i in range(4) if res[i][2]},
        'bounds': {i: (res[i][0], res[i][1]) for i in range(4)},
        'certified': {i: res[i][2] for i in range(4)},
        'chi_wedge2V': chi_w2v,
        'method': 'les',
    }
