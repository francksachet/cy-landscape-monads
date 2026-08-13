"""
koszul_exact.py -- Cohomologie Koszul avec verification de degenerescence E1
et calcul des differentielles quand necessaire.

Le spectral sequence de Koszul a :
  E_1^{-p, q} = H^q(A, ∧^p N* ⊗ L)  pour p = 0, ..., K
  d_1 : E_1^{-p, q} -> E_1^{-p+1, q}  (induit par la differentielle de Koszul)

Converge vers H^{q-p}(X, L|_X).

Degenerescence E1 : si pour chaque n = q-p, au plus un E_1^{-p,q} est non nul,
alors H^n(X,L) = cet unique terme (pas de differentielle possible).

Sinon : il faut calculer E_2 = ker(d1)/im(d1), etc.

Pour le calcul des d1 : d1 est la map sur la cohomologie induite par
la contraction avec les equations definissantes. En general, pour une 
intersection complete GENERIQUE, d1 est de rang maximal (Bott-type argument).

Ref: Green, "Koszul cohomology and the geometry of projective varieties"
     Blumenhagen, Jurke, Rahn, Roschy "cohomCalg" (2010, arXiv:1003.5217)
"""
import numpy as np
from itertools import combinations
from typing import Dict, Tuple

from cy_landscape.core.exact_cohomology import h_ambient


def koszul_E1_table(ambient_dims, config_matrix, charges):
    """
    Calcule la table E_1^{-p, q} complete.

    Retourne:
      e1[p][q] = dim H^q(A, ∧^p N* ⊗ L)
      K = nombre de contraintes
      dim_amb = dimension de l'espace ambiant
    """
    config = np.array(config_matrix)
    K = config.shape[0]
    m = len(ambient_dims)
    dim_amb = sum(ambient_dims)

    # e1[p][q] pour p = 0..K, q = 0..dim_amb
    e1 = np.zeros((K + 1, dim_amb + 1), dtype=np.int64)

    for p in range(K + 1):
        if p == 0:
            subsets = [()]
        else:
            subsets = list(combinations(range(K), p))

        for S in subsets:
            # charges decalees : L ⊗ ∧^p N* -> shift by -sum of rows in S
            shifted = list(charges)
            for a in S:
                shifted = [shifted[r] - int(config[a, r]) for r in range(m)]

            h = h_ambient(ambient_dims, shifted)
            for q in range(dim_amb + 1):
                e1[p, q] += int(h.get(q, 0))

    return e1, K, dim_amb


def check_E1_degeneration(e1, K, dim_amb):
    """
    Verifie si E1 degenere.

    d1: E1^{-p, q} -> E1^{-p+1, q} opere a q FIXE.
    Donc la degenerescence echoue si pour un q donne,
    il y a deux termes E1 non nuls a des p differents
    (source et cible d'un d1 potentiel).

    Retourne (degenere, details).
    """
    issues = []

    for q in range(dim_amb + 1):
        nonzero = []
        for p in range(K + 1):
            if e1[p, q] > 0:
                nonzero.append((p, q, int(e1[p, q])))

        if len(nonzero) > 1:
            issues.append({'q': q, 'terms': nonzero})

    return len(issues) == 0, issues


def koszul_cohomology_exact(ambient_dims, config_matrix, charges):
    """
    Calcule H^n(X, L) avec gestion de la non-degenerescence E1.

    Quand E1 ne degenere pas, utilise l'hypothese de rang maximal
    pour d1 (valide pour des equations generiques).
    
    Retourne: {n: h^n} pour n = 0, 1, 2, 3
              plus 'degenerate': bool et 'warnings': list
    """
    e1, K, dim_amb = koszul_E1_table(ambient_dims, config_matrix, charges)
    degenerate, issues = check_E1_degeneration(e1, K, dim_amb)

    result = {0: 0, 1: 0, 2: 0, 3: 0}
    warnings = []

    if degenerate:
        # Simple: un seul terme par degre total
        for n in range(4):
            for p in range(K + 1):
                q = n + p
                if 0 <= q <= dim_amb:
                    result[n] += int(e1[p, q])
    else:
        # Non-degenere: calculer E2 via rang maximal de d1
        # d1: E1^{-p, q} -> E1^{-p+1, q}  (meme q, p diminue de 1)
        # Pour chaque q, la chaine est :
        #   E1^{-K,q} -d1-> E1^{-K+1,q} -d1-> ... -d1-> E1^{0,q}
        # Avec d1 de rang maximal (equations generiques).
        
        e2 = np.zeros_like(e1)
        
        for q in range(dim_amb + 1):
            # chain[p] = dim E1^{-p, q} pour p = 0, ..., K
            # d1 va de p vers p-1 (ie de E1^{-p,q} vers E1^{-(p-1),q})
            chain = [int(e1[p, q]) for p in range(K + 1)]
            
            # Rang maximal des differentielles
            # d1^{-p}: chain[p] -> chain[p-1], rang = min(source, target apres soustraction)
            remaining = list(chain)
            for p in range(K, 0, -1):  # de p=K a p=1
                rank_d = min(remaining[p], remaining[p-1])
                remaining[p] -= rank_d     # ker(d1) a la position p
                remaining[p-1] -= rank_d   # espace restant apres im(d1)
            
            for p in range(K + 1):
                e2[p, q] = max(0, remaining[p])
        
        # H^n(X, L) = somme des E2^{-p, n+p} pour p = 0..K
        for n in range(4):
            for p in range(K + 1):
                q = n + p
                if 0 <= q <= dim_amb:
                    result[n] += int(e2[p, q])
        
        for issue in issues:
            q_val = issue['q']
            warnings.append(
                f"q={q_val}: {len(issue['terms'])} termes E1 non nuls, "
                f"d1 rang-maximal applique")

    result['degenerate'] = degenerate
    result['warnings'] = warnings
    return result


def validate_quintic():
    """Valide sur le quintique (P4, config=[5]) - resultats connus."""
    ambient = [4]
    config = np.array([[5]])

    print("VALIDATION : Quintique dans P^4")
    print(f"{'O(n)':>6} {'H0':>4} {'H1':>4} {'H2':>4} {'H3':>4} {'chi':>5} {'E1?':>4} {'warn':>5}")

    # Valeurs connues (Hubsch, "Calabi-Yau Manifolds" Table A.1)
    known = {
        -5: (0, 0, 0, 1), -4: (0, 0, 0, 0), -3: (0, 0, 0, 0),
        -2: (0, 0, 0, 0), -1: (0, 0, 0, 5),
        0: (1, 0, 0, 1), 1: (5, 0, 0, 0), 2: (15, 0, 0, 0),
        3: (35, 0, 0, 0), 4: (70, 0, 0, 0), 5: (125, 0, 0, 0),
    }

    all_ok = True
    for n in range(-5, 6):
        h = koszul_cohomology_exact(ambient, config, [n])
        chi = sum((-1)**i * h[i] for i in range(4))
        deg = 'Y' if h['degenerate'] else 'N'
        warn = len(h['warnings'])

        if n in known:
            expected = known[n]
            match = all(h[i] == expected[i] for i in range(4))
            status = "OK" if match else f"FAIL (attendu {expected})"
            if not match: all_ok = False
        else:
            status = ""

        print(f"  O({n:>2}) {h[0]:>4} {h[1]:>4} {h[2]:>4} {h[3]:>4} {chi:>5} {deg:>4} {warn:>5}  {status}")

    return all_ok


if __name__ == "__main__":
    ok = validate_quintic()
    print(f"\nValidation quintique: {'PASSED' if ok else 'FAILED'}")
