"""
parse_oxford.py -- Parseur pour cicylist.txt (Candelas et al., Oxford)
Telecharger depuis : https://www-thphys.physics.ox.ac.uk/projects/CalabiYau/cicylist/cicylist.txt

Usage:
    from cy_landscape.data.parse_oxford import load_oxford_file
    cicys = load_oxford_file("chemin/vers/cicylist.txt")
"""
import re
import numpy as np


def parse_cicy_block(lines):
    """Parse un bloc d'entree CICY (entre deux 'Num :')."""
    result = {}
    matrix_rows = []

    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.strip().startswith('NumPol') and ':' in line:
            result['num_pol'] = int(line.split(':')[1].strip())
        elif line.strip().startswith('NumPs') and ':' in line:
            result['num_ps'] = int(line.split(':')[1].strip())
        elif line.strip().startswith('Num') and ':' in line:
            result['num'] = int(line.split(':')[1].strip())
        elif line.strip().startswith('Eta') and ':' in line:
            result['chi'] = int(line.split(':')[1].strip())
        elif line.strip().startswith('H11') and ':' in line:
            result['h11'] = int(line.split(':')[1].strip())
        elif line.strip().startswith('H21') and ':' in line:
            result['h21'] = int(line.split(':')[1].strip())
        elif line.strip().startswith('C2') and ':' in line:
            c2_str = line.split(':')[1].strip()
            result['c2'] = np.array([int(x) for x in re.findall(r'-?\d+', c2_str)],
                                    dtype=float)
        elif line.startswith('Redun'):
            pass
        elif line.startswith('{') and 'Redun' not in line:
            row = [int(x) for x in re.findall(r'-?\d+', line)]
            if row:
                matrix_rows.append(row)

    if not matrix_rows or 'num_ps' not in result:
        return None

    m = result['num_ps']
    K = result['num_pol']

    if len(matrix_rows) != m:
        return None
    if any(len(r) != K for r in matrix_rows):
        return None

    oxford_matrix = np.array(matrix_rows, dtype=int)   # m x K (espaces x polynomes)
    config = oxford_matrix.T                             # K x m (polynomes x espaces)
    ambient = [int(oxford_matrix[i].sum() - 1) for i in range(m)]

    # Verification CY : sum(n_i) - K == 3
    dim = sum(ambient) - K
    if dim != 3:
        return None

    result['ambient'] = ambient
    result['config'] = config
    return result


def load_oxford_file(filepath):
    """
    Charge le fichier cicylist.txt complet et retourne une liste de dicts.

    Chaque dict contient :
        num      : numero dans la liste de Candelas
        ambient  : [n_1, ..., n_m]  dimensions des espaces projectifs
        config   : np.array K x m   matrice de configuration (polynomes x espaces)
        h11, h21 : nombres de Hodge
        chi      : caracteristique d'Euler = 2*(h11 - h21)
        c2       : np.array   seconds nombres de Chern c2(TY).J_i
    """
    # Lecture robuste (BOM, encodage, fins de ligne Windows)
    with open(filepath, 'r', encoding='utf-8-sig', errors='replace') as f:
        text = f.read()
    
    # Normaliser les fins de ligne
    text = text.replace('\r\n', '\n').replace('\r', '\n')

    # Decouper en blocs : chercher "Num" suivi d'espaces et ":"
    blocks = re.split(r'(?=Num\s*:)', text)
    entries = []
    for block in blocks:
        block = block.strip()
        if not block or not block.startswith('Num'):
            continue
        entry = parse_cicy_block(block.split('\n'))
        if entry is not None:
            entries.append(entry)

    return entries


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python parse_oxford.py chemin/vers/cicylist.txt")
        print("Telecharger depuis:")
        print("  https://www-thphys.physics.ox.ac.uk/projects/CalabiYau/cicylist/cicylist.txt")
        sys.exit(1)

    path = sys.argv[1]
    entries = load_oxford_file(path)
    print(f"Parsed {len(entries)} CICYs from {path}")

    # Stats
    from collections import Counter
    chi_dist = Counter(e['chi'] for e in entries)
    print(f"\nDistribution chi (Euler):")
    for chi in sorted(chi_dist.keys(), key=abs):
        n = chi_dist[chi]
        dh = abs(chi) // 2
        print(f"  chi={chi:>5}  |h11-h21|={dh:>3}  count={n:>4}")

    print(f"\nchi=-6 (3 generations topologiques): "
          f"{'PRESENT' if -6 in chi_dist else 'ABSENT'}")
