#!/usr/bin/env python3
"""
export_checkpoint.py -- Exporte un results.json depuis un checkpoint en cours.

Permet de visualiser les resultats partiels d'un scan main_optimized.py
sans attendre qu'il soit termine.

Usage:
  python export_checkpoint.py                              # dossier par defaut
  python export_checkpoint.py output_optimized              # dossier explicite
  python export_checkpoint.py output_optimized/checkpoint.pkl
"""
import sys, os, json, pickle


def load_checkpoint(path):
    """Charge un fichier checkpoint.pkl."""
    if os.path.isdir(path):
        path = os.path.join(path, 'checkpoint.pkl')

    if not os.path.exists(path):
        print(f"ERREUR : checkpoint introuvable : {path}")
        sys.exit(1)

    with open(path, 'rb') as f:
        data = pickle.load(f)
    return data, path


def export(data, output_path):
    """Ecrit un results.json au meme format que main_optimized.py."""
    all_results = sorted(data['all_results'], key=lambda r: r['score'], reverse=True)
    n_stable = len(all_results)
    elapsed = data.get('elapsed', 0.0)
    params = data.get('params', {})

    export_data = {
        'parameters': params,
        'n_cicys': len(data['done_cicys']),
        'n_stable': n_stable,
        'time_seconds': round(elapsed, 1),
        'results': all_results[:500],
        'partial': True,  # Marque que ce n'est pas un scan complet
    }

    class NpEncoder(json.JSONEncoder):
        def default(self, obj):
            try:
                import numpy as np
                if isinstance(obj, (np.integer,)): return int(obj)
                if isinstance(obj, (np.floating,)): return float(obj)
                if isinstance(obj, np.ndarray): return obj.tolist()
            except ImportError:
                pass
            return super().default(obj)

    with open(output_path, 'w') as f:
        json.dump(export_data, f, indent=2, cls=NpEncoder)

    return export_data


def main():
    arg = sys.argv[1] if len(sys.argv) > 1 else 'output_optimized'
    data, ckpt_path = load_checkpoint(arg)

    out_dir = os.path.dirname(ckpt_path)
    output_path = os.path.join(out_dir, 'results.json')

    export_data = export(data, output_path)

    print(f"{'='*60}")
    print(f"  EXPORT DEPUIS CHECKPOINT (scan partiel)")
    print(f"{'='*60}")
    print(f"  CICYs traitees jusqu'ici : {export_data['n_cicys']}")
    print(f"  Fibres Hoppe-stables     : {export_data['n_stable']}")
    print(f"  Temps deja consomme      : {export_data['time_seconds']/60:.1f} min")
    print(f"  Parametres               : {export_data['parameters']}")
    print(f"\n  Export : {output_path}")
    print(f"\n  ATTENTION : ce fichier ne contient que les CICYs deja")
    print(f"  traitees. Le scan complet donnera davantage de resultats.")
    print(f"  Relance le scan normalement pour continuer, il reprendra")
    print(f"  automatiquement depuis ce checkpoint.")

    # Bonus : generer aussi les graphiques directement
    try:
        from cy_landscape.core.visualize import generate_all_plots
        print(f"\n  Generation des graphiques...")
        generate_all_plots(output_path, out_dir)
    except Exception as e:
        print(f"  Graphiques ignores : {e}")
        print(f"  (Le notebook peut aussi les generer depuis {output_path})")


if __name__ == "__main__":
    main()
