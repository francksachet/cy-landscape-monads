#!/usr/bin/env python3
"""
migrate_checkpoint.py -- Convertit un ancien checkpoint.pkl (format unique)
vers le nouveau format separe progress.pkl + results.jsonl.

Aucune perte de donnees : tous les resultats et toutes les CICYs deja
traitees sont preserves.

Usage:
  python migrate_checkpoint.py                          # dossier par defaut
  python migrate_checkpoint.py output_optimized          # dossier explicite
"""
import sys, os, json, pickle, time


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


def migrate(output_dir):
    old_path = os.path.join(output_dir, 'checkpoint.pkl')

    if not os.path.exists(old_path):
        print(f"Aucun ancien checkpoint.pkl trouve dans {output_dir}")
        print("Rien a migrer -- le nouveau format sera cree automatiquement")
        print("au prochain lancement de main_optimized.py.")
        return

    print(f"Lecture de l'ancien checkpoint : {old_path}")
    with open(old_path, 'rb') as f:
        old_data = pickle.load(f)

    done_cicys = old_data['done_cicys']
    all_results = old_data['all_results']
    elapsed = old_data.get('elapsed', 0.0)
    params = old_data.get('params', None)

    print(f"  {len(done_cicys)} CICYs traitees")
    print(f"  {len(all_results)} resultats")
    print(f"  {elapsed/60:.1f} min deja consommees")
    if params:
        print(f"  Parametres : {params}")
    else:
        print(f"  Parametres : non enregistres (ancien format)")
        print(f"  ATTENTION : il faudra fournir les MEMES parametres au")
        print(f"  prochain lancement (--max-ps, --n-random, etc.) pour")
        print(f"  garantir la coherence des resultats.")

    # --- Ecriture du nouveau progress.pkl (leger) ---
    progress_path = os.path.join(output_dir, 'progress.pkl')
    progress_data = {
        'done_cicys': done_cicys,
        'elapsed': elapsed,
        'saved_at': time.time(),
        'params': params,
    }
    tmp = progress_path + '.tmp'
    with open(tmp, 'wb') as f:
        pickle.dump(progress_data, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, progress_path)
    print(f"\n  Ecrit : {progress_path}")

    # --- Ecriture du nouveau results.jsonl (append-only) ---
    results_path = os.path.join(output_dir, 'results.jsonl')
    if os.path.exists(results_path):
        print(f"  ATTENTION : {results_path} existe deja.")
        resp = input(f"  L'ecraser ? (o/N) : ").strip().lower()
        if resp != 'o':
            print("  Migration annulee pour results.jsonl (progress.pkl deja migre).")
            return

    with open(results_path, 'w', encoding='utf-8') as f:
        for r in all_results:
            f.write(json.dumps(r, cls=NpEncoder) + '\n')
    print(f"  Ecrit : {results_path} ({len(all_results)} lignes)")

    # --- Renommer l'ancien fichier (ne pas le supprimer, par securite) ---
    backup_path = old_path + '.migrated'
    os.rename(old_path, backup_path)
    print(f"\n  Ancien checkpoint renomme en securite : {backup_path}")
    print(f"  (tu peux le supprimer manuellement une fois verifie que tout va bien)")

    print(f"\n{'='*60}")
    print(f"  MIGRATION TERMINEE")
    print(f"{'='*60}")
    print(f"  Tu peux relancer main_optimized.py normalement :")
    print(f"  il reprendra automatiquement depuis {len(done_cicys)} CICYs traitees,")
    print(f"  avec un demarrage rapide grace au nouveau format.")


if __name__ == "__main__":
    output_dir = sys.argv[1] if len(sys.argv) > 1 else 'output_optimized'
    migrate(output_dir)
