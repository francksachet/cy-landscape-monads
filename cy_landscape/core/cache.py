"""
cache.py -- Cache de cohomologie pour eviter les recalculs.

Le goulot principal est le Koszul : chaque appel coute ~1ms,
et le scan appelle koszul_cohomology des milliers de fois
avec les memes arguments (meme CICY, memes charges).

Ce cache est local a une CICY (reinitialise entre les varietes).
"""
from functools import lru_cache
import numpy as np

_current_ambient = None
_current_config = None
_cache = {}


def set_geometry(ambient, config):
    """Definit la geometrie courante et vide le cache."""
    global _current_ambient, _current_config, _cache
    _current_ambient = tuple(ambient)
    _current_config = config.tobytes() if isinstance(config, np.ndarray) else bytes(str(config), 'utf-8')
    _cache = {}


def cached_koszul(ambient, config, charges):
    """Koszul avec cache. ~100x plus rapide sur les appels repetes."""
    key = (tuple(charges),)

    if key in _cache:
        return _cache[key]

    from cy_landscape.core.exact_cohomology import koszul_cohomology
    result = koszul_cohomology(ambient, config, charges)
    _cache[key] = result
    return result


def cache_stats():
    return len(_cache)
