"""Resolve a PDB mmCIF path: use the repo copy if present, else cache a download.

Structures that are not part of the benchmark set (e.g. modern receptor
references used only to derive a numbering offset) are deliberately not kept in
the repository; they are fetched on demand into a gitignored cache so that
every number in data/README.md stays reproducible.
"""
import os, urllib.request

_HERE = os.path.dirname(os.path.abspath(__file__))
COORD = os.path.join(_HERE, '..', 'data', 'coordinates') + os.sep
ANCESTRAL = os.path.join(_HERE, '..', 'data', 'ancestral_coordinates') + os.sep
ALTERNATE = os.path.join(_HERE, '..', 'data', 'alternate_coordinates') + os.sep
CACHE = os.path.join(_HERE, '..', 'data', '.cache') + os.sep


def cif(pdb_id):
    for d in (COORD, ANCESTRAL, ALTERNATE):
        local = d + pdb_id + '.cif'
        if os.path.exists(local):
            return local
    os.makedirs(CACHE, exist_ok=True)
    cached = CACHE + pdb_id + '.cif'
    if not os.path.exists(cached):
        url = f'https://files.rcsb.org/download/{pdb_id}.cif'
        print(f'[fetch] {pdb_id} not in repo; downloading {url}')
        urllib.request.urlretrieve(url, cached)
    return cached
