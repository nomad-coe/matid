"""SBC gallery example: a bulk crystal with a point defect.

Builds a sizable diamond-silicon supercell and substitutes a single interior
atom with a large caesium impurity. SBC recovers the surrounding bulk as one 3D
cluster and leaves the lone impurity atom unclustered, since a single atom
carries no translational symmetry of its own.
"""
import os

import numpy as np
from ase.build import bulk
from ase.io import write

from matid import SBC

OUT = os.path.join(os.path.dirname(__file__), "..", "..", "static", "img", "gallery", "impurity")
ROTATION = "16x,12y,0z"


def render(atoms, name):
    os.makedirs(OUT, exist_ok=True)
    # Reduced radii (scaled covalent radii) open up the diamond channels so the
    # interior impurity stays visible while keeping the relative atom sizes.
    write(os.path.join(OUT, name), atoms, rotation=ROTATION, show_unit_cell=0, scale=24, maxwidth=900, radii=0.6)


# A sizable diamond-silicon crystal. Viewed along [001] the open diamond
# channels let the interior impurity remain visible.
system = bulk("Si", "diamond", a=5.43, cubic=True) * (4, 4, 4)

# Substitute the most central atom with a single, large caesium impurity
center = system.get_positions().mean(axis=0)
impurity_index = int(np.argmin(np.linalg.norm(system.get_positions() - center, axis=1)))
system[impurity_index].symbol = "Cs"

# Cluster and verify
clusters = SBC().get_clusters(system)
print(f"Found {len(clusters)} clusters:")
for c in clusters:
    symbols = "".join(sorted(set(c.get_atoms().get_chemical_symbols())))
    print(f"  {len(c)} atoms, species={symbols}, dimensionality={c.get_dimensionality()}")
assert len(clusters) == 1, f"expected 1 cluster, got {len(clusters)}"
assert set(clusters[0].get_atoms().get_chemical_symbols()) == {"Si"}, "cluster should be pure silicon"
assert impurity_index not in set(clusters[0].indices), "impurity should not be clustered"

# Render the original system and the identified bulk cluster
render(system, "original.png")
for i, cluster in enumerate(clusters):
    render(system[cluster.indices], f"cluster{i}.png")
print(f"Wrote images to {os.path.normpath(OUT)}")
