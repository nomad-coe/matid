"""SBC gallery example: a warped sheet with impurities.

Builds a graphene sheet, applies a sinusoidal ripple to curve it, and decorates
it with a few gold adatoms. Even though the sheet is curved and contains
foreign atoms, SBC recovers the graphene as a single cluster while leaving the
isolated gold impurities unclustered.
"""
import os

import numpy as np
from ase import Atoms
from ase.io import write

from matid import SBC

OUT = os.path.join(os.path.dirname(__file__), "..", "..", "static", "img", "gallery", "curved")
ROTATION = "-75x,0y,8z"


def render(atoms, name):
    os.makedirs(OUT, exist_ok=True)
    write(os.path.join(OUT, name), atoms, rotation=ROTATION, show_unit_cell=0, scale=22, maxwidth=900)


# Graphene sheet
graphene = Atoms(
    symbols=["C", "C"],
    cell=np.array(
        (
            [2.4595121467478055, 0.0, 0.0],
            [-1.2297560733739028, 2.13, 0.0],
            [0.0, 0.0, 20.0],
        )
    ),
    scaled_positions=np.array(([1 / 3, 2 / 3, 0.5], [2 / 3, 1 / 3, 0.5])),
    pbc=[True, True, False],
) * [12, 12, 1]
graphene.set_pbc(False)

# Apply a gentle sinusoidal ripple along x to curve the sheet
pos = graphene.get_positions()
x = pos[:, 0]
x_span = x.max() - x.min()
pos[:, 2] += 2.0 * np.sin(2 * np.pi * 1.5 * (x - x.min()) / x_span)
graphene.set_positions(pos)

# Add a few gold adatoms above the curved sheet as impurities
com = graphene.get_center_of_mass()
gold = Atoms(
    "Au3",
    positions=[
        [com[0] - 6, com[1] - 4, graphene.get_positions()[:, 2].max() + 2.3],
        [com[0] + 3, com[1] + 5, graphene.get_positions()[:, 2].min() - 2.3],
        [com[0] + 8, com[1] - 6, graphene.get_positions()[:, 2].max() + 2.3],
    ],
)
system = graphene + gold
gold_indices = set(range(len(graphene), len(system)))

# Cluster and verify
clusters = SBC().get_clusters(system)
print(f"Found {len(clusters)} clusters:")
for c in clusters:
    symbols = "".join(sorted(set(c.get_atoms().get_chemical_symbols())))
    print(f"  {len(c)} atoms, species={symbols}")
assert len(clusters) == 1, f"expected 1 cluster, got {len(clusters)}"
assert set(clusters[0].get_atoms().get_chemical_symbols()) == {"C"}, "cluster should be pure carbon"
assert not (set(clusters[0].indices) & gold_indices), "gold impurities should not be clustered"

# Render the original system and the identified cluster
render(system, "original.png")
for i, cluster in enumerate(clusters):
    render(system[cluster.indices], f"cluster{i}.png")
print(f"Wrote images to {os.path.normpath(OUT)}")
