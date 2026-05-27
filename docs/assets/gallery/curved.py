"""SBC gallery example: a warped sheet with adsorbates.

Builds a graphene sheet, applies a sinusoidal ripple to curve it, and adsorbs a
few carbon dioxide molecules above it at a physically reasonable distance. Even
though the sheet is curved and decorated with foreign molecules, SBC recovers
the graphene as a single cluster while leaving the adsorbed CO2 unclustered.
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

# Adsorb a few CO2 molecules above the curved sheet. CO2 is linear (O=C=O) with
# a ~1.16 A C-O bond, physisorbed flat ~3.2 A above the *local* sheet height
# (measured under each molecule, so they track the ripple instead of floating).
gpos = graphene.get_positions()


def local_top(x, y, radius=2.0):
    d = np.linalg.norm(gpos[:, :2] - [x, y], axis=1)
    near = gpos[d < radius] if (d < radius).any() else gpos[np.argsort(d)[:3]]
    return near[:, 2].max()


def co2(x, y):
    c = np.array([x, y, local_top(x, y) + 3.2])
    return Atoms("CO2", positions=[c, c + [1.16, 0, 0], c - [1.16, 0, 0]])


com = graphene.get_center_of_mass()
adsorbates = (
    co2(com[0] - 6, com[1] - 4)
    + co2(com[0] + 4, com[1] + 5)
    + co2(com[0] + 9, com[1] - 6)
)
system = graphene + adsorbates
adsorbate_indices = set(range(len(graphene), len(system)))

# Cluster and verify. A looser position tolerance lets the periodic finder
# follow the curvature of the rippled sheet and recover it as a single cluster.
clusters = SBC().get_clusters(system, pos_tol=1.0)
print(f"Found {len(clusters)} clusters:")
for c in clusters:
    symbols = "".join(sorted(set(c.get_atoms().get_chemical_symbols())))
    print(f"  {len(c)} atoms, species={symbols}")
assert len(clusters) == 1, f"expected 1 cluster, got {len(clusters)}"
assert set(clusters[0].get_atoms().get_chemical_symbols()) == {"C"}, "cluster should be pure carbon"
assert not (set(clusters[0].indices) & adsorbate_indices), "CO2 adsorbates should not be clustered"

# Render the original system and the identified cluster
render(system, "original.png")
for i, cluster in enumerate(clusters):
    render(system[cluster.indices], f"cluster{i}.png")
print(f"Wrote images to {os.path.normpath(OUT)}")
