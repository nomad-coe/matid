"""SBC gallery example: a layered heterostructure.

Builds a stack of three different 2D-periodic components — a Cu(111) surface,
a MoS2 monolayer and a graphene sheet — and uses Symmetry-based Clustering to
recover each component as a separate cluster.
"""
import os

import numpy as np
from ase import Atoms
from ase.io import write
from ase.build import fcc111, mx2

from matid import SBC

OUT = os.path.join(os.path.dirname(__file__), "..", "..", "static", "img", "gallery", "stacked")
ROTATION = "30z,-60x,0y"


def render(atoms, name):
    os.makedirs(OUT, exist_ok=True)
    write(os.path.join(OUT, name), atoms, rotation=ROTATION, show_unit_cell=0, scale=25, maxwidth=900)


# Copper (111) surface
surface = fcc111("Cu", size=(6, 6, 3), vacuum=0)
surface.rotate(v=[0, 0, 1], a=-60, rotate_cell=True)
surface.wrap()

# MoS2 monolayer
mos2 = mx2(formula="MoS2", kind="2H", a=3.18, thickness=3.19, size=(5, 5, 1))

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
) * [6, 6, 1]

# Center the components in the horizontal plane
cell_center = 1 / 2 * surface.get_cell().sum(axis=0)
mos2.translate(cell_center - mos2.get_center_of_mass())
graphene.translate(cell_center - graphene.get_center_of_mass())

# Stack the components vertically with small gaps
surface_top = np.max(surface.get_positions()[:, 2])
mos2.translate([0, 0, surface_top - np.min(mos2.get_positions()[:, 2]) + 2])
mos2_top = np.max(mos2.get_positions()[:, 2])
graphene.translate([0, 0, mos2_top - np.max(graphene.get_positions()[:, 2]) + 2])
system = surface + mos2 + graphene

# Fit in a finite cell
system.set_cell(3 * system.get_cell())
system.center()
system.set_pbc(False)

# Cluster and verify
clusters = SBC().get_clusters(system)
found = sorted(
    (frozenset(c.get_atoms().get_chemical_symbols()) for c in clusters), key=sorted
)
expected = sorted(
    [frozenset({"Cu"}), frozenset({"Mo", "S"}), frozenset({"C"})], key=sorted
)
print(f"Found {len(clusters)} clusters:")
for c in clusters:
    symbols = "".join(sorted(set(c.get_atoms().get_chemical_symbols())))
    print(f"  {len(c)} atoms, species={symbols}, dimensionality={c.get_dimensionality()}")
assert len(clusters) == 3, f"expected 3 clusters, got {len(clusters)}"
assert found == expected, f"unexpected cluster species: {found}"

# Render the original system and each identified cluster
render(system, "original.png")
for i, cluster in enumerate(clusters):
    render(system[cluster.indices], f"cluster{i}.png")
print(f"Wrote images to {os.path.normpath(OUT)}")
