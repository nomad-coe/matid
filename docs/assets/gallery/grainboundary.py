"""SBC gallery example: a grain boundary.

Builds a symmetric tilt bicrystal — two grains of the same FCC copper, rotated
with respect to each other about the z axis and joined along a boundary plane.
SBC separates the system into the two grains based on their differing
orientations, something distance- or composition-based clustering cannot do.
"""
import os

import numpy as np
from ase.build import bulk
from ase.io import write

from matid import SBC

OUT = os.path.join(os.path.dirname(__file__), "..", "..", "static", "img", "gallery", "grainboundary")
ROTATION = "0x,0y,0z"

A = 3.6  # Cu lattice constant
W, D, H = 28.0, 20.0, 12.0  # bicrystal box (x, y, z)
TILT = 30.0  # total misorientation angle (degrees)
GAP = 0.8  # small gap at the boundary to avoid unphysical atom overlaps


def render(atoms, name):
    os.makedirs(OUT, exist_ok=True)
    write(os.path.join(OUT, name), atoms, rotation=ROTATION, show_unit_cell=0, scale=26, maxwidth=900)


def carve(atoms, xlo, xhi):
    """Return atoms whose positions fall inside the given box."""
    p = atoms.get_positions()
    mask = (
        (p[:, 0] >= xlo) & (p[:, 0] <= xhi)
        & (np.abs(p[:, 1]) <= D / 2) & (np.abs(p[:, 2]) <= H / 2)
    )
    return atoms[mask]


# A copper crystal large enough to carve a rotated grain from
block = bulk("Cu", "fcc", a=A, cubic=True) * (24, 24, 4)
block.translate(-block.get_center_of_mass())

# Two grains rotated by +/- TILT/2 about z and joined at the x = 0 plane
left = block.copy()
left.rotate(-TILT / 2, "z", center=(0, 0, 0))
grain1 = carve(left, -W / 2, -GAP / 2)

right = block.copy()
right.rotate(TILT / 2, "z", center=(0, 0, 0))
grain2 = carve(right, GAP / 2, W / 2)

system = grain1 + grain2
system.set_pbc(False)
system.center(vacuum=2)

# Cluster and verify
clusters = SBC().get_clusters(system)
print(f"Found {len(clusters)} clusters:")
for c in clusters:
    symbols = "".join(sorted(set(c.get_atoms().get_chemical_symbols())))
    print(f"  {len(c)} atoms, species={symbols}")
assert len(clusters) == 2, f"expected 2 clusters, got {len(clusters)}"
assert all(set(c.get_atoms().get_chemical_symbols()) == {"Cu"} for c in clusters)

# Order the grains left-to-right so cluster0 is grain 1 (left) and cluster1 is grain 2 (right)
clusters = sorted(clusters, key=lambda c: c.get_atoms().get_positions()[:, 0].mean())

# Render the original system and each identified grain
render(system, "original.png")
for i, cluster in enumerate(clusters):
    render(system[cluster.indices], f"cluster{i}.png")
print(f"Wrote images to {os.path.normpath(OUT)}")
