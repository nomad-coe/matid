import numpy as np
from ase import Atoms
from ase.io import write
from ase.build import fcc111, mx2
from ase.visualize import view

from matid import SBC, SymmetryAnalyzer

# Copper surface
surface = fcc111("Cu", size=(6, 6, 3), vacuum=0)
surface.rotate(v=[0, 0, 1], a=-60, rotate_cell=True)
surface.wrap()

# MoS2
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

# Center structures in the horizontal plane
surface_com = surface.get_center_of_mass()
cell_center = 1 / 2 * surface.get_cell().sum(axis=0)
mos2_com = mos2.get_center_of_mass()
mos2.translate(cell_center - mos2_com)
graphene_com = graphene.get_center_of_mass()
graphene.translate(cell_center - graphene_com)

# Stack structures vertically
surface_top = np.max(surface.get_positions()[:, 2], axis=0)
mos2_bottom = np.min(mos2.get_positions()[:, 2], axis=0)
distance_surface_mos2 = 2
mos2.translate([0, 0, surface_top - mos2_bottom + distance_surface_mos2])
mos2_top = np.max(mos2.get_positions()[:, 2], axis=0)
graphene_bottom = np.max(graphene.get_positions()[:, 2], axis=0)
distance_mos2_graphene = 2
graphene.translate([0, 0, mos2_top - graphene_bottom + distance_mos2_graphene])
system = surface + mos2 + graphene

# Fit in cell and set periodicity
scaled_positions = system.get_scaled_positions()
system.set_cell(3 * system.get_cell())
system.center()
system.set_pbc(False)
write("original.eps", system, rotation="30z,-60x,0y", show_unit_cell=0, maxwidth=5000)

# Peform clustering and view results
sbc = SBC()
clusters = sbc.get_clusters(system)
rotations = [
    "240z,-70x,0y",
    "-30z,-70x,0y",
    "240z,-70x,0y",
]

for i, cluster in enumerate(clusters):
    write(f"cluster{i}.eps", system[cluster.indices], rotation="30z,-60x,0y", show_unit_cell=0, maxwidth=3000)
    unit_cell = cluster.get_cell()
    analyzer = SymmetryAnalyzer(unit_cell)
    conv_system = analyzer.get_conventional_system()
    cell = conv_system.get_cell()
    cell[~conv_system.get_pbc(), :] = 0
    conv_system.set_cell(cell)
    write(f"cluster{i}_unit_cell.eps", conv_system, rotation=rotations[i], show_unit_cell=2, maxwidth=3000)
