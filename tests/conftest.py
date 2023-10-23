import numpy as np
from numpy.random import default_rng
from ase import Atoms
from ase.build import surface as ase_surface
import ase.build



def create_graphene():
    system = Atoms(
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
    )
    return system


def create_mos2():
    system = ase.build.mx2(
        formula="MoS2", kind="2H", a=3.18, thickness=3.19, size=(1, 1, 1), vacuum=0
    )
    return system


def create_si(cubic=True):
    system = ase.build.bulk(
        "Si",
        crystalstructure="diamond",
        a=5.430710,
        cubic=cubic,
    )
    return system


def create_fe(cubic=True):
    system = ase.build.bulk(
        "Fe",
        crystalstructure="bcc",
        a=2.834,
        cubic=cubic,
    )
    return system


def create_sic(cubic=True):
    system = ase.build.bulk(
        "SiC",
        crystalstructure="zincblende",
        a=4.329,
        cubic=cubic,
    )
    return system


def rattle(atoms, displacement=0.08):
    rng = default_rng(seed=7)
    noise = rng.random((len(atoms), 3)) - 0.5
    lengths = np.linalg.norm(noise, axis=1)
    noise /= np.expand_dims(lengths, 1)
    noise *= displacement
    atoms_copy = atoms.copy()
    atoms_copy.set_positions(atoms_copy.get_positions() + noise)
    return atoms_copy


def surface(conv_cell, indices, layers=[3, 3, 2], vacuum=10):
    surface = ase_surface(conv_cell, indices, layers[2], vacuum=vacuum, periodic=True)
    surface *= [layers[0], layers[1], 1]
    return surface


def stack(a, b, axis=2, distance=3, vacuum=10):
    a_pos = a.get_positions()[:, axis]
    a_max = np.max(a_pos)
    a_min = np.min(a_pos)
    b_pos = b.get_positions()[:, axis]
    b_max = np.max(b_pos)
    b_min = np.min(b_pos)
    a_shift = np.zeros((len(a), 3))
    a_shift[:, axis] += -a_min
    b_shift = np.zeros((len(b), 3))
    b_shift[:, axis] += -b_min + (a_max - a_min) + distance
    a.translate(a_shift)
    b.translate(b_shift)
    stacked = a + b
    cell = a.get_cell()
    axis_new = cell[axis, :]
    axis_norm = np.linalg.norm(axis_new)
    axis_new = axis_new / axis_norm * (a_max - a_min + b_max - b_min + distance)
    cell[axis, :] = axis_new
    stacked.set_cell(cell)
    ase.build.add_vacuum(stacked, vacuum)
    return stacked


def assert_topology(results, expected):
    # Check that correct clusters are found
    assert len(expected) == len(results)
    cluster_map = {tuple(sorted(x.indices)): x for x in results}
    for cluster_expected in expected:
        cluster = cluster_map[tuple(sorted(cluster_expected.indices))]
        assert cluster.dimensionality() == cluster_expected.dimensionality()
        assert cluster.classification() == cluster_expected.classification()


