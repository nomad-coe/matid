from pathlib import Path

import pytest

import ase.io
from ase import Atoms
from ase.build import bulk

from conftest import surface, stack, rattle, assert_topology
from matid.clustering import SBC, Cluster


surface_fcc = surface(bulk("Cu", "fcc", a=3.6, cubic=True), [1, 0, 0], vacuum=10)
surface_rocksalt = surface(
    bulk("NaCl", "rocksalt", a=5.64, cubic=True), [1, 0, 0], vacuum=10
)
surface_fluorite = surface(bulk("CaF2", "fluorite", a=5.451), [1, 0, 0], vacuum=10)
surface_fluorite_extended = surface_fluorite * [2, 2, 1]
surface_1 = surface(
    Atoms(
        symbols=["O", "C", "C"],
        scaled_positions=[[0, 0, 0], [1 / 3, 0, 0], [2 / 3, 0, 0]],
        cell=[3, 1, 1],
        pbc=True,
    ),
    [0, 0, 1],
    [1, 1, 3],
    vacuum=10,
)
surface_2 = surface(
    Atoms(
        symbols=["O", "N", "N"],
        scaled_positions=[[0, 0, 0], [1 / 3, 0, 0], [2 / 3, 0, 0]],
        cell=[3, 1, 1],
        pbc=True,
    ),
    [0, 0, 1],
    [1, 1, 3],
    vacuum=10,
)
stacked_shared_species = stack(surface_1, surface_2, distance=1)
bulk_one_atom = Atoms(
    symbols=["C"], scaled_positions=[[0, 0, 0]], cell=[2, 2, 2], pbc=True
)
bulk_unwrapped = Atoms(
    symbols=["C"], scaled_positions=[[0, 0, 0]], cell=[2, 2, 2], pbc=True
)
bulk_unwrapped.translate([-10, -10, -10])
broken = Atoms(
    symbols=["O", "O", "O", "O", "O", "O", "O", "O"],
    scaled_positions=[
        [0, 0, 0],
        [0.5, 0, 0],
        [0, 0.5, 0],
        [0.5, 0.5, 0],
        [0, 0.0, 0.5],
        [0.5, 0.0, 0.5],
        [0, 0.5, 0.5],
        [0.5, 0.5, 0.5],
    ],
    cell=[2, 2, 2],
    pbc=True,
)

broken *= [1, 1, 3]
del broken[1:8]
ase.build.add_vacuum(broken, 10)

alternatives = Atoms(
    symbols=["H", "O", "O", "O", "O", "O", "O", "O"],
    scaled_positions=[
        [0, 0, 0],
        [0.5, 0, 0],
        [0, 0.5, 0],
        [0.5, 0.5, 0],
        [0, 0.0, 0.5],
        [0.5, 0.0, 0.5],
        [0, 0.5, 0.5],
        [0.5, 0.5, 0.5],
    ],
    cell=[2, 2, 2],
    pbc=True,
)

alternatives *= [1, 1, 3]
del alternatives[1:8]
ase.build.add_vacuum(alternatives, 10)

sparse = Atoms(symbols=["C"], scaled_positions=[[0, 0, 0]], cell=[4, 4, 4], pbc=True)


@pytest.mark.parametrize(
    "system, clusters_expected, pbc",
    [
        # Test surfaces
        pytest.param(
            surface_fcc,
            [
                Cluster(
                    range(len(surface_fcc)),
                    dimensionality=2,
                )
            ],
            True,
            id="fcc",
        ),
        pytest.param(
            surface_rocksalt,
            [
                Cluster(
                    range(len(surface_rocksalt)),
                    dimensionality=2,
                )
            ],
            True,
            id="rocksalt",
        ),
        pytest.param(
            surface_fluorite,
            [
                Cluster(
                    range(len(surface_fluorite)),
                    dimensionality=2,
                )
            ],
            True,
            id="fluorite surface",
        ),
        # Finite systems
        pytest.param(
            surface_fcc,
            [
                Cluster(
                    range(len(surface_fcc)),
                    dimensionality=0,
                )
            ],
            False,
            id="finite, fcc",
        ),
        pytest.param(
            surface_rocksalt,
            [
                Cluster(
                    range(len(surface_rocksalt)),
                    dimensionality=0,
                )
            ],
            False,
            id="finite, rocksalt",
        ),
        pytest.param(
            surface_fluorite_extended,
            [
                Cluster(
                    range(len(surface_fluorite_extended)),
                    dimensionality=0,
                )
            ],
            False,
            id="finite, fluorite",
        ),
        # Stacked systems
        pytest.param(
            stacked_shared_species,
            [
                Cluster(
                    range(9),
                    dimensionality=2,
                ),
                Cluster(
                    range(9, 18),
                    dimensionality=2,
                ),
            ],
            True,
            id="stacked, shared species",
        ),
        # Bulk
        pytest.param(
            bulk_one_atom,
            [
                Cluster(
                    [0],
                    dimensionality=3,
                )
            ],
            True,
            id="bulk, only one atom in cluster, still a valid cluster",
        ),
        pytest.param(
            bulk_unwrapped,
            [
                Cluster(
                    [0],
                    dimensionality=3,
                )
            ],
            True,
            id="bulk, unwrapped coordinates",
        ),
        # Misc.
        pytest.param(
            sparse, [], True, id="valid region, no clusters due to sparse cell"
        ),
        pytest.param(
            broken,
            [
                Cluster(
                    range(1, 17),
                    dimensionality=2,
                )
            ],
            True,
            id="remove unconnected outliers from region",
        ),
        pytest.param(
            ase.io.read(Path(__file__).parent.parent / "data/system-CVC.extxyz"),
            [
                Cluster(
                    range(123),
                    dimensionality=0,
                ),
                Cluster(
                    range(123, 178),
                    dimensionality=0,
                ),
            ],
            False,
            id="tests that the new CPP implementation does not break things",
        ),
        # pytest.param(
        #     alternatives,
        #     [
        #         Cluster(
        #             range(1, 17),
        #             dimensionality=2,
        #         )
        #     ],
        #     True,
        #     id="cluster that explains the structure better is chosen",
        # ),
    ],
)
@pytest.mark.parametrize("noise", [0, 0.08])
def test_clustering(system, clusters_expected, pbc, noise):
    system = rattle(system, noise)
    system.set_pbc(pbc)
    results = SBC().get_clusters(system)
    assert_topology(results, clusters_expected)


@pytest.mark.parametrize(
    "cell, pbc, error",
    [
        pytest.param(None, False, False, id="three bases missing"),
        pytest.param([0, 0, 25.4], False, False, id="two bases missing"),
        pytest.param([0, 10.8, 25.4], False, False, id="one basis missing"),
        pytest.param([1, 1, 1], False, False, id="too small cell"),
        pytest.param([0, 1, 1], True, True, id="cannot complete with pbc=True"),
    ],
)
def test_completion(cell, pbc, error):
    """Tests that finite systems where the cell does not contain all of the
    atoms are correctly handled.
    """
    system = surface(bulk("Cu", "fcc", a=3.6, cubic=True), [1, 0, 0], vacuum=10)
    system.set_cell(cell)
    system.set_pbc(pbc)
    if error:
        with pytest.raises(ValueError):
            SBC().get_clusters(system)
    else:
        results = SBC().get_clusters(system)
        assert_topology(
            results,
            [
                Cluster(
                    range(len(surface_fcc)),
                    dimensionality=0,
                )
            ],
        )
