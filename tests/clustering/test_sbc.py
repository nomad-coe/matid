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
def test_clustering_default(system, clusters_expected, pbc, noise):
    """Tests that the clustering performs correctly with default settings."""
    system = rattle(system, noise)
    system.set_pbc(pbc)
    results = SBC().get_clusters(system)
    assert_topology(results, clusters_expected)


@pytest.mark.parametrize(
    "system, clusters_expected",
    [
        pytest.param(
            ase.io.read(
                Path(__file__).parent.parent / "data/system-AlOPbSeSiC_PbSe.xyz"
            ),
            [
                Cluster(
                    range(0, 86),
                    dimensionality=0,
                ),
                Cluster(
                    range(86, 122),
                    dimensionality=0,
                ),
                Cluster(
                    range(122, 230),
                    dimensionality=0,
                ),
            ],
            id="layers share a pattern that complicates the choice of basis vectors 1",
        ),
        pytest.param(
            ase.io.read(Path(__file__).parent.parent / "data/system-BNPbSeBN.xyz"),
            [
                Cluster(
                    range(0, 60),
                    dimensionality=0,
                ),
                Cluster(
                    range(96, 156),
                    dimensionality=0,
                ),
                Cluster(
                    range(60, 96),
                    dimensionality=0,
                ),
            ],
            id="layers share a pattern that complicates the choice of basis vectors 2",
        ),
    ],
)
def test_clustering_coincidental_patterns(system, clusters_expected):
    """Tests that clustering works correctly in cases where the environment
    contains coindidental patterns that bump up the selection of some basis
    directions."""
    results = SBC().get_clusters(system, pos_tol=0.1)
    assert_topology(results, clusters_expected)


fcc_clusters = [Cluster(range(len(surface_fcc)), dimensionality=0)]
single_atom_clusters = []


@pytest.mark.parametrize(
    "system, cell, pbc, expected_clusters",
    [
        pytest.param(surface_fcc, None, False, fcc_clusters, id="three bases missing"),
        pytest.param(
            surface_fcc, [0, 0, 25.4], False, fcc_clusters, id="two bases missing"
        ),
        pytest.param(
            surface_fcc, [0, 10.8, 25.4], False, fcc_clusters, id="one basis missing"
        ),
        pytest.param(surface_fcc, [1, 1, 1], False, fcc_clusters, id="too small cell"),
        pytest.param(
            bulk_one_atom,
            [0, 0, 0],
            False,
            single_atom_clusters,
            id="single atom completion",
        ),
    ],
)
def test_completion(system, cell, pbc, expected_clusters):
    """Tests that finite systems where the cell does not contain all of the
    atoms are correctly handled.
    """
    system.set_cell(cell)
    system.set_pbc(pbc)
    clusters = SBC().get_clusters(system)
    assert_topology(clusters, expected_clusters)


@pytest.mark.parametrize(
    "cell, pbc",
    [
        pytest.param([0, 1, 1], True, id="|a| = 0"),
        pytest.param([1, 0, 1], True, id="|b| = 0"),
        pytest.param([1, 0, 1], True, id="|c| = 0"),
        pytest.param([0, 0, 0], True, id="|a| = |b| = |c| = 0"),
    ],
)
def test_completion_error(cell, pbc):
    """Tests that exception is raise when trying to process systems with
    zero-length basis with pbc=True.
    """
    system = surface_fcc
    system.set_cell(cell)
    system.set_pbc(pbc)
    with pytest.raises(ValueError):
        SBC().get_clusters(system)
