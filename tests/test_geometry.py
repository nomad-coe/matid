import numpy as np
from pathlib import Path
import pytest
import ase.lattice.hexagonal
import ase.lattice.cubic
from ase.build import molecule, nanotube, bcc100
from ase import Atoms

import matid.geometry

# 0D system
sys0d = molecule("H2O")
sys0d.set_cell([3, 3, 3])
sys0d.center()

# 1D system
sys1d = nanotube(3, 3, length=6, bond=1.4, symbol="Si")
sys1d.set_cell((10, 10, 15))
sys1d.center()

# 2D system
sys2d = Atoms(
    symbols=[6, 6],
    cell=np.array(
        (
            [2.4595121467478055, 0.0, 0.0],
            [-1.2297560733739028, 2.13, 0.0],
            [0.0, 0.0, 10.0],
        )
    ),
    scaled_positions=np.array(
        (
            [0.3333333333333333, 0.6666666666666666, 0.5],
            [0.6666666666666667, 0.33333333333333337, 0.5],
        )
    ),
)
sys2d = sys2d.repeat((3, 3, 1))

# 3D system
sys3d = ase.lattice.cubic.Diamond(
    size=(1, 1, 1), symbol="Si", pbc=True, latticeconstant=5.430710
)

# Surface that has been split by the cell boundary
split = bcc100("Fe", size=(5, 1, 3), vacuum=8)
split.translate([0, 0, 9])
split.set_pbc(True)
split.wrap(pbc=True)


# Surface with a high amplitude wave
wavy = bcc100("Fe", size=(15, 3, 3), vacuum=8)
pos = wavy.get_positions()
x_len = np.linalg.norm(wavy.get_cell()[0, :])
x = pos[:, 0]
z = pos[:, 2]
z_new = z + 3 * np.sin(4 * (x / x_len) * np.pi)
pos_new = np.array(pos)
pos_new[:, 2] = z_new
wavy.set_positions(pos_new)
wavy.set_pbc(True)

# Graphite
graphite = ase.lattice.hexagonal.Graphite(
    size=(1, 1, 1), symbol="C", pbc=True, latticeconstant=(2.461, 6.708)
)

# Molecule with quite little vacuum between cells
small_vacuum = Atoms(
    symbols=["C", "O"],
    positions=[
        [1.13250529, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ],
    cell=[
        [5.29177211, 0.0, 0.0],
        [0.0, 5.29177211, 0.0],
        [0.0, 0.0, 5.29177211],
    ],
    pbc=True,
)
mx2 = ase.build.mx2(
    formula="MoS2", kind="2H", a=3.18, thickness=3.19, size=(5, 5, 1), vacuum=8
)
mx2 = mx2[[0, 12]]
tolerance = 0.2


@pytest.mark.parametrize(
    "positions, cell, pbc, cutoff, expected_disp, expected_dist, expected_factors",
    [
        pytest.param(
            np.array(
                [
                    [0, 0, 0],
                    [0, 4, 0],
                ]
            ),
            np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]]),
            False,
            10,
            np.array([[[0, 0, 0], [0, -4, 0]], [[0, 4, 0], [0, 0, 0]]]),
            np.array([[0, 4], [4, 0]]),
            np.array([[[0, 0, 0], [0, 0, 0]], [[0, 0, 0], [0, 0, 0]]]),
            id="finite",
        ),
        pytest.param(
            np.array(
                [
                    [0, 0, 0],
                    [0, 4, 0],
                ]
            ),
            np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]]),
            [True, False, False],
            10,
            np.array([[[0, 0, 0], [0, -4, 0]], [[0, 4, 0], [0, 0, 0]]]),
            np.array([[0, 4], [4, 0]]),
            np.array([[[0, 0, 0], [0, 0, 0]], [[0, 0, 0], [0, 0, 0]]]),
            id="periodic-a",
        ),
        pytest.param(
            np.array(
                [
                    [0, 0, 0],
                    [0, 4, 0],
                ]
            ),
            np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]]),
            [False, True, False],
            10,
            np.array([[[0, 0, 0], [0, 1, 0]], [[0, -1, 0], [0, 0, 0]]]),
            np.array([[0, 1], [1, 0]]),
            np.array([[[0, 0, 0], [0, -1, 0]], [[0, 1, 0], [0, 0, 0]]]),
            id="periodic-b",
        ),
        pytest.param(
            np.array(
                [
                    [0, 0, 0],
                    [0, 4, 0],
                ]
            ),
            np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]]),
            [False, False, True],
            10,
            np.array([[[0, 0, 0], [0, -4, 0]], [[0, 4, 0], [0, 0, 0]]]),
            np.array([[0, 4], [4, 0]]),
            np.array([[[0, 0, 0], [0, 0, 0]], [[0, 0, 0], [0, 0, 0]]]),
            id="periodic-c",
        ),
        pytest.param(
            np.array(
                [
                    [0, 0, 0],
                    [0, 4, 0],
                ]
            ),
            np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]]),
            True,
            10,
            np.array([[[0, 0, 0], [0, 1, 0]], [[0, -1, 0], [0, 0, 0]]]),
            np.array([[0, 1], [1, 0]]),
            np.array([[[0, 0, 0], [0, -1, 0]], [[0, 1, 0], [0, 0, 0]]]),
            id="periodic-abc",
        ),
        # pytest.param(
        #     np.array([
        #         [0, 0, 0],
        #         [0, 4, 0],
        #     ]),
        #     np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]]),
        #     True,
        #     0.1,
        #     np.array([
        #         [[0, 0, 0], [np.inf, np.inf, np.inf]],
        #         [[np.inf, np.inf, np.inf], [0, 0, 0]]
        #     ]),
        #     np.array([
        #         [0, np.inf],
        #         [np.inf, 0]
        #     ]),
        #     np.array([[[np.inf, np.inf, np.inf], [0, 0, 0]], [[np.inf, np.inf, np.inf], [0, 0, 0]]]),
        #     id='small-cutoff'
        # ),
        pytest.param(
            np.array(
                [
                    [-1, 0, 0],
                    [6, 0, 0],
                ]
            ),
            np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]]),
            True,
            10,
            np.array([[[0, 0, 0], [-2, 0, 0]], [[2, 0, 0], [0, 0, 0]]]),
            np.array([[0, 2], [2, 0]]),
            np.array([[[0, 0, 0], [-1, 0, 0]], [[1, 0, 0], [0, 0, 0]]]),
            id="outside-cell",
        ),
        pytest.param(
            np.array(
                [
                    [-1, 0, 0],
                    [6, 0, 0],
                ]
            ),
            np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]]),
            True,
            None,
            np.array([[[0, 0, 0], [-2, 0, 0]], [[2, 0, 0], [0, 0, 0]]]),
            np.array([[0, 2], [2, 0]]),
            np.array([[[0, 0, 0], [-1, 0, 0]], [[1, 0, 0], [0, 0, 0]]]),
            id="no-cutoff",
        ),
    ],
)
def test_displacement_tensor(
    positions, cell, pbc, cutoff, expected_disp, expected_dist, expected_factors
):
    disp_tensor, factors, dist_mat = matid.geometry.get_displacement_tensor_old(
        positions,
        positions,
        cell,
        pbc,
        mic=True,
        max_distance=cutoff,
        return_distances=True,
        return_factors=True,
    )
    (
        disp_tensor_ext,
        factors_ext,
        dist_mat_ext,
    ) = matid.geometry.get_displacement_tensor(
        positions,
        cell,
        pbc,
        mic=True,
        cutoff=cutoff,
        return_distances=True,
        return_factors=True,
    )
    assert np.allclose(disp_tensor, expected_disp)
    assert np.allclose(disp_tensor_ext, expected_disp)
    assert np.allclose(dist_mat, expected_dist)
    assert np.allclose(dist_mat_ext, expected_dist)
    assert np.allclose(factors, expected_factors)
    assert np.allclose(factors_ext, expected_factors)


@pytest.mark.parametrize(
    "system, cell, pbc, expected_dimensionality",
    [
        pytest.param(sys0d, None, False, 0, id="0D-finite"),
        pytest.param(sys0d, [10, 10, 10], True, 0, id="0D-periodic-abc"),
        pytest.param(sys1d, None, True, 1, id="1D-periodic-abc"),
        pytest.param(sys1d, None, [False, True, True], 1, id="1D-periodic-bc"),
        pytest.param(sys1d, None, [False, False, True], 1, id="1D-periodic-c"),
        pytest.param(sys2d, None, True, 2, id="2D-periodic-abc"),
        pytest.param(sys2d, None, [True, True, False], 2, id="2D-periodic-ab"),
        pytest.param(sys3d, None, True, 3, id="3D-periodic-abc"),
        pytest.param(
            ase.io.read(
                Path(__file__).parent
                / "data/ROJiORHNwL4q0WTvNUy0mW5s2Buuq+PSX9X4dQR2r1cjQ9kBtuC-wI6MO8B.xyz"
            ),
            None,
            None,
            3,
            id="non-orthogonal-cell",
        ),
        pytest.param(split, None, None, 2, id="split-by-cell-boundary"),
        pytest.param(wavy, None, None, 2, id="wavy-surface"),
        pytest.param(wavy, None, None, 2, id="wavy-surface"),
        pytest.param(graphite, None, None, 3, id="graphite"),
        pytest.param(small_vacuum, None, None, 0, id="small-vacuum"),
    ],
)
def test_dimensionality(system, cell, pbc, expected_dimensionality):
    # Test with two radii setting alternatives
    for radii, cluster_threshold in [
        ("covalent", 2.7),
        ("vdw_covalent", 0.1),
    ]:
        if pbc is not None:
            system.set_pbc(pbc)
        if cell is not None:
            system.set_cell(cell)
        dimensionality = matid.geometry.get_dimensionality(
            system,
            radii=radii,
            cluster_threshold=cluster_threshold,
        )
        assert dimensionality == expected_dimensionality


@pytest.mark.parametrize(
    "system, pbc, cutoff, expected_indices",
    [
        pytest.param(mx2, False, 0, [0, 1], id="finite, zero cutoff"),
        pytest.param(mx2, False, 1, [0, 1], id="finite, nonzero cutoff"),
        pytest.param(mx2, True, 0, [0, 1], id="periodic, zero cutoff"),
        pytest.param(mx2, True, 1, np.tile([0, 1], 27), id="periodic, nonzero cutoff"),
    ],
)
def test_extend_system(system, pbc, cutoff, expected_indices):
    """Test that the correct factor is returned when finding matches that
    are in the neighbouring cells.
    """
    system.set_pbc(pbc)
    extended_system = matid.geometry.get_extended_system(system, cutoff)
    assert np.array_equal(extended_system.indices, expected_indices)


@pytest.mark.parametrize(
    "system, cutoff, position, expected_indices, expected_distances",
    [
        pytest.param(mx2, 0.1, [0, 0, 9.595], [0], [0], id="inside cutoff, zero distance"),
        pytest.param(mx2, 0.1, [0, 0, 9.590], [0], [0.005], id="inside cutoff, nonzero distance"),
        pytest.param(mx2, 0.1, [0, 0, 9.595 + 0.11], [], [], id="outside cutoff"),
    ],
)
def test_cell_list_position(system, cutoff, position, expected_indices, expected_distances):
    """Test that the correct factor is returned when finding matches that
    are in the neighbouring cells.
    """
    cell_list = matid.geometry.get_cell_list(system, cutoff)
    result = cell_list.get_neighbours_for_position(position[0], position[1], position[2])
    assert np.array_equal(result.indices, expected_indices)
    assert np.allclose(result.distances, expected_distances, rtol=0, atol=1e-8)


@pytest.mark.parametrize(
    "system, pbc, position, expected_matches, expected_factors",
    [
        pytest.param(mx2, False, [0, 0, 9.595], [0], (0, 0, 0), id="finite, within tolerance"),
        pytest.param(mx2, False, [0, 0, 9.595 + 1.1 * tolerance], [None], (0, 0, 0), id="finite, out of tolerance"),
        pytest.param(mx2, True, [0, 0, 9.595], [0], (0, 0, 0), id="periodic, within tolerance, same cell",),
        pytest.param(mx2, True, (mx2.get_positions()[0] + mx2.get_cell()[0]), [0], (1, 0, 0), id="orthogonal, match in +a"),
        pytest.param(mx2, True, (mx2.get_positions()[0] - mx2.get_cell()[0]), [0], (-1, 0, 0), id="orthogonal, match in -a"),
        pytest.param(mx2, True, (mx2.get_positions()[0] + mx2.get_cell()[1]), [0], (0, 1, 0), id="orthogonal, match in +b"),
        pytest.param(mx2, True, (mx2.get_positions()[0] - mx2.get_cell()[1]), [0], (0, -1, 0), id="orthogonal, match in -b"),
        pytest.param(mx2, True, (mx2.get_positions()[0] + mx2.get_cell()[2]), [0], (0, 0, 1), id="orthogonal, match in +c"),
        pytest.param(mx2, True, (mx2.get_positions()[0] - mx2.get_cell()[2]), [0], (0, 0, -1), id="orthogonal, match in -c"),
    ],
)
def test_matches(system, pbc, position, expected_matches, expected_factors):
    """Test that the correct factor is returned when finding matches that
    are in the neighbouring cells.
    """
    system.set_pbc(pbc)

    # Old python implementation
    matches, _, _, factors = matid.geometry.get_matches(
        system,
        np.array(position)[None, :],
        numbers=[system.get_atomic_numbers()[0]],
        tolerances=np.array([tolerance]),
    )

    # New CPP implementation
    extended_system = matid.geometry.get_extended_system(system, tolerance)
    cell_list = matid.geometry.get_cell_list(
        extended_system.positions,
        extended_system.indices,
        extended_system.factors,
        tolerance
    )
    matches_ext, _, _, factors_ext = matid.geometry.get_matches_ext(
        system,
        cell_list,
        np.array(position)[None, :],
        numbers=[system.get_atomic_numbers()[0]],
        tolerances=np.array([tolerance]),
    )

    # Make sure that the atom is found in the correct copy
    assert tuple(factors[0]) == expected_factors
    assert tuple(factors_ext[0]) == expected_factors

    # Make sure that the correct atom is found
    assert np.array_equal(matches, expected_matches)
    assert np.array_equal(matches_ext, expected_matches)
