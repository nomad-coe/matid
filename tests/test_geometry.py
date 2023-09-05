import numpy as np
import pytest
import matid.geometry


@pytest.mark.parametrize("positions, cell, pbc, expected", [
    pytest.param(
        np.array([
            [0, 0, 0],
            [0, 4, 0],
        ]),
        np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]]),
        False,
        np.array([
            [[0, 0, 0], [0, -4, 0]],
            [[0, 4, 0], [0, 0, 0]]
        ]),
        id='finite'
    ),
    pytest.param(
        np.array([
            [0, 0, 0],
            [0, 4, 0],
        ]),
        np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]]),
        [True, False, False],
        np.array([
            [[0, 0, 0], [0, -4, 0]],
            [[0, 4, 0], [0, 0, 0]]
        ]),
        id='periodic-a'
    ),
    pytest.param(
        np.array([
            [0, 0, 0],
            [0, 4, 0],
        ]),
        np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]]),
        [False, True, False],
        np.array([
            [[0, 0, 0], [0, 1, 0]],
            [[0, -1, 0], [0, 0, 0]]
        ]),
        id='periodic-b'
    ),
    pytest.param(
        np.array([
            [0, 0, 0],
            [0, 4, 0],
        ]),
        np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]]),
        [False, False, True],
        np.array([
            [[0, 0, 0], [0, -4, 0]],
            [[0, 4, 0], [0, 0, 0]]
        ]),
        id='periodic-c'
    ),
    pytest.param(
        np.array([
            [0, 0, 0],
            [0, 4, 0],
        ]),
        np.array([[5, 0, 0], [0, 5, 0], [0, 0, 5]]),
        True,
        np.array([
            [[0, 0, 0], [0, 1, 0]],
            [[0, -1, 0], [0, 0, 0]]
        ]),
        id='periodic-abc'
    ),
])
def test_displacement_tensor(positions, cell, pbc, expected):
    cutoff = 10
    dist_mat = matid.geometry.get_displacement_tensor(positions, positions, cell, pbc, mic=True, max_distance=cutoff)
    dist_mat_ext = matid.geometry.get_displacement_tensor_ext(positions, cell, pbc, mic=True, cutoff=cutoff)
    assert np.allclose(dist_mat, expected)
    assert np.allclose(dist_mat_ext, expected)
