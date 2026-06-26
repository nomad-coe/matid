import numpy as np


class Distances:
    """Container for all distance information that has been extracted from a
    system.

    Supports two backing representations:

    - Dense: the full ``[n, n]`` / ``[n, n, 3]`` matrices. Produced when
      :func:`matid.geometry.get_distances` is called with an infinite cutoff.
    - Sparse: a minimum-image neighbour list (COO/CSR) holding only pairs within
      a finite cutoff. Produced when a finite cutoff is given. This avoids
      allocating the dense ``O(n^2)`` matrices.

    Consumers that only need local distances should use the ``get_*`` accessor
    methods, which work transparently for both representations. The dense matrix
    properties (``disp_tensor_mic`` etc.) are only available in dense mode; in
    sparse mode they raise, since materializing them would defeat the purpose.
    """

    def __init__(
        self,
        disp_tensor_mic,
        disp_factors,
        dist_matrix_mic,
        dist_matrix_radii_mic,
        cell_list=None,
    ):
        # Dense mode (backwards compatible positional constructor).
        self._dense = True
        self._disp_tensor_mic = disp_tensor_mic
        self._disp_factors = disp_factors
        self._dist_matrix_mic = dist_matrix_mic
        self._dist_matrix_radii_mic = dist_matrix_radii_mic
        self.cell_list = cell_list

    @classmethod
    def from_sparse(cls, n, row, col, distance, displacement, factor, radii):
        """Construct a sparse-backed Distances from a COO neighbour list.

        The minimum-image displacement of an entry ``(row, col)`` is
        ``pos_row - pos_col``, matching the dense ``disp_tensor_mic`` convention.
        Both pair directions are expected to be present in the input.

        Args:
            n: Number of atoms.
            row, col: Integer COO indices, shape (nnz,).
            distance: Pairwise distances, shape (nnz,).
            displacement: Displacement vectors, shape (nnz, 3).
            factor: Periodic image factors, shape (nnz, 3).
            radii: Per-atom radii, shape (n,).
        """
        self = cls.__new__(cls)
        self._dense = False
        self._n = n
        self._radii = np.asarray(radii)
        self.cell_list = None

        # Group the COO entries by row to form a CSR-style structure, so that
        # all neighbours of an atom are a contiguous slice.
        order = np.lexsort((col, row))
        self._col = np.asarray(col)[order]
        self._dist = np.asarray(distance)[order]
        self._disp = np.asarray(displacement)[order]
        self._fac = np.asarray(factor)[order]
        sorted_row = np.asarray(row)[order]
        self._row_ptr = np.searchsorted(sorted_row, np.arange(n + 1))
        return self

    # ------------------------------------------------------------------
    # Dense matrix properties (only valid in dense / infinite-cutoff mode).
    # ------------------------------------------------------------------
    def _require_dense(self, name):
        if not self._dense:
            raise AttributeError(
                "'{}' is a dense O(n^2) matrix and is not available for a "
                "sparse (finite-cutoff) Distances. Use the get_* accessor "
                "methods instead.".format(name)
            )

    @property
    def disp_tensor_mic(self):
        self._require_dense("disp_tensor_mic")
        return self._disp_tensor_mic

    @property
    def disp_factors(self):
        self._require_dense("disp_factors")
        return self._disp_factors

    @property
    def dist_matrix_mic(self):
        self._require_dense("dist_matrix_mic")
        return self._dist_matrix_mic

    @property
    def dist_matrix_radii_mic(self):
        self._require_dense("dist_matrix_radii_mic")
        return self._dist_matrix_radii_mic

    # ------------------------------------------------------------------
    # Accessors that work in both dense and sparse mode.
    # ------------------------------------------------------------------
    def _row_slice(self, i):
        return slice(self._row_ptr[i], self._row_ptr[i + 1])

    def get_displacement_column(self, i):
        """Returns ``disp_tensor_mic[:, i]`` of shape (n, 3), i.e. the vectors
        ``pos_a - pos_i`` for every atom ``a``. Entries beyond the cutoff are
        infinite; the ``i``-th entry is zero."""
        if self._dense:
            return self._disp_tensor_mic[:, i]
        out = np.full((self._n, 3), np.inf)
        s = self._row_slice(i)
        # Stored row entries are pos_i - pos_a; the column needs pos_a - pos_i.
        out[self._col[s]] = -self._disp[s]
        out[i] = 0.0
        return out

    def get_factor_column(self, i):
        """Returns ``disp_factors[i, :]`` of shape (n, 3). Entries beyond the
        cutoff are infinite; the ``i``-th entry is zero."""
        if self._dense:
            return self._disp_factors[i]
        out = np.full((self._n, 3), np.inf)
        s = self._row_slice(i)
        out[self._col[s]] = self._fac[s]
        out[i] = 0.0
        return out

    def get_radii_distance_row(self, i):
        """Returns ``dist_matrix_radii_mic[i, :]`` of shape (n,). Entries beyond
        the cutoff are infinite; the ``i``-th entry is ``-2 * radii[i]``."""
        if self._dense:
            return self._dist_matrix_radii_mic[i]
        out = np.full(self._n, np.inf)
        s = self._row_slice(i)
        cols = self._col[s]
        out[cols] = self._dist[s] - (self._radii[i] + self._radii[cols])
        out[i] = -2.0 * self._radii[i]
        return out

    def get_radii_distance_submatrix(self, indices):
        """Returns ``dist_matrix_radii_mic[np.ix_(indices, indices)]``. Entries
        for pairs beyond the cutoff are infinite."""
        indices = np.asarray(indices)
        if self._dense:
            return self._dist_matrix_radii_mic[np.ix_(indices, indices)]
        k = len(indices)
        out = np.full((k, k), np.inf)
        local = {int(g): l for l, g in enumerate(indices)}
        for li, gi in enumerate(indices):
            s = self._row_slice(gi)
            cols = self._col[s]
            dists = self._dist[s]
            radii_gi = self._radii[gi]
            for c, d in zip(cols, dists):
                lj = local.get(int(c))
                if lj is not None:
                    out[li, lj] = d - (radii_gi + self._radii[c])
        di = np.arange(k)
        out[di, di] = -2.0 * self._radii[indices]
        return out
